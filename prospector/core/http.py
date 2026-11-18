"""HTTP helpers built on urllib only (PyMOL's Python may not ship `requests`).

Every call has a timeout, is cached on disk, and turns low-level failures
into one exception type (NetworkError) with a message a human can act on.
"""
from __future__ import annotations

import json
import os
import socket
import ssl
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

from . import cache

USER_AGENT = "Prospector-PyMOL-plugin/0.1 (python-urllib)"
DEFAULT_TIMEOUT = 20  # seconds per request
ATTEMPTS = 3          # tries per request before giving up


class NetworkError(Exception):
    """Raised when a request fails.  ``offline`` is True when the machine
    seems to have no connection at all; ``status`` holds an HTTP code."""

    def __init__(self, message, offline=False, status=None, url=None):
        Exception.__init__(self, message)
        self.offline = offline
        self.status = status
        self.url = url


# --- notes (e.g. "served from cache while offline"), collected per thread ---
_local = threading.local()


def take_notes():
    notes = getattr(_local, "notes", [])
    _local.notes = []
    return notes


def _note(text):
    notes = getattr(_local, "notes", None)
    if notes is None:
        notes = _local.notes = []
    if text not in notes:
        notes.append(text)


def _host(url):
    return urllib.parse.urlsplit(url).netloc


def _open(request, timeout):
    try:
        return urllib.request.urlopen(request, timeout=timeout)
    except urllib.error.URLError as err:
        # Some bundled Pythons (notably on macOS) lack root certificates.
        # If the optional `certifi` package is present, retry with its bundle.
        if isinstance(getattr(err, "reason", None), ssl.SSLCertVerificationError):
            try:
                import certifi  # optional dependency, never required
            except ImportError:
                raise err
            context = ssl.create_default_context(cafile=certifi.where())
            return urllib.request.urlopen(request, timeout=timeout, context=context)
        raise


def _urllib_transport(url, body, headers, timeout):
    """Send one request.  Returns (status, bytes, headers).  Raises NetworkError."""
    request = urllib.request.Request(url, data=body, headers=headers)
    try:
        with _open(request, timeout) as response:
            headers = dict((k.lower(), v) for k, v in response.headers.items())
            return response.status, response.read(), headers
    except urllib.error.HTTPError as err:
        detail = ""
        try:
            detail = err.read(300).decode("utf-8", "replace")
        except Exception:
            pass
        raise NetworkError("%s answered with HTTP %s. %s" % (_host(url), err.code, detail.strip()),
                           status=err.code, url=url)
    except urllib.error.URLError as err:
        reason = getattr(err, "reason", err)
        if isinstance(reason, ssl.SSLError):
            raise NetworkError(
                "Secure connection to %s failed (%s). PyMOL's Python may be missing "
                "root certificates; installing the optional 'certifi' package usually "
                "fixes this." % (_host(url), reason), url=url)
        if isinstance(reason, socket.timeout) or "timed out" in str(reason):
            raise NetworkError("%s did not answer within %ss." % (_host(url), timeout), url=url)
        raise NetworkError("Could not reach %s (%s)." % (_host(url), reason),
                           offline=True, url=url)
    except (socket.timeout, TimeoutError):
        raise NetworkError("%s did not answer within %ss." % (_host(url), timeout), url=url)
    except OSError as err:
        raise NetworkError("Could not reach %s (%s)." % (_host(url), err), offline=True, url=url)


# Tests replace this with a function that replays recorded responses.
transport = _urllib_transport


def request_json(url, payload=None, ttl=cache.DEFAULT_TTL, timeout=DEFAULT_TIMEOUT,
                 header=None, retry=True):
    """GET (or POST when ``payload`` is given) and decode JSON.

    * a fresh cache entry is returned without touching the network;
    * an empty body (e.g. RCSB's HTTP 204 "no results") returns None;
    * if the network fails, an expired cache entry is used when there is one;
    * with ``header="x-total-results"`` the answer is (data, header value).
    """
    key = cache.make_key(url if header is None else url + "#" + header, payload)
    hit = cache.get(key, ttl)
    if hit is not cache.MISS:
        return _unwrap(hit, header)

    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    body = None
    if payload is not None:
        body = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"

    last_error = None
    for attempt in range(1, ATTEMPTS + 1):  # retries for timeouts and server hiccups
        try:
            status, raw, reply_headers = transport(url, body, headers, timeout)
            text = raw.decode("utf-8", "replace").strip() if raw else ""
            data = json.loads(text) if text else None
            if header is not None:
                data = {"body": data, "header": (reply_headers or {}).get(header.lower())}
            cache.put(key, data, url)
            return _unwrap(data, header)
        except NetworkError as err:
            last_error = err
            retryable = (not err.offline) and (err.status is None or err.status >= 500
                                               or err.status == 429)
            if attempt < ATTEMPTS and retryable and retry:
                time.sleep(0.5 * attempt)
                continue
            break
        except ValueError:
            last_error = NetworkError("%s sent a reply that is not valid JSON." % _host(url),
                                      url=url)
            break

    if last_error.status is None or last_error.status >= 500:
        stale = cache.get(key, ttl=None)
        if stale is not cache.MISS:
            days = (cache.age(key) or 0) / cache.DAY
            _note("Network problem (%s) - showing cached data from %.0f day(s) ago."
                  % (_host(url), days))
            return _unwrap(stale, header)
    raise last_error


def _unwrap(data, header):
    if header is None:
        return data
    return (data or {}).get("body"), (data or {}).get("header")


def download(url, path, timeout=30):
    """Download a file to ``path`` (skipped if it already exists)."""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    tmp = path + ".part"
    try:
        for attempt in range(1, ATTEMPTS + 1):
            try:
                with _open(request, timeout) as response, open(tmp, "wb") as handle:
                    while True:
                        chunk = response.read(1 << 16)
                        if not chunk:
                            break
                        handle.write(chunk)
                break
            except urllib.error.HTTPError:
                raise                      # e.g. 404: retrying will not help
            except (urllib.error.URLError, OSError):
                if attempt == ATTEMPTS:
                    raise
                time.sleep(0.5 * attempt)  # a dropped connection: try again
        os.replace(tmp, path)
        return path
    except urllib.error.HTTPError as err:
        raise NetworkError("%s answered with HTTP %s for %s." % (_host(url), err.code, url),
                           status=err.code, url=url)
    except (urllib.error.URLError, OSError) as err:
        reason = getattr(err, "reason", err)
        timed_out = "timed out" in str(reason)
        raise NetworkError("Download from %s failed (%s)." % (_host(url), reason),
                           offline=not timed_out, url=url)
    finally:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass
