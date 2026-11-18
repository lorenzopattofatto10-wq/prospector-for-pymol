"""Tiny on-disk JSON cache with expiry.  Pure standard library.

Files live in ``~/.prospector/cache`` (override the base folder with the
``PROSPECTOR_HOME`` environment variable).
"""
from __future__ import annotations

import hashlib
import json
import os
import time

DAY = 86400
DEFAULT_TTL = 7 * DAY
MISS = object()

_home_override = None


def set_home(path):
    """Point Prospector at a different base folder (used by the tests)."""
    global _home_override
    _home_override = path


def home():
    base = _home_override or os.environ.get("PROSPECTOR_HOME") or os.path.join(
        os.path.expanduser("~"), ".prospector")
    return base


def cache_dir():
    path = os.path.join(home(), "cache")
    os.makedirs(path, exist_ok=True)
    return path


def structures_dir():
    path = os.path.join(home(), "structures")
    os.makedirs(path, exist_ok=True)
    return path


def make_key(url, payload=None):
    text = url if payload is None else url + "\n" + json.dumps(payload, sort_keys=True)
    return hashlib.sha1(text.encode("utf-8")).hexdigest()


def _path(key):
    return os.path.join(cache_dir(), key + ".json")


def get(key, ttl=DEFAULT_TTL):
    """Return cached data, or MISS.  ``ttl=None`` ignores the expiry time
    (used to serve stale data when the network is down)."""
    try:
        with open(_path(key), "r", encoding="utf-8") as handle:
            record = json.load(handle)
    except (OSError, ValueError):
        return MISS
    if ttl is not None and time.time() - record.get("t", 0) > ttl:
        return MISS
    return record.get("data")


def age(key):
    """Age of a cache record in seconds, or None."""
    try:
        with open(_path(key), "r", encoding="utf-8") as handle:
            return time.time() - json.load(handle).get("t", 0)
    except (OSError, ValueError):
        return None


def put(key, data, url=""):
    tmp = _path(key) + ".tmp"
    try:
        with open(tmp, "w", encoding="utf-8") as handle:
            json.dump({"t": time.time(), "url": url, "data": data}, handle)
        os.replace(tmp, _path(key))
    except OSError:
        pass  # a read-only disk must never break a search


def clear():
    """Delete all cached API responses.  Returns the number of files removed."""
    count = 0
    folder = cache_dir()
    for name in os.listdir(folder):
        if name.endswith(".json") or name.endswith(".tmp"):
            try:
                os.remove(os.path.join(folder, name))
                count += 1
            except OSError:
                pass
    return count


def info():
    folder = cache_dir()
    files = [os.path.join(folder, n) for n in os.listdir(folder) if n.endswith(".json")]
    size = sum(os.path.getsize(f) for f in files)
    return {"path": folder, "files": len(files), "bytes": size}
