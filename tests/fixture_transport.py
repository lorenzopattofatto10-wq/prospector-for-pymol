"""Record real API answers once, replay them in the tests (no network needed)."""
import gzip
import hashlib
import json
import os

from prospector.core import http

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")

# The queries the offline tests replay.  record_fixtures.py runs exactly these.
QUERIES = [
    "hemoglobin",
    "p53",
    "insulin",
    "spike protein sars-cov-2 method=em",
    "lysozime",
    "kinase ligand=STI",
    "MYMK organism=human",
    "asdfqwerty",
    "haemoglobin organism=horse limit=3",
    "4hhb",
    "p53 res<0.5",
]


def tier2_calls():
    """API calls behind about / features / pocket, for record_fixtures.py."""
    from prospector.core import alphafold, card, features, rcsb, sifts
    return [
        ("about 1PPB", lambda: card.entry_card("1PPB")),
        ("about 2LZM", lambda: card.entry_card("2LZM")),
        ("features 1PPB", lambda: features.annotate("1PPB")),
        ("features 4HHB", lambda: features.annotate("4HHB")),
        ("features model P04637", lambda: features.annotate_model("P04637")),
        ("about model P04637", lambda: card.model_card(alphafold.prediction("P04637"))),
        ("detail 1BNA", lambda: rcsb.entry_detail("1BNA")),
        ("detail 1IEP", lambda: rcsb.entry_detail("1IEP")),
        ("sifts 2PTC", lambda: sifts.mapping("2PTC")),
    ]


def tier3_calls():
    from prospector.core import rcsb, similar
    return [
        ("similar 1MBN", lambda: similar.find_similar(
            rcsb.chain_sequences("1MBN")["A"], 20, 6, "1MBN")),
    ]


def _path(url, body):
    digest = hashlib.sha1(url.encode("utf-8") + b"\n" + (body or b"")).hexdigest()
    return os.path.join(FIXTURES, digest + ".json.gz")


def recording(url, body, headers, timeout):
    try:
        status, raw, reply_headers = http._urllib_transport(url, body, headers, timeout)
    except http.NetworkError as error:
        if error.status in (400, 404):       # "not found" is an answer worth keeping
            os.makedirs(FIXTURES, exist_ok=True)
            with gzip.open(_path(url, body), "wt", encoding="utf-8") as handle:
                json.dump({"url": url, "status": error.status, "headers": {},
                           "text": ""}, handle)
        raise
    keep = dict((k, v) for k, v in reply_headers.items() if k == "x-total-results")
    record = {"url": url, "status": status, "headers": keep,
              "text": raw.decode("utf-8", "replace")}
    os.makedirs(FIXTURES, exist_ok=True)
    with gzip.open(_path(url, body), "wt", encoding="utf-8") as handle:
        json.dump(record, handle)
    return status, raw, reply_headers


def replaying(url, body, headers, timeout):
    path = _path(url, body)
    if not os.path.exists(path) and os.environ.get("PROSPECTOR_RECORD_MISSING") == "1":
        # top up the fixtures: record only the answers that are not there yet
        return recording(url, body, headers, timeout)
    if not os.path.exists(path):
        raise http.NetworkError("No recorded answer for %s (the tests never use the "
                                "network)." % url, offline=True, url=url)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        record = json.load(handle)
    if record["status"] >= 400:
        raise http.NetworkError("HTTP %s" % record["status"], status=record["status"], url=url)
    return record["status"], record["text"].encode("utf-8"), record["headers"]


def unplugged(url, body, headers, timeout):
    raise http.NetworkError("Could not reach the network (simulated).", offline=True, url=url)
