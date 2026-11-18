"""Re-record the API fixtures from the live services.

    python -m tests.record_fixtures        (run from the project folder)
"""
import os
import shutil
import tempfile

from prospector.core import cache, http, search
from tests import fixture_transport


def main():
    if os.path.isdir(fixture_transport.FIXTURES):
        shutil.rmtree(fixture_transport.FIXTURES)
    cache.set_home(tempfile.mkdtemp(prefix="prospector-record-"))
    http.transport = fixture_transport.recording
    for query in fixture_transport.QUERIES:
        result = search.find(query)
        print("%-45s -> %s" % (query, result.kind))
    for label, call in fixture_transport.tier2_calls() + fixture_transport.tier3_calls():
        call()
        print("%-45s -> recorded" % label)
    files = os.listdir(fixture_transport.FIXTURES)
    size = sum(os.path.getsize(os.path.join(fixture_transport.FIXTURES, f)) for f in files)
    print("%d fixtures, %.1f MB" % (len(files), size / 1e6))


if __name__ == "__main__":
    main()
