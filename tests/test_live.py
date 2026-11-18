"""Live tests against the real services.  Off by default:

    PROSPECTOR_LIVE=1 python -m unittest tests.test_live -v
"""
import os
import shutil
import tempfile
import unittest

from prospector.core import alphafold, cache, http, rcsb, search


@unittest.skipUnless(os.environ.get("PROSPECTOR_LIVE") == "1",
                     "set PROSPECTOR_LIVE=1 to run tests that use the network")
class LiveTests(unittest.TestCase):
    def setUp(self):
        self.home = tempfile.mkdtemp(prefix="prospector-live-")
        cache.set_home(self.home)
        http.transport = http._urllib_transport

    def tearDown(self):
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def test_services_answer_in_the_expected_shape(self):
        self.assertEqual(rcsb.entries(["4HHB"])["4HHB"].resolution, 1.74)
        self.assertTrue(alphafold.prediction("P04637").cif_url.endswith(".cif"))

    def test_acceptance_queries(self):
        expected = {"hemoglobin": "results", "p53": "results", "insulin": "ambiguous",
                    "lysozime": "suggest", "kinase ligand=STI": "results",
                    "asdfqwerty": "none", "MYMK organism=human": "alphafold"}
        for query, kind in expected.items():
            self.assertEqual(search.find(query).kind, kind, query)


if __name__ == "__main__":
    unittest.main()
