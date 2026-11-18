"""Offline unit tests for the pure-Python core.  Run from the project folder:

    python -m unittest discover -s tests -v
"""
import shutil
import tempfile
import time
import unittest

from prospector.core import cache, filters, fuzzy, http, organisms, rank, report, search, uniprot
from tests import fixture_transport


class CoreCase(unittest.TestCase):
    def setUp(self):
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.replaying

    def tearDown(self):
        http.transport = http._urllib_transport
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)


class FilterTests(unittest.TestCase):
    def test_all_filters(self):
        q = filters.parse("spike protein organism=human method=EM res<2.5 ligand=atp "
                          "year>2015 limit=5")
        self.assertEqual(q.text, "spike protein")
        self.assertEqual(q.organism.taxid, 9606)
        self.assertEqual(q.method, "em")
        self.assertEqual(q.res_max, 2.5)
        self.assertEqual(q.ligands, ["ATP"])
        self.assertEqual(q.year_min, 2016)
        self.assertEqual(q.limit, 5)

    def test_quoted_organism(self):
        q = filters.parse('dnaK organism="e. coli"')
        self.assertEqual((q.text, q.organism.taxid), ("dnaK", 562))

    def test_unknown_filter_suggests(self):
        with self.assertRaises(filters.FilterError) as caught:
            filters.parse("p53 organsim=human")
        self.assertIn("Did you mean 'organism'", str(caught.exception))

    def test_bad_values(self):
        for bad in ("p53 res<high", "p53 method=magic", "p53 limit=0", "p53 ligand=a-b-c-d"):
            with self.assertRaises(filters.FilterError):
                filters.parse(bad)

    def test_filters_round_trip(self):
        q = filters.parse("p53 method=xray res<2 year>2015")
        self.assertEqual(q.filter_words(), ["method=xray", "res<2", "year>2015"])


class FuzzyTests(unittest.TestCase):
    def test_british_spelling(self):
        text, notes = fuzzy.normalize_spelling("haemoglobin")
        self.assertEqual(text, "hemoglobin")
        self.assertTrue(notes)
        self.assertEqual(fuzzy.normalize_spelling("leucine zipper")[0], "leucine zipper")

    def test_typos(self):
        self.assertEqual(fuzzy.suggest("lysozime"), "lysozyme")
        self.assertEqual(fuzzy.suggest("hemoglobbin"), "hemoglobin")
        self.assertIsNone(fuzzy.suggest("asdfqwerty"))
        self.assertIsNone(fuzzy.suggest("lysozyme"))


class OrganismTests(unittest.TestCase):
    def test_lookup_and_extract(self):
        self.assertEqual(organisms.lookup("Human").taxid, 9606)
        self.assertEqual(organisms.lookup("e_coli").taxid, 562)
        self.assertIsNone(organisms.lookup("axolotl").taxid)
        rest, organism = organisms.extract("spike protein sars-cov-2")
        self.assertEqual((rest, organism.taxid), ("spike protein", 2697049))
        self.assertEqual(organisms.extract("horseradish peroxidase")[1], None)
        self.assertEqual(organisms.extract("human")[1], None)   # nothing would be left


class CacheTests(CoreCase):
    def test_expiry(self):
        key = cache.make_key("https://example.org/x")
        cache.put(key, {"a": 1})
        self.assertEqual(cache.get(key, ttl=60), {"a": 1})
        time.sleep(0.05)
        self.assertIs(cache.get(key, ttl=0.01), cache.MISS)
        self.assertEqual(cache.get(key, ttl=None), {"a": 1})
        self.assertEqual(cache.clear(), 1)


class RankTests(unittest.TestCase):
    def test_scores_are_bounded_and_monotonic(self):
        self.assertEqual(rank.resolution_score(0.9), 1.0)
        self.assertGreater(rank.resolution_score(1.5), rank.resolution_score(2.5))
        self.assertEqual(rank.resolution_score(5.0), 0.0)
        self.assertEqual(rank.resolution_score(None), 0.5)
        self.assertGreater(rank.method_score(["X-RAY DIFFRACTION"]),
                           rank.method_score(["SOLUTION NMR"]))
        self.assertEqual(sum(rank.WEIGHTS.values()), 100)

    def test_pdb_ref_coverage_is_a_union(self):
        ref = uniprot.PdbRef("1ABC", "X-ray", 2.0, [(1, 100), (50, 150), (300, 310)])
        self.assertEqual(ref.covered(), 161)


class AcceptanceTests(CoreCase):
    """The Tier 1 acceptance list, replayed from recorded API answers."""

    def test_hemoglobin_classic_human_entry_near_top(self):
        result = search.find("hemoglobin")
        self.assertEqual(result.kind, "results")
        top5 = [hit.pdb_id for hit in result.hits[:5]]
        self.assertTrue({"4HHB", "2HHB", "1A3N"} & set(top5), top5)
        self.assertEqual(result.top.entry.organism_text(result.top.accession), "Homo sapiens")
        text = report.render(result)
        self.assertIn("WHY %s IS #1" % result.top.pdb_id, text)
        # the answer comes first, the details after it
        order = [text.index(part) for part in ("BEST MATCH", "TO LOAD IT",
                                               "HOW I READ YOUR SEARCH", " #  PDB", "WHY ")]
        self.assertEqual(order, sorted(order))
        self.assertIn("TO LOAD IT   grab hemoglobin", text)
        self.assertTrue(all(len(l) <= 80 for l in text.splitlines()),
                        max(text.splitlines(), key=len))
        self.assertIn("Run this:   fetch %s" % result.top.pdb_id, text)

    def test_p53_resolves_through_tp53(self):
        result = search.find("p53")
        lead = result.resolution.group.lead
        self.assertEqual((lead.gene, lead.accession), ("TP53", "P04637"))
        self.assertIn("'p53' is the gene TP53 = UniProt P04637", report.render(result))
        for hit in result.hits:
            self.assertIn("P04637", hit.entry.coverage)

    def test_insulin_is_disambiguated(self):
        result = search.find("insulin")
        self.assertEqual(result.kind, "ambiguous")
        genes = [group.lead.gene for group in result.resolution.choices]
        self.assertIn("INS", genes)
        self.assertIn("INSR", genes)
        self.assertEqual(result.command(), "")          # nothing to run: no guessing
        self.assertIn("find INSR", report.render(result))

    def test_spike_cryo_em(self):
        result = search.find("spike protein sars-cov-2 method=em")
        self.assertEqual(result.kind, "results")
        self.assertEqual(result.resolution.group.lead.accession, "P0DTC2")
        for hit in result.hits:
            self.assertIn("ELECTRON MICROSCOPY", hit.entry.methods)

    def test_typo_gets_a_suggestion(self):
        result = search.find("lysozime")
        self.assertEqual(result.kind, "suggest")
        self.assertEqual(result.message,
                         "No structure found for 'lysozime' - did you mean 'lysozyme'?")
        self.assertEqual(result.command(), "find lysozyme")

    def test_kinase_with_imatinib(self):
        result = search.find("kinase ligand=STI")
        self.assertEqual((result.kind, result.mode), ("results", "text"))
        self.assertTrue(result.hits)
        for hit in result.hits:
            self.assertIn("STI", hit.entry.ligands)

    def test_alphafold_fallback(self):
        result = search.find("MYMK organism=human")
        self.assertEqual(result.kind, "alphafold")
        self.assertEqual(result.model.model_id, "AF-A6NI61-F1")
        self.assertTrue(result.model.cif_url.startswith("https://alphafold.ebi.ac.uk/files/"))
        self.assertEqual(result.command(), "alphafold A6NI61")

    def test_nothing_found(self):
        result = search.find("asdfqwerty")
        self.assertEqual(result.kind, "none")
        self.assertIn("No structure found for 'asdfqwerty'", report.render(result))

    def test_spelling_variant_and_organism_filter(self):
        result = search.find("haemoglobin organism=horse limit=3")
        self.assertEqual(len(result.hits), 3)
        self.assertEqual(result.mode, "protein")
        self.assertTrue(any("hemoglobin" in note for note in result.notes))
        self.assertEqual(result.top.entry.organism_text(result.top.accession), "Equus caballus")

    def test_pdb_id(self):
        result = search.find("4hhb")
        self.assertEqual((result.mode, result.command()), ("id", "fetch 4HHB"))

    def test_filters_that_remove_everything(self):
        result = search.find("p53 res<0.5")
        self.assertEqual(result.kind, "nofilter")
        self.assertIn("none pass your filters (res<0.5)", result.message)

    def test_bad_input_never_raises(self):
        self.assertEqual(search.find("").kind, "error")
        self.assertEqual(search.find("p53 colour=red").kind, "error")


class OfflineTests(CoreCase):
    def test_offline_message(self):
        http.transport = fixture_transport.unplugged
        result = search.find("hemoglobin")
        self.assertEqual(result.kind, "error")
        self.assertTrue(result.offline)
        self.assertIn("You seem to be offline", result.message)

    def test_cached_search_survives_unplugging(self):
        self.assertEqual(search.find("p53").kind, "results")       # fills the cache
        http.transport = fixture_transport.unplugged
        again = search.find("p53")
        self.assertEqual(again.kind, "results")

    def test_expired_cache_is_used_when_offline(self):
        key = cache.make_key("https://example.org/data")
        cache.put(key, {"old": True})
        http.transport = fixture_transport.unplugged
        self.assertEqual(http.request_json("https://example.org/data", ttl=0), {"old": True})
        self.assertTrue(any("cached data" in note for note in http.take_notes()))


if __name__ == "__main__":
    unittest.main()
