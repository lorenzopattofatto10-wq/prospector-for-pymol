"""Offline tests for Tier 3: similar, history/favourites, compare, look, snap, batch."""
import io
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout

from prospector.core import cache, http, rcsb, similar, store
from tests import fixture_transport

try:
    from pymol import cmd
    HAVE_PYMOL = True
except ImportError:
    HAVE_PYMOL = False


class Tier3Case(unittest.TestCase):
    def setUp(self):
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.replaying

    def tearDown(self):
        http.transport = http._urllib_transport
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)


class SimilarTests(Tier3Case):
    def test_myoglobin_relatives(self):
        sequence = rcsb.chain_sequences("1MBN")["A"]
        self.assertEqual(len(sequence), 153)
        result = similar.find_similar(sequence, 20, 6, "1MBN")
        self.assertEqual(len(result.families), 6)
        identities = [family.identity for family in result.families]
        self.assertEqual(identities, sorted(identities, reverse=True))
        self.assertGreaterEqual(identities[0], 0.99)
        organisms = set(family.best[2]["organism"] for family in result.families)
        self.assertGreater(len(organisms), 3)       # relatives, not 6 copies of one protein
        text = similar.render(result, "1MBN chain A")
        self.assertIn("distinct proteins", text)
        self.assertIn("Equus caballus", text)

    def test_short_sequence(self):
        self.assertIn("at least 20", similar.find_similar("ACDEFGHIK").message)


class StoreTests(Tier3Case):
    def test_history_persists_and_clears(self):
        store.add_history("find", "hemoglobin", "results", "1A3N")
        store.add_history("grab", "p53", "results", "1TUP")
        store.add_history("find", "   ", "error")            # blank: ignored
        items = store.history()
        self.assertEqual([i["query"] for i in items], ["hemoglobin", "p53"])
        self.assertEqual(store.history(1)[0]["top"], "1TUP")
        self.assertTrue(os.path.exists(os.path.join(self.home, "history.json")))
        self.assertEqual(store.clear_history(), 2)
        self.assertEqual(store.history(), [])

    def test_history_is_capped(self):
        for number in range(store.HISTORY_LIMIT + 20):
            store.add_history("find", "q%d" % number, "results")
        self.assertEqual(len(store.history()), store.HISTORY_LIMIT)

    def test_favourites(self):
        self.assertTrue(store.add_favorite("4HHB", "Deoxyhaemoglobin", "classic"))
        self.assertFalse(store.add_favorite("4hhb", note="updated"))     # no duplicate
        self.assertEqual([(f["id"], f["note"]) for f in store.favorites()],
                         [("4HHB", "updated")])
        self.assertFalse(store.remove_favorite("1ABC"))
        self.assertTrue(store.remove_favorite("4hhb"))
        self.assertEqual(store.favorites(), [])

    def test_corrupt_file_does_not_crash(self):
        with open(os.path.join(self.home, "favorites.json"), "w") as handle:
            handle.write("{not json")
        self.assertEqual(store.favorites(), [])


@unittest.skipUnless(HAVE_PYMOL, "PyMOL is not installed in this Python")
class PyMOLTier3Tests(Tier3Case):
    def setUp(self):
        Tier3Case.setUp(self)
        from prospector.pymol_layer import commands, tier3
        self.commands, self.tier3 = commands, tier3
        commands.register()
        cmd.delete("all")
        self.cwd = os.getcwd()

    def tearDown(self):
        os.chdir(self.cwd)
        cmd.delete("all")
        Tier3Case.tearDown(self)

    def run_line(self, line):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            cmd.do(line, echo=0)
        return buffer.getvalue()

    def peptides(self):
        sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
        cmd.fab(sequence, "pepA", ss=1)
        cmd.fab(sequence[:16] + "W" + sequence[17:], "pepB", ss=1)     # one mutation
        cmd.translate([15, 3, -7], "pepB")
        cmd.rotate("x", 40, "pepB")

    def test_compare_argument_forms(self):
        split = self.tier3._split_compare
        self.assertEqual(split("1HHO 4HHB")[0], ["1HHO", "4HHB"])
        self.assertEqual(split("human myoglobin , horse myoglobin")[0],
                         ["human myoglobin", "horse myoglobin"])
        self.assertEqual(split("p53 vs MDM2 method=super"),
                         (["p53", "MDM2"], {"method": "super"}))
        self.assertEqual(self.tier3._split_chain("4HHB:A"), ("4HHB", "A"))

    def test_compare_superposes_and_reports(self):
        self.peptides()
        text = self.run_line("compare pepA pepB")
        self.assertIn("Sequence identity: 97%", text)
        self.assertIn("Method: align", text)
        self.assertIn("RMSD: 0.0", text)                 # identical backbones
        self.assertIn("RESULT       pepA and pepB: essentially the same conformation", text)
        order = [text.index(part) for part in ("RESULT", "ON SCREEN", "HOW IT WAS DONE")]
        self.assertEqual(order, sorted(order))        # verdict first, method last
        self.assertLess(cmd.rms_cur("pepB and name CA", "pepA and name CA"), 0.1)
        self.assertIn("cealign", self.run_line("compare pepA pepB method=cealign"))
        self.assertIn("exactly two structures", self.run_line("compare pepA"))
        self.assertIn("method must be", self.run_line("compare pepA pepB method=magic"))

    def test_look_presets_and_unknown_preset(self):
        self.peptides()
        for preset in self.tier3.PRESETS:
            text = self.run_line("look %s pepA" % preset)
            self.assertIn("Applied '%s' to pepA" % preset, text)
        listing = self.run_line("look sparkly")
        self.assertIn("no preset called 'sparkly'", listing)
        self.assertIn("look publication", listing)

    def test_viewas_lists_options_and_changes_the_drawing(self):
        self.peptides()
        menu = self.run_line("viewas")
        for name in self.tier3.STYLE_NAMES:
            self.assertRegex(menu, r"(?m)^   %s\s" % name)
        self.assertIn("viewas <protein> <style>", menu)

        def shown(rep):
            return cmd.count_atoms("pepA and rep %s" % rep)

        text = self.run_line("viewas pepA sticks")              # protein, then style
        self.assertIn("NOW VIEWING  pepA as sticks", text)
        self.assertIn("ON SCREEN every bond as a stick", " ".join(text.split()))
        self.assertTrue(shown("sticks") and not shown("cartoon"))
        self.run_line("viewas surface pepA")                    # style, then protein
        self.assertTrue(shown("surface") and not shown("sticks"))
        self.run_line("viewas pepA ball and stick")
        self.assertTrue(shown("sticks") and shown("spheres"))
        self.run_line("viewas pepA cartoon+lines")              # combination
        self.assertTrue(shown("cartoon") and shown("lines") and not shown("spheres"))
        self.run_line("viewas pepA balls")                      # alias
        self.assertTrue(shown("spheres") and not shown("cartoon"))
        self.assertEqual(shown("spheres"), cmd.count_atoms("pepA"))   # full size again
        self.assertEqual(cmd.count_atoms("pepB and rep spheres"), 0)  # only pepA changed
        typo = self.run_line("viewas pepA stiks")
        self.assertIn("Did you mean 'sticks'?", typo)
        self.assertIn("WAYS TO VIEW A STRUCTURE", typo)
        for line in (menu + text + typo).splitlines():
            self.assertLessEqual(len(line), 80, line)

    def test_snap_never_overwrites(self):
        self.peptides()
        os.chdir(self.home)
        first = self.run_line("snap figure width=300")
        second = self.run_line("snap figure width=300")
        self.assertTrue(os.path.getsize(os.path.join(self.home, "figure.png")) > 0)
        self.assertTrue(os.path.exists(os.path.join(self.home, "figure_2.png")))
        self.assertIn("figure.png", first)
        self.assertIn("figure_2.png", second)
        self.assertIn("width must be", self.run_line("snap figure width=5"))
        cmd.delete("all")
        self.assertIn("Nothing is loaded", self.run_line("snap empty"))

    def test_history_and_fav_commands(self):
        self.run_line("find lysozime")
        self.run_line("find p53 res<0.5")
        text = self.run_line("history")
        self.assertIn("find lysozime -> (suggest)", text)
        self.assertIn("find p53 res<0.5", text)
        self.assertIn("Added 4HHB", self.run_line("fav add 4HHB my favourite"))
        self.assertIn("my favourite", self.run_line("fav"))
        self.assertIn("Removed 4HHB", self.run_line("fav rm 4HHB"))
        self.assertIn("Removed 2 history entries", self.run_line("history clear"))

    def test_batch(self):
        mini = os.path.join(self.home, "mini.pdb")
        with open(mini, "w") as handle:
            handle.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00"
                         "           C\nEND\n")
        original = self.commands._download_pdb
        self.commands._download_pdb = lambda pdb_id: mini      # no network in tests
        listing = os.path.join(self.home, "list.txt")
        with open(listing, "w") as handle:
            handle.write("# my proteins\n4hhb\n\ninsulin\nlysozime\n")
        try:
            text = self.run_line("batch %s" % listing)
        finally:
            self.commands._download_pdb = original
        self.assertIn("4HHB", cmd.get_names())
        self.assertIn("SKIPPED: ambiguous name", text)
        self.assertIn("did you mean 'lysozyme'", text)
        self.assertIn("Batch finished: 1 of 3 loaded", text)
        self.assertIn("cannot find the file", self.run_line("batch nothere.txt"))


if __name__ == "__main__":
    unittest.main()
