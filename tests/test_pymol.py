"""Tests of the PyMOL layer (skipped when PyMOL is not installed)."""
import shutil
import tempfile
import unittest

try:
    from pymol import cmd
    HAVE_PYMOL = True
except ImportError:
    HAVE_PYMOL = False

from prospector.core import cache, http
from tests import fixture_transport


@unittest.skipUnless(HAVE_PYMOL, "PyMOL is not installed in this Python")
class PyMOLTests(unittest.TestCase):
    def setUp(self):
        from prospector.pymol_layer import commands
        self.commands = commands
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.replaying
        commands.register()
        self.results = []
        self.listener = lambda name, result: self.results.append((name, result))
        commands._listeners.append(self.listener)

    def tearDown(self):
        self.commands._listeners.remove(self.listener)
        http.transport = http._urllib_transport
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def test_commands_are_registered(self):
        for name in ("find", "grab", "alphafold", "prospector", "p_find", "p_grab", "info"):
            self.assertIn(name, cmd.keyword)

    def test_info_lists_every_command(self):
        text = self._capture("info")
        for name in ("find", "grab", "alphafold", "about", "features", "pocket", "compare",
                     "similar", "batch", "look", "snap", "history", "fav", "tour", "next",
                     "focus", "viewas",
                     "prev", "quiz", "surprise", "prospector teach", "info"):
            self.assertRegex(text, r"(?m)^   .*\b%s\b" % name)
        self.assertIn(" GET A STRUCTURE", self._capture("prospector help"))
        self.assertLess(len(text.splitlines()), 70)        # short enough to read
        for line in text.splitlines():
            self.assertLessEqual(len(line), 80, line)
        detail = self._capture("info scout")
        self.assertIn("COMMAND  scout", detail)
        self.assertIn("HOW TO TYPE IT", detail)
        self.assertIn("scout p53 histidine", detail)
        self.assertIn("organism=human", self._capture("info filters"))
        self.assertIn("Did you mean: scout?", self._capture("info scuot"))

    def test_no_output_line_is_wider_than_the_console(self):
        self._mini()
        for line in ("find hemoglobin", "find insulin", "find kinase ligand=STI", "info",
                     "about", "features 1PPB", "pocket 1PPB", "find p53 colour=red"):
            for row in self._capture(line).splitlines():
                self.assertLessEqual(len(row), 80, "%s -> %r" % (line, row))

    def test_get_is_left_alone(self):
        self.assertIs(cmd.keyword["get"][0], cmd.get)

    def test_arguments_survive_pymols_parser(self):
        cmd.do("find kinase ligand=STI")
        cmd.do("find p53 res<0.5")
        cmd.do("find organism=horse haemoglobin, limit=3")
        raws = [result.raw for _, result in self.results]
        self.assertEqual(raws, ["kinase ligand=STI", "p53 res<0.5",
                                "organism=horse haemoglobin limit=3"])
        self.assertEqual([r.kind for _, r in self.results], ["results", "nofilter", "results"])

    def test_ambiguous_grab_loads_nothing(self):
        before = cmd.get_names()
        cmd.do("grab insulin")
        self.assertEqual(self.results[-1][1].kind, "ambiguous")
        self.assertEqual(cmd.get_names(), before)

    # ---- Tier 2 (uses a tiny hand-made "structure", so no download is needed)
    MINI = """\
ATOM      1  N   HIS H  57      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  HIS H  57      11.400  10.000  10.000  1.00 20.00           C
ATOM      3  C   HIS H  57      12.000  11.400  10.000  1.00 20.00           C
ATOM      4  O   HIS H  57      11.300  12.400  10.000  1.00 20.00           O
ATOM      5  N   ASP H 102      20.000  10.000  10.000  1.00 20.00           N
ATOM      6  CA  ASP H 102      21.400  10.000  10.000  1.00 20.00           C
ATOM      7  C   ASP H 102      22.000  11.400  10.000  1.00 20.00           C
ATOM      8  O   ASP H 102      21.300  12.400  10.000  1.00 20.00           O
ATOM      9  N   ALA H 195      30.000  10.000  10.000  1.00 20.00           N
ATOM     10  CA  ALA H 195      31.400  10.000  10.000  1.00 20.00           C
ATOM     11  C   ALA H 195      32.000  11.400  10.000  1.00 20.00           C
ATOM     12  O   ALA H 195      31.300  12.400  10.000  1.00 20.00           O
HETATM   13  N1  0G6 H   1      11.300  15.200  10.000  1.00 20.00           N
HETATM   14  C1  0G6 H   1      12.500  16.000  10.000  1.00 20.00           C
HETATM   15 NA    NA H 301      50.000  50.000  50.000  1.00 20.00          NA
END
"""

    def _mini(self):
        cmd.delete("all")
        cmd.read_pdbstr(self.MINI, "1PPB")

    def _capture(self, line):
        import io
        from contextlib import redirect_stdout
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            cmd.do(line, echo=0)
        return buffer.getvalue()

    def test_features_creates_selections_and_checks_residues(self):
        self._mini()
        text = self._capture("features 1PPB types=active")
        self.assertIn("1PPB_active", cmd.get_names("selections"))
        atoms = cmd.get_model("1PPB_active and name CA").atom
        self.assertEqual([(a.chain, a.resi) for a in atoms],
                         [("H", "57"), ("H", "102"), ("H", "195")])
        order = [text.index(part) for part in ("HIGHLIGHTED ON 1PPB", "CHECK", "DETAILS")]
        self.assertEqual(order, sorted(order))        # summary first, details after
        # the mini model has ALA at 195 where UniProt expects Ser: must be reported
        self.assertIn("2 of 3 mapped residues", text)
        self.assertIn("H/195 is ALA (UniProt 568: S)", " ".join(text.split()))

    def test_features_types_with_commas(self):
        self._mini()
        text = self._capture("features 1PPB, types=active,disulfide")
        self.assertIn("Active site", text)
        self.assertIn("types=active,disulfide", text)
        self.assertNotIn("Unknown feature type", text)

    def test_pocket(self):
        self._mini()
        text = self._capture("pocket 1PPB")
        self.assertIn("the pocket around 0G6 in chain H", text)   # not the Na+ ion
        order = [text.index(part) for part in ("NOW SHOWING", "POCKET RESIDUES",
                                               "POSSIBLE HYDROGEN BONDS", "ALL LIGANDS IN")]
        self.assertEqual(order, sorted(order))        # summary first, details after
        self.assertIn("HIS57", text)
        self.assertIn("0G6 N1   - HIS57 O", text)
        self.assertEqual(cmd.count_atoms("1PPB_pocket and name CA"), 1)
        self.assertNotIn("select 1PPB_pocket, byres", text)       # teach mode is off
        from prospector.pymol_layer import teach
        teach.set_enabled(True)
        try:
            taught = self._capture("pocket 1PPB")
        finally:
            teach.set_enabled(False)
        self.assertIn("Teach mode - the PyMOL commands behind this", taught)
        self.assertIn("select 1PPB_pocket, byres", taught)
        self.assertIn("# mode=2 draws polar contacts", taught)
        self.assertIn("no ligand 'XYZ'", self._capture("pocket 1PPB XYZ"))

    def test_about_and_unknown_objects(self):
        self._mini()
        self.assertIn("Prothrombin", self._capture("about"))
        cmd.set_name("1PPB", "mystery_object")
        self.assertIn("don't know which database entry", self._capture("about"))
        # a word that is not an object is treated as a protein name and searched
        self.assertIn("'nothere' is not a loaded object, so I searched for it",
                      self._capture("features nothere"))
        cmd.delete("all")
        self.assertIn("Nothing is loaded yet", self._capture("pocket"))

    def test_offline_grab_does_not_crash(self):
        http.transport = fixture_transport.unplugged
        cmd.do("grab hemoglobin")
        self.assertTrue(self.results[-1][1].offline)


if __name__ == "__main__":
    unittest.main()


@unittest.skipUnless(HAVE_PYMOL and __import__("os").environ.get("PROSPECTOR_QT_TEST") == "1",
                     "set PROSPECTOR_QT_TEST=1 (and QT_QPA_PLATFORM=offscreen on a "
                     "machine without a display) to test the background-thread mode")
class BackgroundModeTests(unittest.TestCase):
    """With a Qt event loop, commands run in a worker thread - but in order."""

    def test_commands_wait_their_turn(self):
        import time
        from pymol.Qt import QtWidgets
        from prospector.pymol_layer import commands, worker
        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(home)
        http.transport = fixture_transport.replaying
        commands.register()
        worker.set_async(True)
        cmd.delete("all")
        events = []
        original = commands._header
        commands._header = lambda text: events.append(text)
        try:
            # a slow command followed immediately by one that needs its result
            worker.submit(lambda: time.sleep(0.3),
                          lambda r, e: cmd.read_pdbstr(PyMOLTests.MINI, "1PPB"))
            cmd.do("features 1PPB types=active", echo=0)
            self.assertGreater(worker.pending(), 0)          # returned at once
            self.assertNotIn("1PPB", cmd.get_names())        # ...before the "download"
            deadline = time.time() + 20
            while worker.pending() and time.time() < deadline:
                app.processEvents()
                time.sleep(0.01)
            self.assertIn("1PPB_active", cmd.get_names("selections"))
        finally:
            commands._header = original
            http.transport = http._urllib_transport
            cache.set_home(None)
            shutil.rmtree(home, ignore_errors=True)
