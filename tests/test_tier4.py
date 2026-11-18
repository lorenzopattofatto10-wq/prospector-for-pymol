"""Offline tests for Tier 4: teach mode, tour, quiz, surprise."""
import io
import json
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout

from prospector.core import cache, http, rank, store
from tests import fixture_transport
from tests import test_pymol

HAVE_PYMOL = test_pymol.HAVE_PYMOL

if HAVE_PYMOL:
    from pymol import cmd

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "prospector", "core", "data")


class QuizDataTests(unittest.TestCase):
    def test_every_question_has_an_answer(self):
        with open(os.path.join(DATA, "quiz.json"), encoding="utf-8") as handle:
            entries = json.load(handle)
        self.assertGreaterEqual(len(entries), 10)
        self.assertEqual(len(set(e["pdb"] for e in entries)), len(entries))
        for entry in entries:
            self.assertIn(entry["fold"], ("mainly alpha", "mainly beta", "alpha-beta"))
            self.assertEqual(len(entry["hints"]), 3)
            self.assertTrue(entry["accept"] and entry["facts"] and entry["ligand_name"])
            for key in entry["accept"]:
                self.assertEqual(key, key.lower())
            # no hint may contain the answer
            for hint in entry["hints"][:2]:
                for key in entry["accept"]:
                    if len(key) > 3:
                        self.assertNotIn(key, hint.lower(), (entry["pdb"], hint))


@unittest.skipUnless(HAVE_PYMOL, "PyMOL is not installed in this Python")
class Tier4Tests(unittest.TestCase):
    def setUp(self):
        from prospector.pymol_layer import commands, teach, tier4
        self.commands, self.teach, self.tier4 = commands, teach, tier4
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.replaying
        commands.register()
        cmd.delete("all")
        self.mini = os.path.join(self.home, "mini.pdb")
        with open(self.mini, "w") as handle:
            handle.write(test_pymol.PyMOLTests.MINI)
        self.original = commands._download_pdb
        commands._download_pdb = lambda pdb_id: self.mini       # no network in tests

    def tearDown(self):
        self.commands._download_pdb = self.original
        http.transport = http._urllib_transport
        cmd.delete("all")
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def run_line(self, line):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            cmd.do(line, echo=0)
        return buffer.getvalue()

    def test_teach_mode_is_remembered_and_gates_output(self):
        self.assertFalse(self.teach.enabled())
        self.assertIn("teach mode on", self.run_line("prospector teach on"))
        self.assertTrue(store.get_setting("teach"))              # written to disk
        self.assertIn("nothing in PyMOL changed",
                      " ".join(self.run_line("find lysozime").split()))
        cmd.read_pdbstr(test_pymol.PyMOLTests.MINI, "1PPB")
        text = self.run_line("look chains 1PPB")
        self.assertIn("util.cbc 1PPB", text)
        self.assertIn("# colour by chain", text)
        self.assertIn("teach mode off", self.run_line("prospector teach off"))
        self.assertNotIn("util.cbc 1PPB\n", self.run_line("look chains 1PPB"))

    def test_tour_builds_scenes_and_steps(self):
        self.assertIn("No tour is running", self.run_line("next"))
        cmd.read_pdbstr(test_pymol.PyMOLTests.MINI, "1PPB")
        text = self.run_line("tour 1PPB")
        self.assertIn("1. Overall fold", text)
        self.assertIn("Functional sites", text)
        self.assertIn("Ligand", text)
        scenes = [s for s in cmd.get_scene_list() if s.startswith("tour_")]
        self.assertEqual(len(scenes), len(self.tier4._tour["scenes"]))
        self.assertGreaterEqual(len(scenes), 3)
        self.assertIn("[1/%d]" % len(scenes), text)
        self.assertIn("[2/%d]" % len(scenes), self.run_line("next"))
        self.assertIn("[1/%d]" % len(scenes), self.run_line("prev"))
        self.assertIn("first stop", self.run_line("prev"))
        for _ in scenes:
            last = self.run_line("next")
        self.assertIn("last stop", last)
        # a second tour replaces the first one's scenes instead of piling up
        self.run_line("tour 1PPB")
        self.assertEqual(len([s for s in cmd.get_scene_list() if s.startswith("tour_")]),
                         len(scenes))

    def test_tour_by_protein_name_loads_it_first(self):
        grabbed = self.run_line("grab hemoglobin")
        self.assertLess(grabbed.index("LOADED"), grabbed.index("WHAT NEXT"))
        self.assertIn("From now on you can call it 'hemoglobin'", grabbed)
        cmd.delete("all")
        text = self.run_line("look chains hemoglobin")     # a name, nothing loaded
        self.assertIn("Loaded the best match for 'hemoglobin'", text)
        self.assertIn("Applied 'chains'", text)
        cmd.delete("all")
        text = self.run_line("tour hemoglobin")           # nothing is loaded yet
        self.assertIn("Loaded the best match for 'hemoglobin'", text)
        loaded = [n for n in cmd.get_names("objects") if n in ("1A3N", "2HHB", "4HHB")]
        self.assertEqual(len(loaded), 1)
        self.assertIn("Tour of %s" % loaded[0], text)
        # ambiguous names and typos are explained, never guessed
        self.assertIn("could mean several different proteins", self.run_line("tour insulin"))
        self.assertIn("did you mean 'lysozyme'", self.run_line("pocket lysozime"))
        # the name now stands for the loaded object: no second search, no "no object"
        self.assertEqual(self.commands._known_object("Hemoglobin"), loaded[0])
        again = self.run_line("tour hemoglobin")
        self.assertIn("Tour of %s" % loaded[0], again)
        self.assertNotIn("Loaded the best match", again)

    def test_name_used_in_grab_works_in_later_commands(self):
        self.run_line("grab kinase ligand=STI")
        obj = [n for n in cmd.get_names("objects")][0]
        for line in ("about kinase", "tour kinase", "features kinase", "pocket kinase"):
            text = self.run_line(line)
            self.assertNotIn("no object called", text, line)
            self.assertNotIn("I searched for it", text, line)
        self.assertIn("Tour of %s" % obj, self.run_line("tour kinase"))
        self.assertEqual(self.commands._pick_object("KINASE")[0], obj)
        # ... and in the Tier 3 commands
        self.assertIn("Applied 'chains' to %s" % obj, self.run_line("look chains kinase"))
        self.assertIn("Added %s" % obj, self.run_line("fav add kinase my note"))
        cmd.delete(obj)                        # once deleted, the name is free again
        self.assertIsNone(self.commands._known_object("kinase"))

    def test_about_lists_key_places_and_focus_zooms_to_them(self):
        cmd.read_pdbstr(test_pymol.PyMOLTests.MINI, "1PPB")
        about = self.run_line("about 1PPB")
        self.assertIn("KEY PLACES in 1PPB", about)
        self.assertLess(about.index("EXPERIMENT"), about.index("KEY PLACES"))
        self.assertRegex(about, r"1  Active site\s+H/HIS57, H/ASP102, H/ALA195\s+focus 1")
        self.assertRegex(about, r"Ligand 0G6 .* focus 2")
        self.assertIn("Whole structure", about)
        # by number, by name, by residue, with a typo, and back out again
        text = self.run_line("focus 1")
        self.assertIn("NOW SHOWING  active site of 1PPB", text)
        self.assertIn("ON SCREEN    3 residues as red sticks", text)
        self.assertIn("WHAT IT IS   the residues that carry out", text)
        labels = []
        cmd.iterate("1PPB and name CA", "labels.append(label)", space={"labels": labels})
        self.assertEqual(sorted(l for l in labels if l), ["ALA195", "ASP102", "HIS57"])
        self.assertIn("ligand 0G6", self.run_line("focus ligand"))
        self.assertIn("active site", self.run_line("focus activ site"))
        residue = self.run_line("focus H/57")
        self.assertIn("HIS = histidine, residue number 57 of chain H", residue)
        self.assertIn("the whole structure", self.run_line("focus all"))
        # a surface from an earlier preset must not hide the close-up
        self.run_line("look hydrophobicity 1PPB")
        self.assertGreater(cmd.count_atoms("1PPB and rep surface"), 0)
        covered = self.run_line("focus ligand")
        self.assertEqual(cmd.count_atoms("1PPB and rep surface"), 0)
        self.assertIn("the surface was hidden so you can see inside", covered)
        self.assertIn("KEY PLACES in 1PPB", self.run_line("focus"))
        unknown = self.run_line("focus banana")
        self.assertIn("I don't know a place called 'banana'", unknown)
        self.assertIn("focus 1", unknown)                 # ... and shows what exists
        for line in (about + text + residue).splitlines():
            self.assertLessEqual(len(line), 80, line)

    def test_tour_says_what_is_on_screen(self):
        cmd.read_pdbstr(test_pymol.PyMOLTests.MINI, "1PPB")
        first = self.run_line("tour 1PPB")
        for label in ("You see", "Shapes", "Colours", "Facts", "Try"):
            self.assertIn("   %-12s" % label, first)
        self.assertIn("spirals = alpha-helices, flat arrows = beta-strands", first)
        second = " ".join(self.run_line("next").split())
        self.assertIn("each domain painted in its own colour", second)
        second = " ".join(self.run_line("next").split())
        self.assertIn("FUNCTIONAL SITES", second)
        self.assertIn("red sticks = active site (HIS57, ASP102, ALA195)", second)
        self.assertIn("blue = nitrogen, red = oxygen", second)
        third = " ".join(self.run_line("next").split())
        self.assertIn("the bound molecule 0G6 in the centre, as yellow sticks", third)
        self.assertIn("yellow dashed lines = likely hydrogen bonds", third)

    def test_tour_without_annotations(self):
        cmd.fab("ACDEFGHIKLMNPQRSTVWY", "peptide", ss=1)
        text = self.run_line("tour peptide")
        self.assertIn("not named after a PDB ID", text)
        self.assertIn("Tour of peptide: 1 stops", text)
        self.assertIn("mostly alpha-helices", text)

    def test_quiz_round(self):
        self.assertIn("No quiz is running", self.run_line("quiz hint"))
        text = self.run_line("quiz")
        self.assertIn("mystery", cmd.get_names())
        self.assertIn("1. What protein is this?", text)
        entry = self.tier4._quiz["entry"]
        self.assertNotIn(entry["pdb"], text)                     # no spoilers
        self.assertNotIn(entry["name"], text)
        self.assertIsNone(self.commands._identify("mystery"))    # `about` cannot cheat
        self.assertIn("Hint 1 of 3", self.run_line("quiz hint"))
        self.assertIn("Not that one", self.run_line("quiz guess banana"))
        self.assertIn("Correct", self.run_line("quiz guess the %s" % entry["accept"][0]))
        reveal = self.run_line("quiz reveal")
        self.assertIn(entry["name"], reveal)
        self.assertIn(entry["fold_note"], " ".join(reveal.split()))
        self.assertIn(entry["pdb"], cmd.get_names())             # renamed after the reveal
        self.assertNotIn("mystery", cmd.get_names())

    def test_surprise_loads_a_classic_and_survives_missing_card(self):
        text = self.run_line("surprise")       # the info card is not in the fixtures
        loaded = [n for n in cmd.get_names() if n in rank.classics()]
        self.assertEqual(len(loaded), 1)
        self.assertIn("Surprise! %s" % loaded[0], text)
        self.assertIn(rank.classics()[loaded[0]], text)

    def test_new_commands_are_registered(self):
        for name in ("tour", "next", "prev", "quiz", "surprise", "compare", "look"):
            self.assertIn(name, cmd.keyword)


if __name__ == "__main__":
    unittest.main()
