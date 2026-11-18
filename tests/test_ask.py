"""Offline tests for the plain-language front door: ask, start, status, glossary."""
import io
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout

from prospector.core import cache, http
from tests import fixture_transport, test_pymol

HAVE_PYMOL = test_pymol.HAVE_PYMOL
if HAVE_PYMOL:
    from pymol import cmd


@unittest.skipUnless(HAVE_PYMOL, "PyMOL is not installed in this Python")
class AskTests(unittest.TestCase):
    def setUp(self):
        from prospector.pymol_layer import ask, commands
        self.ask, self.commands = ask, commands
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        os.environ.pop("ANTHROPIC_API_KEY", None)      # never call the AI here
        http.transport = fixture_transport.replaying
        commands.register()
        cmd.delete("all")
        del commands._undo[:]
        cmd.fab("MKHACDEFGHIKLMNPQRSTVWYHCAAC", "pep", ss=1)
        commands._remember("mypep", "pep")

    def tearDown(self):
        http.transport = http._urllib_transport
        cmd.delete("all")
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def run_line(self, line):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            cmd.do(line, echo=0)
        return buffer.getvalue()

    def test_sentences_become_the_right_commands(self):
        expected = {
            "load haemoglobin": "grab haemoglobin",
            "what is mypep": "about mypep",
            "show the active site": "focus active site",
            "what is the active site of this protein": "focus active site",
            "what is the active site": "focus active site",
            "what is the active site of hemoglobin": "focus active site of hemoglobin",
            "where is the active site in mypep": "focus active site of mypep",
            "which residues are the binding sites of p53": "focus binding sites of p53",
            "what is an active site": "glossary active site",
            "what does active site mean": "glossary active site",
            "zoom to the ligand": "focus ligand",
            "where does the drug bind in mypep": "pocket mypep",
            "show it as sticks": "viewas sticks",
            "show mypep as ball and stick": "viewas mypep ballstick",
            "colour by chain": "look chains",
            "show hydrophobicity": "look hydrophobicity",
            "highlight the histidines": "scout histidines",
            "where are the disulfides in mypep": "scout mypep disulfides",
            "show salt bridges": "scout salt bridges",
            "find motif CxxC in mypep": "scout mypep motif CxxC",
            "change arginine 248 to glutamine": "whatif R248Q",
            "change arg 248 to gln in chain A": "whatif A/R248Q",
            "what if K2E in mypep": "whatif mypep K2E",
            "how far is residue 3 from residue 10": "measure 3 10",
            "explain residue A/248": "explain A/248",
            "zoom to residue 5": "focus 5",
            "show chain B": "focus chain B",
            "give me a tour": "tour",
            "what does RMSD mean": "glossary rmsd",
            "what is resolution": "glossary resolution",
            "undo": "back",
            "zoom out": "focus all",
            "what can you do": "info",
            "where am i": "status",
            "i'm new": "start",
            "search for insulin": "find insulin",
            "sequence of mypep": "sequence mypep",
            "quiz me": "quiz",
            "show the pymol commands": "prospector teach on",
            "mypep": "about mypep",
            "insulin": "find insulin",
        }
        for sentence, command in expected.items():
            got = self.ask._interpret(sentence)[0]
            self.assertEqual(" ".join((got or "").split()), command, sentence)
        compare = self.ask._interpret("compare 1HHO and 4HHB")[0]
        self.assertEqual(compare.lower(), "compare 1hho, 4hhb")
        self.assertIsNone(self.ask._interpret("make me a sandwich")[0])

    def test_several_requests_in_one_sentence(self):
        split = self.ask._split_steps
        self.assertEqual(split("load p53 and show it as sticks"),
                         ["load p53", "show it as sticks"])
        self.assertEqual(split("compare 1HHO and 4HHB"), ["compare 1HHO and 4HHB"])
        self.assertEqual(len(split("show it as ball and stick")), 1)
        self.assertEqual(len(split("load p53 then colour by chain then zoom out")), 3)

    def test_ask_shows_the_command_and_runs_it(self):
        text = self.run_line("ask show mypep as sticks")
        self.assertIn("I UNDERSTOOD  draw mypep as sticks", text)
        self.assertIn("RUNNING     viewas mypep sticks", text)
        self.assertIn("NOW VIEWING  pep as sticks", text)        # it really ran
        self.assertTrue(cmd.count_atoms("pep and rep sticks"))
        both = self.run_line("ask show it as spheres and highlight the histidines")
        self.assertIn("1. I UNDERSTOOD", both)
        self.assertIn("2. I UNDERSTOOD", both)
        self.assertIn("SPOTTED      3 histidine residues in pep", both)

    def test_how_do_i_only_explains(self):
        cmd.show_as("cartoon", "pep")
        text = self.run_line("ask how do I show it as sticks")
        self.assertIn("TYPE        viewas sticks", text)
        self.assertIn("WHAT IT DOES", text)
        self.assertIn("Nothing was run", text)
        self.assertFalse(cmd.count_atoms("pep and rep sticks"))   # really nothing ran

    def test_unclear_requests_do_nothing_and_offer_choices(self):
        cmd.show_as("cartoon", "pep")
        text = self.run_line("ask make me a sandwich")
        self.assertIn("I'm not sure what you mean by 'make me a sandwich', so I did "
                      "nothing", " ".join(text.split()))
        self.assertIn("THINGS I UNDERSTAND", text)
        # if one part of a sentence is unclear, none of it runs
        mixed = self.run_line("ask show it as sticks and then frobnicate the widget")
        self.assertIn("so I did nothing", " ".join(mixed.split()))
        self.assertFalse(cmd.count_atoms("pep and rep sticks"))
        self.assertIn("ask load haemoglobin", self.run_line("ask"))

    def test_ask_runs_before_lines_typed_after_it(self):
        """In background mode, what `ask` starts must not be overtaken."""
        from prospector.pymol_layer import worker
        order = []
        with worker.next_in_line():
            worker._jobs.put("first")
            worker._jobs.put("second")
        worker._jobs.put("typed later")
        worker._jobs._items.insert(0, None) if False else None
        self.assertEqual(worker._jobs._items[-3:], ["first", "second", "typed later"])
        del worker._jobs._items[-3:]
        # something already waiting is overtaken by what ask starts
        worker._jobs.put("waiting")
        with worker.next_in_line():
            worker._jobs.put("from ask")
        self.assertEqual(worker._jobs._items[-2:], ["from ask", "waiting"])
        del worker._jobs._items[-2:]

    def test_show_me_something_on_screen_never_loads_another_protein(self):
        cmd.delete("all")
        cmd.read_pdbstr(test_pymol.PyMOLTests.MINI, "1PPB")
        self.commands._last[0] = "1PPB"
        self.assertEqual(self.ask._interpret("show me the 0G6")[0], "focus 0G6")
        command, reason = self.ask._interpret("show me the heme")
        self.assertIsNone(command)                     # no heme here: do not guess
        self.assertIn("grab heme", reason)
        text = " ".join(self.run_line("ask show me the heme").split())
        self.assertIn("so I did nothing", text)
        self.assertIn("To load a protein called that: grab heme", text)
        self.assertEqual(cmd.get_names("objects"), ["1PPB"])
        # an explicit "load" is still a load
        self.assertEqual(self.ask._interpret("load heme")[0], "grab heme")

    def test_asking_about_the_active_site_highlights_it(self):
        cmd.delete("all")
        cmd.read_pdbstr(test_pymol.PyMOLTests.MINI, "1PPB")
        self.commands._last[0] = "1PPB"
        text = self.run_line("ask what is the active site of this protein")
        self.assertIn("RUNNING     focus active site", text)
        self.assertIn("NOW SHOWING  active site of 1PPB", text)
        self.assertIn("H/HIS57, H/ASP102, H/ALA195", text)
        self.assertTrue(cmd.count_atoms("1PPB and resi 57 and rep sticks"))
        labels = []
        cmd.iterate("1PPB and name CA", "labels.append(label)", space={"labels": labels})
        self.assertIn("HIS57", labels)
        # the same, naming the protein
        named = self.run_line("focus active site of 1PPB")
        self.assertIn("NOW SHOWING  active site of 1PPB", named)

    def test_other_wordings_plurals_and_typos(self):
        """The same request said differently must give the same command."""
        self.commands._remember("p53", "pep")
        same = {
            "focus active site of p53": [
                "what are the p53 active sites", "p53 active sites",
                "active sites of p53 please", "where are the catalytic residues of p53",
                "point out the functional residues in p53"],
            "focus binding sites of p53": ["show me the p53 binding sites"],
            "focus variants of p53": [
                "show the disease mutations of p53", "mutation hotspots in p53"],
            "focus active site": [
                "which residues make up the active site", "get closer to the active site"],
            "scout histidines": ["highlight all the histidnes"],
            "scout disulfides": ["where are the sulfur bridges"],
            "scout salt bridges": ["show ionic bonds"],
            "scout hydrophobic": ["show me the oily residues"],
            "scout exposed": ["which residues are on the outside"],
            "scout buried": ["show the hidden residues"],
            "viewas sticks": ["make it look like sticks"],
            "look hydrophobicity": ["show hydrofobicity"],
            "compare 1hho, 4hhb": ["how different are 1HHO and 4HHB", "1HHO vs 4HHB"],
            "measure 248 273": ["how far apart are residues 248 and 273"],
            "whatif R248Q": ["what happens if arginine 248 becomes glutamine"],
            "glossary rmsd": ["what does rmsd stand for"],
            "glossary ligand": ["whats a ligand"],
            "similar p53": ["proteins that resemble p53"],
            "grab myoglobin": ["i want to see myoglobin"],
        }
        for command, sentences in same.items():
            for sentence in sentences:
                self.assertEqual(self.ask._interpret(sentence)[0], command, sentence)
        for nonsense in ("make me a sandwich", "do my homework"):
            self.assertIsNone(self.ask._interpret(nonsense)[0], nonsense)

    def test_two_hundred_ways_of_asking(self):
        """tests/ask_phrases.txt: every sentence must give the command next to it."""
        import os
        self.commands._remember("p53", "pep")
        path = os.path.join(os.path.dirname(__file__), "ask_phrases.txt")
        wrong = []
        for line in open(path, encoding="utf-8"):
            if "=>" not in line or line.startswith("#"):
                continue
            sentence, expected = [part.strip() for part in line.split("=>", 1)]
            got = " ; ".join(str(self.ask._interpret(step)[0]).strip()
                             for step in (self.ask._split_steps(sentence) or [sentence]))
            if got != expected:
                wrong.append("%s: expected %r, got %r" % (sentence, expected, got))
        self.assertEqual(wrong, [])

    def test_a_change_without_a_target_is_completed(self):
        tools = self.ask.tools
        words, notes = tools._fill_in_changes("pep", ["first"], None)
        self.assertEqual(words, ["M1A"])          # methionine 1 -> alanine
        self.assertTrue(any("alanine" in note for note in notes))
        self.assertEqual(tools._fill_in_changes("pep", ["3W"], None)[0], ["H3W"])
        self.assertEqual(tools._fill_in_changes("pep", ["last"], None)[0], ["C28A"])
        self.assertEqual(tools._fill_in_changes("pep", ["4"], None)[0], ["A4G"])
        # two halves of one change are left for the "did you mean" hint
        self.assertEqual(tools._fill_in_changes("pep", ["1V", "1A"], None)[0], ["1V", "1A"])
        out = self.run_line("whatif mypep first")
        self.assertIn("M1A", out)

    def test_clear_empties_the_window_and_back_restores_it(self):
        from pymol import cmd
        out = self.run_line("clear")
        self.assertIn("CLEARED", out)
        self.assertEqual(cmd.get_names("all"), [])
        self.assertIn("already empty", self.run_line("clear"))
        out = self.run_line("back")
        self.assertIn("WENT BACK", out)
        self.assertIn("pep", cmd.get_names("objects"))
        self.assertEqual(self.commands._known_object("mypep"), "pep")
        self.assertIn("clear", self.run_line("info"))

    def test_the_rules_step_aside_when_a_word_is_left_over(self):
        refused = [
            "show me p53 normal and p53 most common mutated version in breast cancer",
            "zooma sul sito attivo della proteina e evidenzia le cariche negative",
            "show it as a surface colored by charge", "show electrostatics",
            "dock the ligand", "blast this sequence", "which residues are conserved",
            "colour chain B green and chain A red", "show side chains of the active site",
            "show surface with 50% transparency", "make me a sandwich"]
        self.commands._remember("p53", "pep")
        for sentence in refused:
            self.assertIsNone(self.ask._interpret(sentence)[0], sentence)
        # ... and why: the word nobody accounted for
        self.assertEqual(self.ask._unaccounted("viewas surface",
                                               "show it as a surface colored by charge"),
                         ["colored", "charges"])
        self.assertEqual(self.ask._unaccounted("whatif p53 R248Q",
                                               "what if arginine 248 becomes glutamine "
                                               "in p53"), [])
        self.assertTrue(self.ask._name_like("green fluorescent", "x", False))
        self.assertTrue(self.ask._name_like("atp", "show ATP", False))
        self.assertFalse(self.ask._name_like("side chains", "x", False))
        self.assertFalse(self.ask._name_like("50 transparency", "x", False))

    def test_new_requests_really_run(self):
        from pymol import cmd
        self.run_line("ask white background")
        self.run_line("ask hide the sticks")
        out = self.run_line("ask make the sticks thicker")
        self.assertIn("set stick_radius, 0.4", out)
        self.run_line("viewas mypep sticks")
        self.run_line("viewas mypep spheres")
        self.run_line("look rainbow mypep")
        out = self.run_line("ask undo the last two steps")
        self.assertIn("RUNNING     back 2", out)
        self.assertEqual(out.count("WENT BACK"), 2)
        out = self.run_line("ask show residues 3 to 6")
        self.assertIn("select chosen, pep and (resi 3-6)", out)
        self.assertEqual(cmd.count_atoms("chosen and name CA"), 4)
        out = self.run_line("ask delete mypep")
        self.assertIn("RUNNING     delete pep", out)
        self.assertNotIn("pep", cmd.get_names("objects"))

    def test_glossary(self):
        self.assertIn("root-mean-square deviation", self.run_line("glossary RMSD"))
        self.assertIn("B-FACTOR", self.run_line("glossary b factor"))        # alias
        self.assertIn("imidazole", self.run_line("glossary histidine"))
        listing = self.run_line("glossary")
        self.assertIn("resolution", listing)
        self.assertIn("Did you mean: resolution?", self.run_line("glossary resolutoin"))
        for line in listing.splitlines():
            self.assertLessEqual(len(line), 80, line)

    def test_status_and_start(self):
        text = self.run_line("status")
        self.assertIn("LOADED       1 structure", text)
        self.assertIn("you can also call it: mypep", text)
        self.assertIn("drawn as", text)
        self.assertIn("YOU COULD", text)
        cmd.delete("all")
        self.assertIn("NOTHING IS LOADED YET", self.run_line("status"))
        start = self.run_line("start")
        for phrase in ("grab haemoglobin", "ask ", "left-drag = rotate", "back", "info"):
            self.assertIn(phrase, start)
        for line in (text + start).splitlines():
            self.assertLessEqual(len(line), 80, line)

    def test_registered(self):
        for name in ("ask", "start", "status", "glossary"):
            self.assertIn(name, cmd.keyword)
            self.assertRegex(self.run_line("info"), r"(?m)^   %s\b" % name)


if __name__ == "__main__":
    unittest.main()
