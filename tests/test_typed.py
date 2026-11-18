"""What is typed after a command: every usual spelling and word order."""
import io
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pymol import cmd                                          # noqa: E402

from prospector.core import cache, http                        # noqa: E402
from tests import fixture_transport                            # noqa: E402


class TypedTests(unittest.TestCase):
    def setUp(self):
        from prospector.pymol_layer import commands, typed
        self.commands, self.typed = commands, typed
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        os.environ.pop("ANTHROPIC_API_KEY", None)
        http.transport = fixture_transport.unplugged      # any lookup would fail loudly
        commands.register()
        cmd.delete("all")
        del commands._undo[:]
        commands._objects.clear()
        commands._aliases.clear()
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

    def test_residues_in_every_spelling(self):
        tidy = self.typed.tidy
        for text in ("3", "H3", "his3", "histidine 3", "residue 3", "the histidine 3"):
            self.assertEqual(tidy("explain", text), "3", text)
        self.assertEqual(tidy("explain", "mypep H3"), "mypep 3")
        self.assertEqual(tidy("measure", "3 and 10"), "3 10")
        self.assertEqual(tidy("measure", "his 3 to his 10"), "3 10")
        self.assertEqual(tidy("measure", "residue 3 and residue 10"), "3 10")
        self.assertEqual(tidy("measure", "A/3 to A/10"), "A/3 A/10")

    def test_changes_in_every_spelling(self):
        tidy = self.typed.tidy
        for text in ("H3W", "his3trp", "H3 to W", "3 W", "histidine 3 to tryptophan",
                     "change 3 to tryptophan", "his 3 trp"):
            self.assertIn(tidy("whatif", text), ("H3W", "3W"), text)
        self.assertEqual(tidy("whatif", "mypep H3W"), "mypep H3W")
        self.assertEqual(tidy("whatif", "H3W K2A"), "H3W K2A")

    def test_protein_names_that_look_like_residues_are_kept(self):
        tidy = self.typed.tidy
        # the AI's answer to "p53 normal vs breast cancer side by side"
        self.assertEqual(tidy("whatif", "p53 R175H"), "p53 R175H")
        self.assertEqual(tidy("whatif", "R175H p53"), "R175H p53")
        self.assertEqual(tidy("mutate", "p53 R248Q"), "p53 R248Q")
        self.assertEqual(tidy("explain", "p53 248"), "p53 248")
        self.assertEqual(tidy("measure", "p53 176 238"), "p53 176 238")
        self.assertEqual(tidy("whatif", "p21 K16A"), "p21 K16A")
        self.assertEqual(tidy("whatif", "kras G12D"), "kras G12D")
        # a residue of the structure on screen is still read as one
        self.assertEqual(tidy("whatif", "H3 to W"), "H3W")
        self.assertEqual(tidy("whatif", "H3"), "3")

    def test_filler_words_and_word_order(self):
        tidy = self.typed.tidy
        self.assertEqual(tidy("scout", "cysteines in mypep"), "cysteines mypep")
        self.assertEqual(tidy("scout", "show the histidines"), "histidines")
        self.assertEqual(tidy("scout", "salt bridges"), "salt bridges")
        self.assertEqual(tidy("look", "colour by chain"), "chain")
        self.assertEqual(tidy("compare", "1HHO with 4HHB"), "1HHO , 4HHB")
        self.assertEqual(tidy("compare", "1HHO, 4HHB"), "1HHO, 4HHB")
        # commands with free text are left alone
        self.assertEqual(tidy("grab", "protein of the week"), "protein of the week")
        self.assertEqual(tidy("choose", "histidines in chain A"), "histidines in chain A")

    def test_keywords_and_numbers_are_never_looked_up_as_proteins(self):
        not_a_name = self.typed.not_a_name
        for text in ("rainbow", "and", "to", "Q", "by chain", "248", "sticks", "A/248"):
            self.assertIsNotNone(not_a_name(text), text)
        for text in ("histidnes", "zinc", "arginine"):       # typo, bound molecule, amino acid
            self.assertIsNotNone(not_a_name(text, strict=True), text)
        for text in ("haemoglobin", "p53", "green fluorescent protein", "4HHB", "thrombin"):
            self.assertIsNone(not_a_name(text, strict=True), text)
        fetched = self.commands._fetch_query("rainbow")
        self.assertIn("did not search", fetched["error"])

    def test_commands_work_however_they_are_typed(self):
        # the transport is unplugged: any attempt to look something up would show
        for line, expected in (
                ("look mypep rainbow", "Applied 'rainbow' to pep"),
                ("look rainbow mypep", "Applied 'rainbow' to pep"),
                ("look ranbow", "Did you mean 'rainbow'"),
                ("scout cysteines in mypep", "cysteine residues in pep"),
                ("scout histidnes", "Did you mean 'histidines'"),
                ("scout lysines and histidines", "one kind of thing at a time"),
                ("explain H3", "HIS3"),
                ("explain histidine 3", "HIS3"),
                ("measure 3 and 10", "DISTANCE"),
                ("measure his 3 to his 10", "DISTANCE"),
                ("whatif H3 to W", "H3W"),
                ("whatif his3trp", "H3W"),
                ("viewas stiks", "Did you mean 'sticks'")):
            out = self.run_line(line)
            self.assertIn(expected, out, line)
            self.assertNotIn("Could not reach the network", out, line)
            self.assertNotIn("Loaded the best match", out, line)
        self.assertEqual([n for n in cmd.get_names("objects") if n not in
                          ("pep", "pep_H3W") and not n.startswith(("measure", "whatif"))],
                         [])

    def test_a_number_means_the_protein_residue_not_a_water_with_that_number(self):
        from prospector.pymol_layer import tools
        cmd.pseudoatom("pep", name="O", resn="HOH", resi="3", chain="W", hetatm=1,
                       pos=[30, 30, 30])
        selection = tools._spec_selection("pep", "3")
        names = set()
        cmd.iterate(selection, "names.add(resn)", space={"names": names})
        self.assertEqual(names, {"HIS"})

    def test_session_save_has_a_default_name(self):
        out = self.run_line("session save")
        self.assertIn("SAVED        'pep'", out)


if __name__ == "__main__":
    unittest.main()
