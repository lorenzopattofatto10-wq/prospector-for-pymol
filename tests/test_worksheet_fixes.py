"""Fixes found while testing the haemoglobin worksheet."""
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


class WorksheetFixes(unittest.TestCase):
    def setUp(self):
        from prospector.pymol_layer import build, commands, typed
        self.build, self.commands, self.typed = build, commands, typed
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        os.environ.pop("ANTHROPIC_API_KEY", None)
        http.transport = fixture_transport.unplugged
        commands.register()
        cmd.delete("all")
        del commands._undo[:]
        commands._objects.clear()
        commands._aliases.clear()
        cmd.fab("MKHACDEFGHIKLMNPQRSTVWYHCAAC", "pep", ss=1, chain="A")
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

    def test_address_then_new_amino_acid(self):
        tidy = self.typed.tidy
        self.assertEqual(tidy("whatif", "A/3 W"), "A/3W")
        self.assertEqual(tidy("whatif", "A/3 to tryptophan"), "A/3W")
        self.assertEqual(tidy("whatif", "mypep A/his3 trp"), "mypep A/H3W")
        self.assertEqual(tidy("whatif", "A/H3W"), "A/H3W")          # already complete
        self.assertEqual(tidy("measure", "A/3 A/10"), "A/3 A/10")   # not a change

    def test_measure_to_a_metal_inside_a_molecule(self):
        x, y, z = cmd.get_coords("pep and chain A and resi 3 and name NE2")[0]
        cmd.pseudoatom("pep", name="FE", resn="HEM", resi="100", chain="A", elem="Fe",
                       pos=[x + 2.1, y, z], hetatm=1)
        out = self.run_line("measure pep A/3 to the iron")
        self.assertIn("DISTANCE     2.10 A", out)
        self.assertIn("HEM100 FE", out)

    def test_pocket_oxygen_means_the_heme(self):
        self.assertEqual(self.commands._ligand_word("oxygen"), "OXY")
        self.assertEqual(self.commands._ligand_word("the heme"), "HEM")
        self.assertEqual(self.commands._ligand_word("carbon monoxide"), "CMO")
        self.assertIsNone(self.commands._ligand_word("histidine"))

    def test_assembly_chains_get_plain_names_copy_by_copy(self):
        cmd.delete("all")
        for chain, seq in (("A", "MVLSPA"), ("A-2", "MVLSPA"), ("B", "MVHLTP"),
                           ("B-2", "MVHLTP")):
            cmd.fab(seq, "part", chain=chain)
            cmd.create("asm", "asm or part" if "asm" in cmd.get_names() else "part")
            cmd.delete("part")
        renamed = self.build._plain_chain_names("asm")
        chains = {}
        cmd.iterate("asm and name CA", "chains[chain] = resn",
                    space={"chains": chains})
        self.assertEqual(sorted(chains), ["A", "B", "C", "D"])
        # copy 1 = A, B; copy 2 = C, D - alpha, beta, alpha, beta like a normal entry
        self.assertEqual(chains["A"], chains["C"])
        self.assertEqual(chains["B"], chains["D"])
        self.assertIn(("A-2", "C"), renamed)

    def test_explain_suggests_the_files_numbering(self):
        import re
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "prospector", "pymol_layer", "tools.py")
        with open(path) as handle:
            source = handle.read()
        self.assertIn('replace(">", "%s" % resi)', source)
        self.assertTrue(re.search(r"model one with:  whatif", source))


if __name__ == "__main__":
    unittest.main()
