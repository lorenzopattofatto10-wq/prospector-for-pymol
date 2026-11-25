"""Offline tests for  piece : part of a structure as a separate object."""
import io
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout

from prospector.core import cache, http, sifts
from tests import fixture_transport, test_pymol

HAVE_PYMOL = test_pymol.HAVE_PYMOL
if HAVE_PYMOL:
    from pymol import cmd


class _Protein(object):
    def __init__(self, chains):
        self.chains = chains


class _Annotation(object):
    def __init__(self, chains):
        self.proteins = [_Protein(chains)]


@unittest.skipUnless(HAVE_PYMOL, "PyMOL is not installed in this Python")
class PieceTests(unittest.TestCase):
    def setUp(self):
        from prospector.pymol_layer import commands, piece, places
        self.commands, self.piece, self.places = commands, piece, places
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        os.environ.pop("ANTHROPIC_API_KEY", None)
        http.transport = fixture_transport.replaying
        commands.register()
        cmd.delete("all")
        del commands._undo[:]
        # a two-chain test structure: chain A residues 1-30, chain B residues 101-120
        cmd.fab("ACDEFGHIKLMNPQRSTVWYACDEFGHIKL", "partA", ss=1, chain="A")
        cmd.fab("MKHACDEFGHIKLMNPQRST", "partB", ss=1, chain="B", resi=101)
        cmd.create("TEST", "partA or partB")
        cmd.delete("partA or partB")
        commands._remember("testprotein", "TEST")
        self.fetch = places.fetch_annotation

    def tearDown(self):
        self.places.fetch_annotation = self.fetch
        self.commands._objects.pop("TEST", None)
        http.transport = http._urllib_transport
        cmd.delete("all")
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def run_line(self, line):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            cmd.do(line, echo=0)
        return buffer.getvalue()

    def residues(self, obj):
        found = set()
        cmd.iterate("%s and name CA" % obj, "found.add((chain, int(resi)))",
                    space={"found": found})
        return found

    def test_the_usual_ways_of_writing_it(self):
        cases = {
            "piece A 19-20": ("TEST_A19-20", {("A", 19), ("A", 20)}),
            "piece testprotein A/5-8": ("TEST_A5-8", {("A", n) for n in range(5, 9)}),
            "piece residues 10 to 25 of chain A": ("TEST_A10-25",
                                                    {("A", n) for n in range(10, 26)}),
            "piece chain B": ("TEST_B", {("B", n) for n in range(101, 121)}),
            "piece /TEST//A/3": ("TEST_A3", {("A", 3)}),
            "piece A:2-4+6-7": ("TEST_A2-4_A6-7", {("A", n) for n in (2, 3, 4, 6, 7)}),
            "piece A 1-3 B 101-102": ("TEST_A1-3_B101-102",
                                      {("A", 1), ("A", 2), ("A", 3), ("B", 101),
                                       ("B", 102)}),
            "piece 110-112": ("TEST_B110-112", {("B", n) for n in (110, 111, 112)}),
            "piece chains A and B": ("TEST_A_B", None),
        }
        for line, (name, wanted) in cases.items():
            out = self.run_line(line)
            self.assertIn("MADE", out, line)
            self.assertIn(name, cmd.get_names("objects"), line)
            if wanted is not None:
                self.assertEqual(self.residues(name), wanted, line)
            self.assertIn("color red, %s" % name, out)
        # the original keeps every atom
        self.assertEqual(len(self.residues("TEST")), 50)

    def test_colour_and_short_pieces(self):
        out = self.run_line("piece A 19-40 color=blue")
        name = "TEST_A19-30"                        # chain A ends at 30
        self.assertIn(name, cmd.get_names("objects"))
        self.assertEqual(cmd.count_atoms("%s and name CA" % name),
                         cmd.count_atoms("%s and name CA and color blue" % name))
        self.assertIn("only residues 19-30 of 19-40", out)
        out = self.run_line("piece A 2-3 red")
        self.assertIn("sticks", out)
        self.assertGreater(cmd.count_atoms("TEST_A2-3 and rep sticks"), 0)
        self.assertEqual(cmd.count_atoms("TEST and chain A and resi 2-3 and rep cartoon"), 0)
        self.assertIn("COLOURED     red", out)

    def test_back_removes_the_piece(self):
        self.run_line("piece chain B")
        self.assertIn("TEST_B", cmd.get_names("objects"))
        self.run_line("back")
        self.assertNotIn("TEST_B", cmd.get_names("objects"))

    def test_helpful_errors(self):
        out = self.run_line("piece Z 1-5")
        self.assertIn("no chain 'Z'", out)
        self.assertIn("A (1-30), B (101-120)", out)
        out = self.run_line("piece A 200-210")
        self.assertIn("Chain A of TEST has no residues 200-210", out)
        self.assertIn("A (1-30)", out)
        out = self.run_line("piece")
        self.assertIn("Which part of TEST?", out)
        self.assertEqual(cmd.get_names("objects"), ["TEST"])
        cmd.delete("all")
        out = self.run_line("piece A 1-5")
        self.assertIn("Nothing is loaded", out)

    def test_uniprot_numbers(self):
        # UniProt numbers chain A from 51: file residue 1 = UniProt 51
        chain_map = sifts.ChainMap("P00000", "A", "1", [str(n) for n in range(1, 31)])
        chain_map.segments = [(51, 80, 1)]
        self.places.fetch_annotation = lambda identity: (_Annotation([chain_map]), {})
        self.commands._objects["TEST"] = ("pdb", "TEST")
        out = self.run_line("piece A 60-62")          # not in the file: UniProt numbers
        self.assertIn("UniProt numbering", out)
        self.assertEqual(self.residues("TEST_A10-12"), {("A", 10), ("A", 11), ("A", 12)})
        out = self.run_line("piece A 5-6")            # in the file: say what UniProt says
        self.assertIn("UniProt numbers the same residues 55-56", out)
        out = self.run_line("piece A 100-110")
        self.assertIn("UniProt residues 51-80 of P00000", out)

    def test_ask_makes_a_piece(self):
        out = self.run_line("ask make residues 19 to 22 of chain A a separate object")
        self.assertIn("piece A/19-22", out)
        self.assertIn("TEST_A19-22", cmd.get_names("objects"))
