"""`build`: DNA, RNA, peptides, sheets, membranes, micelles, assemblies.
All offline: the network is unplugged except for recorded PDB answers."""
import io
import math
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


class BuildTests(unittest.TestCase):
    def setUp(self):
        from prospector.pymol_layer import build, commands
        self.build, self.commands = build, commands
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        os.environ.pop("ANTHROPIC_API_KEY", None)
        http.transport = fixture_transport.replaying
        commands.register()
        cmd.delete("all")
        del commands._undo[:]
        commands._objects.clear()
        commands._aliases.clear()

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

    # ---- nucleic acids ---------------------------------------------------------
    def test_dna_double_helix(self):
        out = self.run_line("build dna ATGCGTACGTTA")
        self.assertIn("double-stranded DNA, 12 base pairs: dna_1", out)
        self.assertIn("IDEALISED MODEL", out)
        self.assertIn("3' TACGCATGCAAT 5'", out)            # the partner strand
        self.assertEqual(cmd.get_chains("dna_1"), ["A", "B"])
        self.assertEqual(cmd.count_atoms("dna_1 and name P and chain A"), 12)
        # rise of B-DNA: about 3.4 A per base pair
        first = cmd.get_coords("dna_1 and chain A and resi 2 and name C1'")[0]
        last = cmd.get_coords("dna_1 and chain A and resi 11 and name C1'")[0]
        self.assertTrue(25 < math.dist(first, last) < 37)
        self.assertIn("A-form", self.run_line("build dna ATGCGTACGTTA form=A"))
        self.assertIn("'U' is not one of them", self.run_line("build dna ATGCU"))
        self.assertIn("every U was read as T", self.run_line("build dna AUGCUU"))
        self.assertIn("Give the sequence", self.run_line("build dna"))

    def test_rna_is_a_single_strand(self):
        out = self.run_line("build rna GGCUAGCUAC")
        self.assertIn("a single RNA, 10 bases", out)
        self.assertIn("cannot build double-stranded RNA", out)
        self.assertEqual(cmd.get_chains("rna_1"), ["A"])
        self.assertIn("every T was read as U", self.run_line("build rna GGCTAGC"))

    # ---- peptides --------------------------------------------------------------
    def test_helix_strand_and_peptide(self):
        out = self.run_line("build helix AKLEELLKKLAELLKK")
        self.assertIn("16-residue peptide as an alpha-helix", out)
        cmd.dss("helix_1")
        helical = cmd.count_atoms("helix_1 and name CA and ss H")
        self.assertGreaterEqual(helical, 12)
        self.assertIn("helix_1_hbonds", cmd.get_names("all"))
        self.assertIn("extended beta-strand", self.run_line("build strand KVTVKVE"))
        self.assertIn("stretched-out chain", self.run_line("build peptide ACDEFGHIK"))
        self.assertIn("not one of the 20", self.run_line("build helix ABXZ"))

    def test_sheet_strands_lie_side_by_side_and_bond(self):
        out = self.run_line("build sheet KVTVKVE TVKVTVE KVEVKVT")
        self.assertIn("3-stranded beta-sheet", out)
        self.assertEqual(cmd.get_chains("sheet_1"), ["A", "B", "C"])
        for first, second in (("A", "B"), ("B", "C")):
            bonds = cmd.count_atoms("(sheet_1 and chain %s and name N) within 3.3 of "
                                    "(sheet_1 and chain %s and name O)" % (second, first))
            self.assertGreaterEqual(bonds, 2, (first, second))
        # flat: the third strand lies beyond the second, not on top of the first
        a = cmd.get_coords("sheet_1 and chain A and name CA").mean(axis=0)
        c = cmd.get_coords("sheet_1 and chain C and name CA").mean(axis=0)
        self.assertGreater(math.dist(a, c), 7.5)
        bumps = cmd.count_atoms("(sheet_1 and chain C and not hydro) within 2.4 of "
                                "(sheet_1 and chain A and not hydro)")
        self.assertEqual(bumps, 0)

    # ---- membranes and micelles ------------------------------------------------------
    def test_membrane_is_two_leaflets_tails_inward(self):
        out = self.run_line("build membrane popc 6")
        self.assertIn("POPC bilayer patch", out)
        self.assertIn("IDEALISED MODEL", out)
        self.assertIn("HOW REAL", out)
        upper = cmd.get_coords("membrane_1 and chain U and elem P")
        lower = cmd.get_coords("membrane_1 and chain L and elem P")
        self.assertEqual(len(upper), len(lower))
        self.assertTrue((upper[:, 2] > 15).all() and (lower[:, 2] < -15).all())
        tails = cmd.get_coords("membrane_1 and chain U and elem C")
        self.assertLess(tails[:, 2].min(), upper[:, 2].min())          # tails below heads
        # no two lipids overlap
        bumps = cmd.count_atoms("(membrane_1 and resi 1 and chain U) within 2.5 of "
                                "(membrane_1 and not (resi 1 and chain U))")
        self.assertEqual(bumps, 0)
        out = self.run_line("measure membrane_1 U/1/P L/1/P")
        self.assertIn("DISTANCE", out)

    def test_membrane_with_cholesterol_and_bad_words(self):
        out = self.run_line("build membrane dopc 6 cholesterol=30")
        self.assertIn("of them cholesterol", out)
        self.assertTrue(cmd.count_atoms("membrane_1 and resn CLR") > 0)
        self.assertIn("not one of them", self.run_line("build membrane banana"))
        self.assertIn("2 to 20 nm", self.run_line("build membrane popc 50"))

    def test_micelle_heads_out_tails_in(self):
        out = self.run_line("build micelle sds 40")
        self.assertIn("micelle of 40 SDS molecules", out)
        heads = cmd.get_coords("micelle_1 and elem S")
        tails = cmd.get_coords("micelle_1 and elem C")
        centre = cmd.get_coords("micelle_1").mean(axis=0)
        self.assertEqual(len(heads), 40)
        head_r = sum(math.dist(h, centre) for h in heads) / len(heads)
        tail_r = sum(math.dist(t, centre) for t in tails) / len(tails)
        self.assertGreater(head_r, tail_r + 5)
        self.assertIn("ldao", self.run_line("build micelle banana"))

    # ---- assemblies --------------------------------------------------------------------
    def test_operator_expressions(self):
        ops = self.build._operators
        self.assertEqual(ops("1"), [["1"]])
        self.assertEqual(ops("1,2,3"), [["1", "2", "3"]])
        self.assertEqual(ops("(1-3)"), [["1", "2", "3"]])
        self.assertEqual(ops("(1-2)(5,6)"), [["1", "2"], ["5", "6"]])
        self.assertEqual(ops("(X0)(1-2)"), [["X0"], ["1", "2"]])

    def test_coarse_view_applies_the_symmetry(self):
        cmd.fab("ACDEFGHIKLMNPQRSTVWY", "src", ss=1)
        path = os.path.join(self.home, "src.cif")
        cmd.alter("src", "segi='A'")
        cmd.save(path, "src")
        cmd.delete("src")
        turn = {"id": "2", "matrix_1_1": -1.0, "matrix_1_2": 0.0, "matrix_1_3": 0.0,
                "matrix_2_1": 0.0, "matrix_2_2": -1.0, "matrix_2_3": 0.0, "matrix_3_1": 0.0,
                "matrix_3_2": 0.0, "matrix_3_3": 1.0, "vector_1": 50.0, "vector_2": 0.0,
                "vector_3": 0.0}
        same = dict(turn, id="1", matrix_1_1=1.0, matrix_2_2=1.0, vector_1=0.0)
        plan = {"path": path, "operators": {
            "pdbx_struct_assembly_gen": [{"asym_id_list": ["A"], "oper_expression": "1,2"}],
            "pdbx_struct_oper_list": [same, turn]}}
        copies, step, balls = self.build._coarse(plan, "coarse")
        self.assertEqual((copies, step, balls), (2, 1, 40))
        one = cmd.get_coords("coarse and chain A1")
        two = cmd.get_coords("coarse and chain A2")
        # copy 2 = copy 1 turned 180 degrees about z and moved 50 A along x
        self.assertAlmostEqual(two[0][0], 50.0 - one[0][0], places=2)
        self.assertAlmostEqual(two[0][1], -one[0][1], places=2)

    def test_assembly_already_whole_and_ambiguous_names(self):
        mini = os.path.join(self.home, "4HHB.cif")
        cmd.fab("ACDEFGHIK", "tmp")
        cmd.save(mini, "tmp")
        cmd.delete("tmp")
        saved = self.commands._download_pdb
        self.commands._download_pdb = lambda pdb_id: mini
        try:
            out = self.run_line("build assembly 4HHB")
        finally:
            self.commands._download_pdb = saved
        self.assertIn("ALREADY WHOLE  4HHB", out)
        self.assertIn("tetrameric, 4 chains", out)
        out = self.run_line("build assembly insulin")
        self.assertIn("could mean several proteins, so I did not guess", out)
        self.assertIn("build assembly INS", out)

    # ---- the rest of Prospector ----------------------------------------------------------
    def test_built_objects_in_other_commands(self):
        self.run_line("build membrane popc 4")
        out = self.run_line("about membrane_1")
        self.assertIn("POPC bilayer patch", out)
        self.assertIn("you built, not a protein", self.run_line("tour membrane_1"))
        self.run_line("build helix AKLEELLKKLAE")
        self.assertIn("no database record", self.run_line("features helix_1"))
        self.assertIn("SPOTTED", self.run_line("scout helix_1 leucines"))
        self.assertIn("helix_1", cmd.get_names("objects"))
        self.run_line("back 2")         # undoes the scout, then the build itself
        self.assertNotIn("helix_1", cmd.get_names("objects"))
        self.assertIn("WHAT I CAN BUILD", self.run_line("build"))
        self.assertIn("build <what>", self.run_line("info"))

    def test_ask_understands_build_requests(self):
        from prospector.pymol_layer import ask
        cases = {
            "build a dna double helix of ATGCGTACGT": "build dna ATGCGTACGT",
            "make an A-form DNA helix GGCCAATTGG": "build dna GGCCAATTGG form=A",
            "create an alpha helix with the sequence AKLEELLKKLAE": "build helix AKLEELLKKLAE",
            "build a beta sheet from KVTVKVE and TVKVTVE": "build sheet KVTVKVE TVKVTVE",
            "build a lipid bilayer with 30% cholesterol": "build membrane cholesterol=30",
            "make a 10 nm popc membrane": "build membrane popc 10",
            "build a micelle of sds": "build micelle sds",
            "show the whole ferritin cage": "build assembly ferritin",
            "build the biological assembly of 1AON": "build assembly 1aon",
        }
        for sentence, command in cases.items():
            self.assertEqual(ask._interpret(sentence)[0], command, sentence)
        self.assertIsNone(ask._interpret("build me a sandwich")[0])


if __name__ == "__main__":
    unittest.main()
