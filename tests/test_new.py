"""Small molecules, plain-word selections, the fuller undo and short mode."""
import io
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from prospector.core import cache, chemcomp, http, report, search      # noqa: E402
from tests import fixture_transport                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


class Base(unittest.TestCase):
    def setUp(self):
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        os.environ.pop("ANTHROPIC_API_KEY", None)
        http.transport = fixture_transport.replaying

    def tearDown(self):
        http.transport = http._urllib_transport
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)


class MoleculeSearch(Base):
    def test_a_molecule_name_gives_the_molecule(self):
        result = search.find("cholesterol")
        self.assertEqual(result.kind, "molecule")
        molecule = result.molecule
        self.assertEqual(molecule.comp_id, "CLR")
        self.assertEqual(molecule.formula, "C27 H46 O")
        self.assertEqual(molecule.kind(), "steroid")
        self.assertTrue(molecule.exact)
        self.assertGreater(molecule.entries, 100)
        text = report.render(result)
        self.assertIn("BEST MATCH   Cholesterol  -  steroid, PDB code CLR", text)
        self.assertIn("grab protein ligand=CLR", text)
        self.assertIn("find cholesterol kind=protein", text)
        for line in text.split("\n"):
            self.assertLessEqual(len(line), 80, line)

    def test_everyday_names_and_drugs(self):
        glucose = search.find("glucose")
        self.assertEqual(glucose.molecule.comp_id, "BGC")
        self.assertEqual(glucose.molecule.kind(), "sugar")
        self.assertIn("official name: beta-D-glucopyranose", report.render(glucose))
        caffeine = search.find("caffeine")
        self.assertEqual(caffeine.molecule.comp_id, "CFF")
        self.assertEqual(search.find("ATP").molecule.comp_id, "ATP")
        # "Anti-Inflammatory Agents, Non-Steroidal" must not make ibuprofen a steroid
        self.assertNotEqual(chemcomp.lookup("ibuprofen").kind(), "steroid")

    def test_proteins_stay_proteins(self):
        self.assertEqual(search.find("p53").kind, "results")
        self.assertEqual(search.find("insulin").kind, "ambiguous")
        self.assertEqual(search.find("hemoglobin").kind, "results")
        self.assertEqual(search.find("4hhb").kind, "results")
        # a class word is not hijacked by a drug that merely mentions it
        self.assertNotEqual(search.find("kinase").kind, "molecule")

    def test_the_kind_filter_decides(self):
        self.assertEqual(search.find("cholesterol kind=protein").kind, "ambiguous")
        forced = search.find("HEM kind=molecule")
        self.assertEqual(forced.kind, "molecule")
        self.assertEqual(forced.molecule.comp_id, "HEM")
        missing = search.find("florbglorb kind=molecule")
        self.assertEqual(missing.kind, "none")
        self.assertIn("No small molecule", missing.message)
        self.assertEqual(search.find("p53 kind=banana").kind, "error")


class PubChemSource(Base):
    """PubChem is asked only when the PDB's dictionary does not know a name."""

    def test_a_molecule_that_is_only_in_pubchem(self):
        result = search.find("limonene kind=molecule")
        self.assertEqual(result.kind, "molecule")
        molecule = result.molecule
        self.assertEqual(molecule.source, "PubChem")
        self.assertTrue(molecule.has_shape())
        self.assertEqual(molecule.formula, "C10 H16")
        self.assertTrue(molecule.identifier.startswith("CID"))
        self.assertIn("record_type=3d", molecule.sdf_url())
        text = report.render(result)
        self.assertIn("from PubChem", text)
        self.assertIn("not in the PDB", text)
        self.assertNotIn("ligand=", text)            # it has no PDB code
        again = chemcomp.get(molecule.identifier)
        self.assertEqual(again.cid, molecule.cid)

    def test_the_pdb_dictionary_still_comes_first(self):
        self.assertEqual(search.find("cholesterol").molecule.source, "PDB")
        self.assertEqual(search.find("caffeine kind=molecule").molecule.source, "PDB")

    def test_proteins_are_not_hijacked_by_compounds_of_the_same_name(self):
        for name in ("ubiquitin", "calmodulin", "amylase", "ras", "hemoglobin"):
            self.assertNotEqual(search.find(name).kind, "molecule", name)
        # an enzyme is named after it: the protein keeps the name, with a hint
        limonene = search.find("limonene")
        self.assertEqual(limonene.kind, "results")
        self.assertIn("molecule limonene", report.render(limonene))

    def test_no_shape_and_not_found_are_said_plainly(self):
        missing = search.find("florbglorb kind=molecule")
        self.assertEqual(missing.kind, "none")
        self.assertIn("nor in PubChem", missing.message)
        shapeless = search.find("hemin kind=molecule")
        self.assertEqual(shapeless.kind, "none")
        self.assertIn("no 3D shape", shapeless.message)

    def test_only_an_organism_is_a_valid_search(self):
        result = search.find("a protein organism=ecoli")
        self.assertEqual(result.kind, "results")
        self.assertIn("Escherichia coli", result.top.entry.organism_text(None))


class PyMOLSide(Base):
    def setUp(self):
        Base.setUp(self)
        from pymol import cmd
        from prospector.core import store
        from prospector.pymol_layer import choose, commands
        self.cmd, self.commands, self.choose, self.store = cmd, commands, choose, store
        commands.register()
        cmd.delete("all")
        del commands._undo[:]
        commands._objects.clear()
        commands._aliases.clear()
        self.saved = commands._download_molecule
        commands._download_molecule = lambda molecule: os.path.join(
            HERE, "data", "CFF_ideal.sdf")

    def tearDown(self):
        self.commands._download_molecule = self.saved
        self.cmd.delete("all")
        Base.tearDown(self)

    def run_line(self, line):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            self.cmd.do(line, echo=0)
        return buffer.getvalue()

    # ---- small molecules ---------------------------------------------------
    def test_grab_loads_a_molecule_and_names_it(self):
        out = self.run_line("grab caffeine")
        self.assertIn("LOADED       Caffeine  -  small molecule, PDB code CFF", out)
        self.assertIn("ON SCREEN", out)
        self.assertIn("CFF", self.cmd.get_names("objects"))
        self.assertEqual(self.cmd.count_atoms("CFF and elem C"), 8)
        self.assertEqual(self.commands._known_object("caffeine"), "CFF")
        self.assertTrue(self.cmd.count_atoms("CFF and rep sticks") > 0)

    def test_molecule_command_and_card(self):
        out = self.run_line("molecule caffeine")
        self.assertIn("LOADED", out)
        card = self.run_line("about caffeine")
        self.assertIn("C8 H10 N4 O2", card)
        self.assertIn("MOLECULE", card)
        self.assertIn("idealised geometry", card)
        self.assertNotIn("KEY PLACES", card)

    def test_protein_only_commands_say_so(self):
        self.run_line("grab caffeine")
        for command in ("tour caffeine", "sequence caffeine", "scout caffeine histidines",
                        "whatif caffeine R248Q", "pocket caffeine", "features caffeine"):
            out = self.run_line(command)
            self.assertIn("is a small molecule, not a protein", out, command)
            self.assertIn("grab protein ligand=CFF", out, command)
        self.assertIn("spheres", self.run_line("viewas caffeine spheres"))

    def test_a_name_that_is_a_molecule_is_loaded_by_any_command(self):
        out = self.run_line("viewas caffeine sticks")
        self.assertIn("Loaded the small molecule 'caffeine'", out)
        self.assertIn("CFF", self.cmd.get_names("objects"))

    # ---- back removes what was added ----------------------------------------
    def test_back_takes_away_what_was_added(self):
        self.run_line("grab caffeine")
        self.cmd.fab("MKHACDEFGHIKLMNPQRSTVWYHCAAC", "pep", ss=1)
        self.commands._remember("mypep", "pep")
        self.run_line("measure mypep 3 10")
        self.assertIn("measure_1", self.cmd.get_names("all"))
        out = self.run_line("back")
        self.assertIn("Removed what was added since: measure_1", out)
        self.assertNotIn("measure_1", self.cmd.get_names("all"))
        self.assertIn("pep", self.cmd.get_names("objects"))

    def test_back_undoes_a_load(self):
        self.run_line("grab caffeine")
        self.assertIn("CFF", self.cmd.get_names("objects"))
        out = self.run_line("back")
        self.assertIn("loading CFF", out)
        self.assertEqual(self.cmd.get_names("objects"), [])
        self.assertIsNone(self.commands._known_object("caffeine"))

    # ---- choose --------------------------------------------------------------
    def test_translate(self):
        cases = {
            "histidines within 5 of the zinc":
                "byres ((X and (resn HIS)) within 5 of (X and (resn ZN)))",
            "charged residues in chain A": "X and (resn ASP+GLU+LYS+ARG and chain A)",
            "residues 100 to 120 except glycines":
                "X and (resi 100-120 and not (resn GLY))",
            "histidines and cysteines": "X and (resn HIS+CYS)",
            "the helices of chain B": "X and (ss h and chain B)",
            "chains A and B": "X and (chain A+B)",
            "the sulfur atoms of cysteines": "X and (resn CYS and elem S)",
            "water near the ligand":
                "byres ((X and (solvent)) within 5 of (X and (organic)))",
        }
        for words, expected in cases.items():
            self.assertEqual(self.choose.translate(words, {"ZN"}, "X")[0], expected, words)
        with self.assertRaises(self.choose.Unknown):
            self.choose.translate("florbs in chain A", {"ZN"}, "X")

    def test_choose_makes_a_selection_and_shows_the_pymol_line(self):
        self.cmd.fab("MKHACDEFGHIKLMNPQRSTVWYHCAAC", "pep", ss=1)
        self.commands._remember("mypep", "pep")
        out = self.run_line("choose histidines and cysteines in mypep as mine")
        self.assertIn("PYMOL        select mine, pep and (resn HIS+CYS)", out)
        self.assertIn("READ IT", out)
        self.assertIn("resn HIS = residues named HIS (histidine)", out)
        self.assertIn("mine", self.cmd.get_names("selections"))
        self.assertEqual(self.cmd.count_atoms("mine and name CA"), 6)
        out = self.run_line("choose florbs")
        self.assertIn("I don't know what 'florbs' means", out)
        self.run_line("back")
        self.run_line("back")
        self.assertNotIn("mine", self.cmd.get_names("selections"))
        self.assertIn("EXAMPLE", self.run_line("choose"))

    def test_choose_and_about_load_a_named_structure_by_themselves(self):
        out = self.run_line("choose the carbons in caffeine")
        self.assertIn("Loaded the small molecule 'caffeine'", out)
        self.assertIn("PYMOL        select chosen, CFF and (elem C)", out)
        self.assertEqual(self.cmd.count_atoms("chosen"), 8)
        self.cmd.delete("all")
        self.commands._objects.clear()
        self.commands._aliases.clear()
        card = self.run_line("about caffeine")
        self.assertIn("C8 H10 N4 O2", card)
        self.assertIn("CFF", self.cmd.get_names("objects"))      # loaded, not only described
        # "in chain A" is part of the description, never a protein to load
        self.cmd.delete("all")
        out = self.run_line("choose charged residues in chain A")
        self.assertIn("Nothing is loaded", out)

    # ---- clicking -------------------------------------------------------------
    def test_this_residue_means_the_one_clicked(self):
        from prospector.pymol_layer import ask
        self.cmd.fab("MKHACDEFGHIKLMNPQRSTVWYHCAAC", "pep", ss=1)
        self.commands._remember("mypep", "pep")
        command, reason = ask._interpret("why is this residue important")
        self.assertIsNone(command)
        self.assertIn("Click a residue", reason)
        self.cmd.select("sele", "pep and resi 5")           # what a click does
        self.assertEqual(self.commands._clicked_residue()[2:], ("5", "CYS"))
        self.assertEqual(ask._interpret("why is this residue important")[0], "explain 5")
        self.assertEqual(ask._interpret("zoom to the selected residue")[0], "focus 5")
        self.cmd.select("sele", "pep and resi 5-9")         # several: not "this residue"
        self.assertIsNone(self.commands._clicked_residue())

    # ---- short mode ----------------------------------------------------------
    def test_short_mode_prints_the_answer_and_more_prints_the_rest(self):
        self.cmd.fab("MKHACDEFGHIKLMNPQRSTVWYHCAAC", "pep", ss=1)
        self.commands._remember("mypep", "pep")
        full = self.run_line("scout mypep histidines")
        self.assertIn("WHAT IT IS", full)
        self.assertIn("short mode on", self.run_line("prospector short on"))
        short = self.run_line("scout mypep histidines")
        self.assertIn("SPOTTED", short)
        self.assertNotIn("WHAT IT IS", short)
        self.assertIn("more )", short)
        self.assertLess(len(short.strip().split("\n")), len(full.strip().split("\n")))
        again = self.run_line("more")
        self.assertIn("WHAT IT IS", again)
        self.assertIn("GO CLOSER", again)
        # cards and lists are never shortened
        self.assertIn("EXPERIMENT", self.run_line("info"))
        self.run_line("prospector short off")
        self.assertIn("WHAT IT IS", self.run_line("scout mypep histidines"))

    def test_shorten_keeps_notes(self):
        lines = ["Loaded the best match for 'x':  1ABC", "NOTE         important thing",
                 "NOW SHOWING  binding sites", "             A/1, A/2",
                 "ON SCREEN    4 residues", "             name tags", "             ribbon",
                 "WHAT IT IS   residues", "GO BACK      focus all", "MOUSE        drag"]
        kept, hidden = self.commands._shorten(lines)
        self.assertIn("NOTE         important thing", kept)
        self.assertIn("NOW SHOWING  binding sites", kept)
        self.assertNotIn("MOUSE        drag", kept)
        self.assertEqual(hidden, len(lines) - len(kept))


if __name__ == "__main__":
    unittest.main()
