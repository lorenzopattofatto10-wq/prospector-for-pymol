"""Offline tests for scout, explain, measure, mutate, sequence, back, session,
report, tour jumping and pick=N."""
import io
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout

from prospector.core import cache, filters, http, report, search, sifts
from tests import fixture_transport, test_pymol

HAVE_PYMOL = test_pymol.HAVE_PYMOL
if HAVE_PYMOL:
    from pymol import cmd


class CoreAdditions(unittest.TestCase):
    def setUp(self):
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.replaying

    def tearDown(self):
        http.transport = http._urllib_transport
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def test_pick_chooses_another_row(self):
        self.assertEqual(filters.parse("p53 pick=2").pick, 2)
        self.assertEqual(filters.parse("p53 pick=12").limit, 12)     # list is long enough
        with self.assertRaises(filters.FilterError):
            filters.parse("p53 pick=0")
        first = search.find("p53")
        second = search.find("p53 pick=2")
        self.assertEqual(second.top.pdb_id, first.hits[1].pdb_id)
        self.assertEqual(second.command(), "fetch %s" % first.hits[1].pdb_id)
        self.assertIn("Just type its number, e.g.  2", report.render(first))
        self.assertEqual(report.choices(first)[1], "grab %s" % first.hits[1].pdb_id)

    def test_only_modelled_residues_count_as_coverage(self):
        result = search.find("p53")
        for hit in result.hits:
            self.assertLessEqual(hit.coverage, 1.0)
            self.assertGreater(hit.entry.modelled_fraction, 0)
        self.assertIn("is actually modelled", report.render(result))

    def test_textbook_structure_comes_first_in_a_text_search(self):
        result = search.find("kinase ligand=STI")
        self.assertTrue(result.top.textbook)
        self.assertIn("textbook structure of what you searched for", report.render(result))

    def test_reverse_residue_mapping(self):
        mapping = sifts.mapping("1PPB", prefer="rcsb")
        heavy = [c for c in mapping.chains if c.chain == "H"][0]
        self.assertEqual(heavy.uniprot_position("57"), 406)
        self.assertEqual(heavy.uniprot_position("60G"), 416)
        self.assertIsNone(heavy.uniprot_position("9999"))


@unittest.skipUnless(HAVE_PYMOL, "PyMOL is not installed in this Python")
class ToolTests(unittest.TestCase):
    SEQUENCE = "MKHACDEFGHIKLMNPQRSTVWYHCAAC"

    def setUp(self):
        from prospector.pymol_layer import commands, tools
        self.commands, self.tools = commands, tools
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.replaying
        commands.register()
        cmd.delete("all")
        del commands._undo[:]
        cmd.fab(self.SEQUENCE, "pep", ss=1)
        self.cwd = os.getcwd()

    def tearDown(self):
        os.chdir(self.cwd)
        http.transport = http._urllib_transport
        cmd.delete("all")
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def run_line(self, line):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            cmd.do(line, echo=0)
        return buffer.getvalue()

    def spotted(self):
        found = []
        cmd.iterate("pep_scout and name CA", "found.append(resn + resi)",
                    space={"found": found})
        return found

    def test_scout_amino_acids_groups_and_shapes(self):
        text = self.run_line("scout pep histidine")
        self.assertIn("SPOTTED      3 histidine residues in pep", text)
        self.assertEqual(self.spotted(), ["HIS3", "HIS10", "HIS24"])
        self.assertIn("WHAT IT IS   histidine (HIS", text)
        self.assertIn("ON SCREEN", text)
        self.assertGreater(cmd.count_atoms("pep and resn HIS and rep sticks"), 0)
        self.run_line("scout pep cysteines")                       # plural
        self.assertEqual(self.spotted(), ["CYS5", "CYS25", "CYS28"])
        self.assertEqual(cmd.count_atoms("pep and resn HIS and rep sticks"), 0)  # cleaned
        self.run_line("scout pep aromatic")
        self.assertEqual(sorted(set(r[:3] for r in self.spotted())),
                         ["HIS", "PHE", "TRP", "TYR"])
        self.run_line("scout pep negative")
        self.assertEqual(self.spotted(), ["ASP6", "GLU7"])
        self.assertIn("alpha-helices", self.run_line("scout pep helices"))
        self.assertGreater(len(self.spotted()), 20)                # fab built one helix
        self.run_line("scout pep termini")
        self.assertEqual(self.spotted(), ["MET1", "CYS28"])

    def test_scout_motifs_computed_things_and_help(self):
        self.run_line("scout pep motif CxxC")
        self.assertEqual(self.spotted(), ["CYS25", "ALA26", "ALA27", "CYS28"])
        self.assertIn("free cysteines", self.run_line("scout pep free cysteines"))
        self.assertEqual(self.spotted(), ["CYS5", "CYS25", "CYS28"])
        none = self.run_line("scout pep disulfides")
        self.assertIn("(0 bonds)", " ".join(none.split()))
        self.assertIn("none: pep has no metal ions", self.run_line("scout pep metals"))
        self.run_line("scout pep surface")
        self.assertGreater(len(self.spotted()), 5)
        menu = self.run_line("scout")
        for word in ("Amino acids", "By property", "Bonds", "Not protein", "Sequence"):
            self.assertIn(word, menu)
        typo = self.run_line("scout pep histidne")
        self.assertIn("Did you mean 'histidine'?", typo)
        for line in (menu + typo).splitlines():
            self.assertLessEqual(len(line), 80, line)

    def test_explain_measure_and_sequence(self):
        text = " ".join(self.run_line("explain pep 3").split())
        self.assertIn("RESIDUE HIS3", text)
        self.assertIn("AMINO ACID histidine - polar, sometimes charged", text)
        self.assertIn("WHERE in an alpha-helix", text)
        self.assertIn("TOUCHES", text)
        self.assertIn("no residue '99'", self.run_line("explain pep 99"))
        measured = self.run_line("measure pep 3 10")
        self.assertRegex(measured, r"DISTANCE     \d+\.\d\d A  between HIS3 and HIS10")
        self.assertIn("WHAT IT MEANS", measured)
        self.assertIn("measure_1", cmd.get_names())
        self.assertIn("needs two things", self.run_line("measure pep 3"))
        self.assertIn("can't find '99'", self.run_line("measure pep 3 99"))
        printed = self.run_line("sequence pep")
        self.assertIn(self.SEQUENCE, printed)
        self.assertIn("H" * 20, printed)
        self.assertIn("H = alpha-helix, E = beta-strand", printed)

    def test_mutate_leaves_the_original_untouched(self):
        atoms_before = cmd.count_atoms("pep")
        text = " ".join(self.run_line("mutate pep H3A").split())
        self.assertIn("MUTATION H3A in pep: histidine -> alanine", text)
        self.assertIn("CAUTION this is a quick model", text)
        self.assertEqual(cmd.count_atoms("pep"), atoms_before)
        names = []
        cmd.iterate("pep_H3A and resi 3 and name CA", "names.append(resn)",
                    space={"names": names})
        self.assertEqual(names, ["ALA"])
        wrong = self.run_line("mutate pep K3A")
        self.assertIn("Residue 3 is HIS in pep, not LYS (K)", wrong)
        self.assertIn("original letter + position + new letter", self.run_line("mutate pep"))

    def test_whatif_shows_two_panels_and_explains(self):
        text = self.run_line("whatif pep K2E")
        flat = " ".join(text.split())
        self.assertIn("WHAT IF K2E in pep", flat)
        self.assertIn("LEFT the original", flat)
        self.assertIn("RIGHT the changed version", flat)
        self.assertIn("Charge positive (+) -> negative (-) (the charge is reversed)", flat)
        self.assertIn("Size 5 -> 5 side-chain atoms", flat)
        self.assertIn("IN SHORT K2E reverses a charge", flat)
        self.assertIn("CAUTION a sketch, not a prediction", flat)
        self.assertEqual(str(cmd.get("grid_mode")), "1")
        for group in ("whatif_original", "whatif_changed"):
            self.assertEqual(cmd.get_type(group), "object:group")
        new, old = [], []
        cmd.iterate("pep_K2E and resi 2 and name CA", "new.append(resn)", space={"new": new})
        cmd.iterate("pep and resi 2 and name CA", "old.append(resn)", space={"old": old})
        self.assertEqual((old, new), (["LYS"], ["GLU"]))       # original untouched
        for line in text.splitlines():
            self.assertLessEqual(len(line), 80, line)
        self.assertIn("Back to one panel", self.run_line("whatif off"))
        self.assertEqual(str(cmd.get("grid_mode")), "0")
        self.assertNotIn("whatif_original", cmd.get_names("all"))
        self.assertIn("pep", cmd.get_names("objects"))
        # several changes at once, and polite refusals
        both = " ".join(self.run_line("whatif pep K2E H3A").split())
        self.assertIn("K2E, H3A in pep", both)
        self.assertIn("pep_K2E_H3A", cmd.get_names("objects"))
        self.run_line("back")
        self.assertEqual(str(cmd.get("grid_mode")), "0")       # back also leaves panels
        self.assertIn("Residue 2 is LYS in pep, not ARG (R)", self.run_line("whatif pep R2E"))
        self.assertIn("changes nothing", self.run_line("whatif pep K2K"))
        self.assertIn("original letter + position + new letter", self.run_line("whatif pep"))
        # a change typed the wrong way round is explained, never searched for
        confused = " ".join(self.run_line("whatif pep 2K 2E").split())
        self.assertIn("I can't read '2K 2E' as a change", confused)
        self.assertIn("Did you mean whatif pep K2E ?", confused)
        self.assertNotIn("Loaded the best match", confused)

    def test_changes_work_on_residues_with_two_modelled_shapes(self):
        # give residue 3 an alternate conformation, as many real entries have
        cmd.alter("pep and resi 3", "alt='A'")
        cmd.create("copyB", "pep and resi 3", zoom=0)
        cmd.alter("copyB", "alt='B'")
        cmd.create("pep", "pep or copyB", zoom=0)
        cmd.delete("copyB")
        self.run_line("mutate pep H3A")
        names = set()
        cmd.iterate("pep_H3A and resi 3", "names.add(resn)", space={"names": names})
        self.assertEqual(names, set(["ALA"]))

    def test_back_undoes_views_step_by_step(self):
        self.assertIn("Nothing to undo", self.run_line("back"))
        cmd.show_as("cartoon", "pep")
        self.run_line("viewas pep sticks")
        self.run_line("viewas pep spheres")
        self.assertTrue(cmd.count_atoms("pep and rep spheres"))
        text = self.run_line("back")
        self.assertIn("WENT BACK    to how it looked before:  viewas pep spheres", text)
        self.assertTrue(cmd.count_atoms("pep and rep sticks"))
        self.assertFalse(cmd.count_atoms("pep and rep spheres"))
        self.run_line("back")
        self.assertFalse(cmd.count_atoms("pep and rep sticks"))
        self.assertTrue(cmd.count_atoms("pep and rep cartoon"))
        self.assertIn("Nothing to undo", self.run_line("back"))

    def test_session_save_open_and_backup(self):
        self.assertIn("No saved sessions yet", self.run_line("session"))
        self.assertIn("SAVED        'homework'", self.run_line("session save homework"))
        self.assertIn("kept as a backup", self.run_line("session save homework"))
        self.assertIn("homework", self.run_line("session"))
        cmd.delete("all")
        cmd.fab("GGG", "other")
        opened = self.run_line("session open homework")
        self.assertIn("OPENED       'homework'", opened)
        self.assertIn("pep", cmd.get_names())
        self.assertIn("saved first as 'before_last_open'", opened)   # nothing was lost
        self.assertIn("Did you mean 'homework'?", self.run_line("session open homewrk"))
        cmd.delete("all")
        self.assertIn("nothing to save", self.run_line("session save empty"))

    def test_report_writes_one_html_file_and_restores_the_view(self):
        os.chdir(self.home)
        cmd.show_as("sticks", "pep")
        view = cmd.get_view()
        text = self.run_line("report pep")
        path = os.path.join(self.home, "pep_report.html")
        self.assertTrue(os.path.exists(path), text)
        page = open(path, encoding="utf-8").read()
        self.assertIn("<h2>Key places</h2>", page)
        self.assertGreaterEqual(page.count("data:image/png;base64,"), 2)
        self.assertTrue(cmd.count_atoms("pep and rep sticks"))       # view put back
        self.assertFalse(cmd.count_atoms("pep and rep surface"))
        self.assertEqual([round(v, 2) for v in cmd.get_view()], [round(v, 2) for v in view])
        self.run_line("report pep")
        self.assertTrue(os.path.exists(os.path.join(self.home, "pep_report_2.html")))

    def test_tour_jumps_to_a_stop(self):
        cmd.read_pdbstr(test_pymol.PyMOLTests.MINI, "1PPB")
        self.run_line("tour 1PPB")
        self.assertIn("[3/", self.run_line("tour 3"))
        self.assertIn("[2/", self.run_line("prev"))
        self.assertIn("This tour has stops 1 to", self.run_line("tour 9"))

    def test_the_last_structure_is_the_default(self):
        cmd.fab("ACDHHH", "second")
        self.run_line("scout pep histidine")
        text = self.run_line("measure 3 10")          # no protein named: uses pep
        self.assertIn("between HIS3 and HIS10", text)

    def test_new_commands_are_registered_and_listed(self):
        info = self.run_line("info")
        for name in ("scout", "explain", "measure", "mutate", "sequence", "back",
                     "session", "report", "whatif"):
            self.assertIn(name, cmd.keyword)
            self.assertIn("p_" + name, cmd.keyword)
            self.assertRegex(info, r"(?m)^   %s\b" % name)
        self.assertIs(cmd.keyword["undo"][0], cmd.undo)      # PyMOL's own is untouched


if __name__ == "__main__":
    unittest.main()
