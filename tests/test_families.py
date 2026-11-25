"""Family names give a list; experiment molecules and partial files are explained."""
import io
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pymol import cmd                                          # noqa: E402
from prospector.core import cache, http, rcsb, report, resolve, search  # noqa: E402
from tests import fixture_transport                            # noqa: E402


class _Protein(object):
    def __init__(self, name):
        self.name = name


class _Group(object):
    def __init__(self, *names):
        self.members = [_Protein(n) for n in names]


class Families(unittest.TestCase):
    def test_family_members_are_named_after_the_word(self):
        aquaporins = [_Group("Aquaporin-1"), _Group("Aquaporin Z"),
                      _Group("Lens fiber major intrinsic protein")]
        self.assertEqual(search._family_members("aquaporin", aquaporins), 2)
        ribosome = [_Group("Ribosome biogenesis protein ERB1"),
                    _Group("mRNA turnover protein 4 homolog")]
        self.assertEqual(search._family_members("ribosome", ribosome), 0)
        self.assertEqual(search._family_members("myosin", [_Group("Myosin-7"),
                                                           _Group("Myosin-9")]), 2)

    def test_a_family_word_must_head_the_name(self):
        heads = resolve._heads
        self.assertTrue(heads("Caspase-3", "caspase"))
        self.assertTrue(heads("Tubulin alpha-1B chain", "tubulin"))
        self.assertTrue(heads("Actin, cytoplasmic 1", "actin"))
        self.assertTrue(heads("Histone H3.1", "histone"))
        self.assertTrue(heads("Cyclin-dependent kinase 2", "kinase"))
        self.assertTrue(heads("G protein-coupled receptor 151", "G protein-coupled receptor"))
        self.assertTrue(heads("Dihydroorotate dehydrogenase (quinone), mitochondrial",
                              "dehydrogenase"))
        self.assertFalse(heads("Cyclin-dependent kinase 2", "cyclin"))
        self.assertFalse(heads("Histone deacetylase 1", "histone"))
        self.assertFalse(heads("Nuclear receptor coactivator 2", "nuclear receptor"))
        self.assertFalse(heads("Glycogen synthase kinase-3 beta", "synthase"))
        self.assertFalse(heads("dTDP-glucose 4,6-dehydratase", "glucose"))

    def test_how_a_protein_belongs_to_a_family(self):
        class Entry(object):
            name, families, keywords = "Adenosine receptor A2a", [], []
        receptor = Entry()
        phrases = resolve.family_phrases("GPCR")
        self.assertEqual(phrases, ["G protein-coupled receptor"])
        self.assertIsNone(resolve.family_member(receptor, phrases))
        receptor.families = ["Belongs to the G-protein coupled receptor 1 family"]
        self.assertEqual(resolve.family_member(receptor, phrases), "family")
        receptor.families, receptor.keywords = [], ["G-protein coupled receptor"]
        self.assertEqual(resolve.family_member(receptor, phrases), "keyword")
        self.assertEqual(resolve.family("ribosome"), [])        # a complex, not a family
        self.assertEqual(resolve.family("zinc finger"), [])     # a part of many proteins
        self.assertEqual(resolve.family("ubiquitin"), [])       # one protein, four genes

    def test_detergents_are_known_additives(self):
        for code in ("BNG", "BOG", "LDA", "LMT"):
            self.assertIn(code, rcsb.DETERGENTS)
            self.assertIn(code, rcsb.ADDITIVES)


class FamilyLists(unittest.TestCase):
    """Recorded answers: a family name gives a list, a protein loads, a filter loads."""

    def setUp(self):
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.replaying

    def tearDown(self):
        http.transport = http._urllib_transport
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def members(self, result):
        return [g.lead.gene for g in result.resolution.choices]

    def test_family_names_give_their_best_studied_members(self):
        kinase = search.find("kinase")
        self.assertEqual(kinase.kind, "ambiguous")
        self.assertTrue(kinase.resolution.family)
        self.assertEqual(self.members(kinase)[0], "CDK2")
        text = report.render(kinase, "grab")
        self.assertIn("'kinase' is a family of many proteins", text)
        self.assertIn("grab CDK2", text)
        self.assertIn("find kinase mode=text", text)
        self.assertNotIn("molecule kinase", text)
        self.assertEqual(report.choices(kinase)[0], "grab P24941")
        gpcr = search.find("GPCR")                     # not the plant protein GCR2
        self.assertEqual(gpcr.kind, "ambiguous")
        self.assertIn("ADORA2A", self.members(gpcr))
        self.assertIn("ADRB2", self.members(gpcr))
        globin = search.find("globin")                 # few human globins: all organisms
        self.assertEqual(globin.kind, "ambiguous")
        self.assertIn("MB", self.members(globin))
        self.assertIn("all organisms", " ".join(globin.notes))

    def test_proteins_complexes_and_filters_still_load(self):
        self.assertEqual(search.find("hemoglobin").kind, "results")
        self.assertEqual(search.find("p53").kind, "results")
        filtered = search.find("kinase ligand=STI")
        self.assertEqual(filtered.kind, "results")
        self.assertIn("STI", filtered.top.entry.ligands)
        nucleosome = search.find("nucleosome")
        self.assertEqual(nucleosome.kind, "results")
        self.assertEqual(nucleosome.top.pdb_id, "1KX5")


class ExperimentLines(unittest.TestCase):
    def setUp(self):
        from prospector.pymol_layer import commands
        self.commands = commands
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.unplugged
        commands.register()
        cmd.delete("all")
        cmd.fab("MKHACDEFG", "memb", chain="A")

    def tearDown(self):
        http.transport = http._urllib_transport
        cmd.delete("all")
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def test_detergent_water_and_partial_file_are_explained(self):
        cmd.pseudoatom("memb", name="C1", resn="BNG", resi="801", chain="A", hetatm=1)
        cmd.pseudoatom("memb", name="O", resn="HOH", resi="900", chain="A", hetatm=1)
        text = "\n".join(self.commands._experiment_lines("memb", (4, 1, "tetrameric")))
        self.assertIn("DETERGENT    BNG (beta-nonylglucoside)", text)
        self.assertIn("1 water molecule ", text)
        self.assertIn("1 of the 4 chains", text)
        self.assertIn("build assembly memb", text)
        rebuilt = "\n".join(self.commands._experiment_lines("memb", (4, 1, "tetrameric"),
                                                            assembled=True))
        self.assertIn("WHOLE MOLECULE", rebuilt)
        self.assertIn("this is the complete molecule", rebuilt)
        self.assertIn("unit memb file", rebuilt)

    def test_nothing_to_say_for_a_plain_complete_protein(self):
        self.assertEqual(self.commands._experiment_lines("memb", (1, 1, "monomeric")), [])


if __name__ == "__main__":
    unittest.main()


class WholeMolecule(unittest.TestCase):
    def setUp(self):
        from prospector.pymol_layer import build, commands, places
        self.build, self.commands, self.places = build, commands, places
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.unplugged
        commands.register()
        cmd.delete("all")

    def tearDown(self):
        http.transport = http._urllib_transport
        cmd.delete("all")
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def test_first_copy_keeps_its_chain_names(self):
        for chain in ("H", "L", "H-2", "L-2"):
            cmd.fab("MVLSPA", "part", chain=chain)
            cmd.create("asm", "asm or part" if "asm" in cmd.get_names() else "part")
            cmd.delete("part")
        renamed = dict(self.build._plain_chain_names("asm"))
        self.assertEqual(sorted(cmd.get_chains("asm")), ["A", "B", "H", "L"])
        self.assertEqual(renamed, {"H-2": "A", "L-2": "B"})

    def test_copies_get_the_annotation_of_their_original(self):
        from prospector.core import features as feat, sifts
        annotation = feat.Annotation("1ABC", "test")
        mapping = sifts.ChainMap("P00001", "A", 1, ["1", "2", "3"])
        feature = feat.Feature("site", "Active site", 2, 2, "", "P00001")
        feature.placed = {"A": ["2"]}
        annotation.proteins.append(feat.ProteinFeatures(None, {}, [mapping], [feature]))
        copied = self.places._onto_copies(annotation, {"C": "A", "D": "A"})
        self.assertEqual(sorted(m.chain for m in copied.proteins[0].chains), ["A", "C", "D"])
        self.assertEqual(copied.proteins[0].features[0].placed,
                         {"A": ["2"], "C": ["2"], "D": ["2"]})
        self.assertEqual(feature.placed, {"A": ["2"]})       # the cached one is untouched

    def test_ions_are_not_called_fake(self):
        cmd.fab("MKHACDEFG", "chan", chain="A")
        cmd.pseudoatom("chan", name="K", resn="K", resi="301", chain="A", elem="K", hetatm=1)
        text = "\n".join(self.commands._experiment_lines("chan", None))
        self.assertIn("IONS", text)
        self.assertIn("part of how the protein works", text)
        self.assertNotIn("not real ligands", text)

    def test_letter_o_and_digit_zero(self):
        self.assertEqual(self.build._id_variants("8Y80"), ["8Y8O"])
        self.assertIn("1HH0", self.build._id_variants("1HHO"))


class UnitViews(unittest.TestCase):
    def setUp(self):
        from prospector.pymol_layer import commands
        self.commands = commands
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.unplugged
        commands.register()
        cmd.delete("all")
        commands._objects.clear()
        for chain in ("A", "B", "C", "D"):
            cmd.fab("MKHACDEFG", "part", chain=chain)
            cmd.create("TET", "TET or part" if "TET" in cmd.get_names() else "part")
            cmd.delete("part")
        commands._objects["TET"] = ("pdb", "TET")
        commands._last[0] = "TET"

    def tearDown(self):
        http.transport = http._urllib_transport
        cmd.delete("all")
        self.commands._assembled.pop("TET", None)
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)

    def run_line(self, line):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            cmd.do(line, echo=0)
        return buffer.getvalue()

    def test_one_chain_and_back(self):
        out = self.run_line("unit one")
        self.assertIn("one chain: chain A of TET", out)
        self.assertEqual(cmd.get_names("objects", enabled_only=1), ["TET_one"])
        self.assertEqual(cmd.get_chains("TET_one"), ["A"])
        out = self.run_line("unit whole")
        self.assertIn("4 chains", out)
        self.assertEqual(cmd.get_names("objects", enabled_only=1), ["TET"])

    def test_file_view_when_the_file_is_already_whole(self):
        out = " ".join(self.run_line("unit TET file").split())
        self.assertIn("already is exactly the PDB file", out)

    def test_sentences(self):
        from prospector.pymol_layer import ask
        for text, want in (("show the monomer", "unit one"), ("just one chain", "unit one"),
                           ("back to the tetramer", "unit whole"),
                           ("show the pdb file as it is", "unit file")):
            self.assertEqual(ask._interpret(text)[0], want, text)
        self.assertEqual(ask._interpret("show the whole ferritin cage")[0],
                         "build assembly ferritin")
