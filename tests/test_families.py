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
from prospector.core import cache, http, rcsb, search          # noqa: E402
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

    def test_detergents_are_known_additives(self):
        for code in ("BNG", "BOG", "LDA", "LMT"):
            self.assertIn(code, rcsb.DETERGENTS)
            self.assertIn(code, rcsb.ADDITIVES)


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
        self.assertIn("fetch memb", rebuilt)

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
