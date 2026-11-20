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
        self.assertIn("BNG (beta-nonylglucoside): detergent", text)
        self.assertIn("1 water molecules", text)
        self.assertIn("1 of the 4 chains", text)
        self.assertIn("build assembly memb", text)

    def test_nothing_to_say_for_a_plain_complete_protein(self):
        self.assertEqual(self.commands._experiment_lines("memb", (1, 1, "monomeric")), [])


if __name__ == "__main__":
    unittest.main()
