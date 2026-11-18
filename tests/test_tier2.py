"""Offline tests for Tier 2: residue mapping (SIFTS), features, info cards."""
import os
import shutil
import tempfile
import unittest

from prospector.core import alphafold, cache, card, features, http, rcsb, sifts
from tests import fixture_transport


class Tier2Case(unittest.TestCase):
    def setUp(self):
        self.home = tempfile.mkdtemp(prefix="prospector-test-")
        cache.set_home(self.home)
        http.transport = fixture_transport.replaying

    def tearDown(self):
        http.transport = http._urllib_transport
        cache.set_home(None)
        shutil.rmtree(self.home, ignore_errors=True)


class MappingTests(Tier2Case):
    def chain(self, pdb_id, chain, prefer="rcsb"):
        mapping = sifts.mapping(pdb_id, prefer=prefer)
        return [c for c in mapping.chains if c.chain == chain][0]

    def test_thrombin_offset_and_insertion_codes(self):
        """1PPB uses chymotrypsin numbering: the catalytic triad His57/Asp102/
        Ser195 is His406/Asp462/Ser568 in UniProt P00734."""
        heavy = self.chain("1PPB", "H")
        self.assertEqual([heavy.resi(p) for p in (406, 462, 568)], ["57", "102", "195"])
        self.assertEqual(heavy.resi(364), "16")       # first residue of the heavy chain
        self.assertEqual(heavy.resi(416), "60G")      # inside the 60A-60I insertion loop
        self.assertIsNone(heavy.resi(100))            # Gla domain: not in this crystal
        self.assertIn("irregular numbering", heavy.describe_numbering())
        light = self.chain("1PPB", "L")
        self.assertEqual(light.resi(328), "1H")
        self.assertEqual(light.resi(336), "1")

    def test_trypsin(self):
        chain = self.chain("2PTC", "E")
        self.assertEqual([chain.resi(p) for p in (63, 107, 200)], ["57", "102", "195"])

    def test_hemoglobin_constant_offset(self):
        alpha = self.chain("4HHB", "A")
        self.assertEqual(alpha.accession, "P69905")
        self.assertEqual(alpha.resi(88), "87")        # proximal histidine
        self.assertEqual(alpha.describe_numbering(), "PDB number = UniProt number - 1")
        self.assertIsNone(alpha.resi(1))              # initiator Met is removed

    def test_pdbe_and_rcsb_agree(self):
        """Both sources derive from SIFTS; they must give identical numbers.
        (If PDBe was unreachable when the fixtures were recorded, the PDBe
        request falls back to RCSB and the test still holds.)"""
        for pdb_id in ("1PPB", "4HHB"):
            first = sifts.mapping(pdb_id, prefer="pdbe")
            second = sifts.mapping(pdb_id, prefer="rcsb")
            self.assertEqual(sorted((c.accession, c.chain) for c in first.chains),
                             sorted((c.accession, c.chain) for c in second.chains))
            for a in first.chains:
                b = [c for c in second.chains
                     if (c.accession, c.chain) == (a.accession, a.chain)][0]
                for start, end, _ in a.segments:
                    for position in (start, (start + end) // 2, end):
                        self.assertEqual(a.resi(position), b.resi(position),
                                         (pdb_id, a.chain, position))

    def test_alphafold_identity(self):
        chain = sifts.identity("P04637", 1, 393).chains[0]
        self.assertEqual((chain.resi(1), chain.resi(248), chain.resi(394)),
                         ("1", "248", None))


class FeatureTests(Tier2Case):
    def test_thrombin_features(self):
        annotation = features.annotate("1PPB")
        protein = annotation.proteins[0]
        self.assertEqual(protein.protein.accession, "P00734")
        active = [f for f in protein.features if f.group == "active"]
        self.assertEqual([f.placed for f in active],
                         [{"H": ["57"]}, {"H": ["102"]}, {"H": ["195"]}])
        # interchain disulfide: one cysteine in each chain
        bridge = [f for f in protein.features if f.group == "disulfide"
                  and (f.start, f.end) == (336, 482)][0]
        self.assertEqual(bridge.placed, {"L": ["1"], "H": ["122"]})
        self.assertEqual(bridge.pdb_text(), "H/122 + L/1")
        # features of the Gla / kringle domains are honestly reported as absent
        self.assertTrue(any(not f.placed for f in protein.features if f.group == "domain"))
        text = features.render(annotation)
        self.assertIn("H/57", text)
        self.assertIn("Not in this structure", text)

    def test_selection_plan_and_checks(self):
        annotation = features.annotate("1PPB")
        plan = dict((item[0], item) for item in annotation.selections())
        self.assertEqual(plan["active"][3], {"H": ["57", "102", "195"]})
        self.assertEqual(features.selection_text("1PPB", plan["active"][3]),
                         "1PPB and ((chain H and resi 57+102+195))")
        self.assertIn("dom_Peptidase_S1", plan)
        self.assertIn(("H", "57", "H", 406), annotation.checks())
        self.assertIn(("H", "195", "S", 568), annotation.checks())

    def test_two_proteins_in_one_entry(self):
        annotation = features.annotate("4HHB", features.parse_types("binding"))
        accessions = sorted(p.protein.accession for p in annotation.proteins)
        self.assertEqual(accessions, ["P68871", "P69905"])
        merged = annotation.selections()[0][3]
        self.assertEqual(sorted(merged), ["A", "B", "C", "D"])
        self.assertIn("87", merged["A"])      # proximal His of the alpha chain
        self.assertIn("92", merged["B"])      # proximal His of the beta chain

    def test_model_features_use_uniprot_numbers(self):
        annotation = features.annotate_model("P04637", ["binding"])
        zinc = [f.placed for f in annotation.proteins[0].features]
        self.assertEqual(zinc, [{"A": ["176"]}, {"A": ["179"]}, {"A": ["238"]},
                                {"A": ["242"]}])

    def test_types_and_errors(self):
        self.assertEqual(features.parse_types("active, Binding+ss"),
                         ["active", "binding", "disulfide"])
        self.assertNotIn("site", features.parse_types(""))
        self.assertIn("site", features.parse_types("all"))
        with self.assertRaises(features.FeatureError):
            features.parse_types("colour")
        with self.assertRaises(features.FeatureError) as caught:
            features.annotate("1BNA")              # DNA only
        self.assertIn("no chain linked to UniProt", str(caught.exception))

    def test_disease_variant_filter(self):
        self.assertTrue(features._is_disease("in LFS; germline mutation; dbSNP:rs28934578"))
        self.assertFalse(features._is_disease("in dbSNP:rs1800371"))
        self.assertFalse(features._is_disease(""))

    def test_negative_residue_numbers_are_escaped(self):
        self.assertEqual(features.resi_text(["-2", "-1", "1", "60A"]), "\\-2+\\-1+1+60A")


class CardTests(Tier2Case):
    def test_entry_card(self):
        text = card.entry_card("1PPB")
        lines = text.splitlines()
        for label, value in (("Method", "X-ray diffraction"), ("Resolution", "1.92 A"),
                             ("R-work", "0.156"), ("R-free", "not reported"),
                             ("Name", "Prothrombin"), ("Gene", "F2"),
                             ("Organism", "Homo sapiens"), ("PubMed", "2583108"),
                             ("Year", "1989")):
            self.assertIn("   %-13s%s" % (label, value), lines)
        for heading in (" EXPERIMENT", " PROTEIN", " CHAINS", " LIGANDS", " CITATION"):
            self.assertIn(heading, lines)
        self.assertIn("P00734  (622 amino acids)", text)
        chain_h = [l for l in lines if l.startswith("   H  ")][0]
        self.assertIn("259 aa", chain_h)
        self.assertIn("364-622", chain_h)
        self.assertTrue(any(l.startswith("   0G6 ") for l in lines))
        self.assertTrue(all(len(l) <= 80 for l in lines), max(lines, key=len))
        self.assertIn("https://doi.org/10.1016/0022-2836(87)90636-X", card.entry_card("2LZM"))

    def test_model_card(self):
        text = card.model_card(alphafold.prediction("P04637"))
        self.assertIn("a computed prediction, not an experiment", text)
        self.assertIn("Cellular tumor antigen p53", text)
        self.assertNotIn("no experimental structure exists", text)
        self.assertIn(" MODEL", text.splitlines())

    def test_ligand_names(self):
        names = rcsb.ligand_names(rcsb.entry_detail("1IEP"))
        self.assertIn("STI", names)


if __name__ == "__main__":
    unittest.main()
