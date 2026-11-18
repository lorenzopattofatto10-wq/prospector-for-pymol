"""Small molecules: steroids, lipids, sugars, nucleotides, drugs, cofactors.

Source: the PDB's Chemical Component Dictionary (RCSB).  It holds every small
molecule that has ever been seen in a PDB structure - about 45,000 of them -
each with a name, synonyms, formula and a set of ideal 3D coordinates.

    Search API v2   return_type "mol_definition", service "text_chem"
    Data API        GraphQL  chem_comps(comp_ids: [...])
    Files           https://files.rcsb.org/ligands/download/<ID>_ideal.sdf

A molecule that was never crystallised with a protein or nucleic acid is not
in this dictionary.  Membranes and bilayers are not single molecules and are
not covered either.
"""
from __future__ import annotations

import re

from . import cache, http

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
SDF_URL = "https://files.rcsb.org/ligands/download/%s_ideal.sdf"

# Everyday names -> the dictionary's ID, for molecules whose official name is
# something a student would not type ("glucose" is filed as
# "beta-D-glucopyranose").  Each ID is checked against the live dictionary in
# tests/test_live.py.
COMMON = {
    "glucose": "BGC", "d-glucose": "BGC", "beta-glucose": "BGC", "alpha-glucose": "GLC",
    "fructose": "FRU", "galactose": "GAL", "mannose": "MAN",
    "ribose": "RIB", "glycerol": "GOL", "urea": "URE",
    "atp": "ATP", "adp": "ADP", "amp": "AMP", "gtp": "GTP", "gdp": "GDP", "camp": "CMP",
    "nad": "NAD", "nad+": "NAD", "nadh": "NAI", "nadp": "NAP", "nadph": "NDP",
    "fad": "FAD", "fmn": "FMN", "coa": "COA", "coenzyme a": "COA", "acetyl-coa": "ACO",
    "heme": "HEM", "haem": "HEM", "heme b": "HEM", "chlorophyll": "CLA",
    "chlorophyll a": "CLA", "retinal": "RET", "biotin": "BTN", "water": "HOH",
    "cholesterol": "CLR", "testosterone": "TES", "estradiol": "EST", "oestradiol": "EST",
    "cortisol": "HCY", "progesterone": "STR", "popc": "POV", "palmitate": "PLM",
    "palmitic acid": "PLM", "oleic acid": "OLA", "stearic acid": "STE",
    "aspirin": "AIN", "caffeine": "CFF", "ethanol": "EOH", "citrate": "CIT",
    "pyruvate": "PYR", "dopamine": "LDP", "serotonin": "SRO", "adrenaline": "ALE",
    "glutathione": "GSH", "vitamin c": "ASC", "ascorbic acid": "ASC",
}

_FIELDS = """{ chem_comps(comp_ids: %s) {
  chem_comp { id name formula formula_weight type pdbx_formal_charge }
  rcsb_chem_comp_info { atom_count atom_count_heavy atom_count_chiral bond_count_aromatic }
  rcsb_chem_comp_synonyms { name type }
  rcsb_chem_comp_descriptor { SMILES InChIKey }
  drugbank { drugbank_info { description drug_categories drug_groups } }
} }"""


class Molecule(object):
    def __init__(self, data):
        comp = data.get("chem_comp") or {}
        info = data.get("rcsb_chem_comp_info") or {}
        drug = ((data.get("drugbank") or {}).get("drugbank_info") or {})
        self.comp_id = comp.get("id", "")
        self.name = (comp.get("name") or "").strip()
        self.formula = comp.get("formula") or ""
        self.weight = comp.get("formula_weight")
        self.type = comp.get("type") or ""
        self.charge = comp.get("pdbx_formal_charge") or 0
        self.atoms = info.get("atom_count") or 0
        self.heavy_atoms = info.get("atom_count_heavy") or 0
        self.chiral_centres = info.get("atom_count_chiral") or 0
        self.aromatic_bonds = info.get("bond_count_aromatic") or 0
        self.synonyms = [s.get("name", "") for s in data.get("rcsb_chem_comp_synonyms") or []
                         if s.get("type") != "Systematic Name" and s.get("name")]
        self.smiles = (data.get("rcsb_chem_comp_descriptor") or {}).get("SMILES") or ""
        # DrugBank marks its sources as [L8983,A1234]; they mean nothing here
        text = re.sub(r"\s*\[[A-Z]\d+(?:\s*,\s*[A-Z]\d+)*\]", "",
                      drug.get("description") or "")
        self.description = " ".join(text.split())
        self.categories = drug.get("drug_categories") or []
        self.groups = drug.get("drug_groups") or []
        self.exact = False           # the typed name is one of its own names
        self.how = ""                # how the name led here, in words
        self.alternatives = []       # other Molecules the name could mean
        self.entries = None          # number of PDB structures that contain it
        self.nickname = ""           # the everyday name that led here ("glucose")
        self.strong = False          # matched by its official or everyday name
        self.source = "PDB"          # "PDB" dictionary or "PubChem"
        self.identifier = self.comp_id   # what `get` needs to find it again

    def title(self):
        """A readable name: 'Cholesterol' rather than 'CHOLESTEROL'."""
        if self.nickname:
            return self.nickname[:1].upper() + self.nickname[1:] if len(self.nickname) > 4 \
                else self.nickname.upper()
        name = self.name
        if len(name) > 36 or re.search(r"\d[,\)]|\(\d", name):     # a systematic name
            short = sorted((s for s in self.synonyms if 3 < len(s) <= 36
                            and not re.search(r"\d[,\)]|\(\d|~", s)), key=len)
            if short:
                name = short[0]
        return name.capitalize() if name.isupper() and len(name) > 4 else name

    def all_names(self):
        return [self.name] + self.synonyms

    def object_name(self):
        """What the loaded PyMOL object is called."""
        return self.comp_id

    def sdf_url(self):
        return SDF_URL % self.comp_id

    def kind(self):
        """A plain word for what sort of molecule this is."""
        text = " | ".join(self.categories).lower() + " | " + self.type.lower()
        for pattern, label in ((r"steroids?", "steroid"), (r"fatty acids?", "fatty acid"),
                               (r"(phospho|sphingo|glyco)?lipids?", "lipid"),
                               (r"(mono|di|oligo)?saccharides?|carbohydrates?", "sugar"),
                               (r"nucleotides?", "nucleotide"),
                               (r"nucleosides?", "nucleoside"),
                               (r"amino acids?|l-peptide linking", "amino acid"),
                               (r"vitamins?", "vitamin"), (r"hormones?", "hormone")):
            # "Anti-Inflammatory Agents, Non-Steroidal" is not a steroid
            if re.search(r"(?<![a-z\-])%s(?![a-z])" % pattern, text):
                return label
        if self.heavy_atoms == 1:
            return "single atom or ion"
        return "small molecule"


def get(identifier):
    """The molecule for a stored identifier: a PDB code (CLR) or CID12345."""
    if re.match(r"^CID\d+$", identifier or ""):
        from . import pubchem
        return pubchem.by_cid(int(identifier[3:]))
    found = by_ids([identifier])
    return found[0] if found else None


def _norm(text):
    return re.sub(r"[^a-z0-9+]", "", text.lower())


def looks_like_id(text):
    return bool(re.match(r"^[A-Za-z0-9]{1,3}$|^[A-Za-z][A-Za-z0-9]{4}$", text.strip()))


def _search(attribute, operator, value, rows=12):
    payload = {"query": {"type": "terminal", "service": "text_chem",
                         "parameters": {"attribute": attribute, "operator": operator,
                                        "value": value}},
               "return_type": "mol_definition",
               "request_options": {"paginate": {"start": 0, "rows": rows}}}
    data = http.request_json(SEARCH_URL, payload, ttl=7 * cache.DAY)
    if not data:                 # HTTP 204: nothing found
        return []
    return [row["identifier"] for row in data.get("result_set", [])]


def by_ids(comp_ids):
    """[Molecule] for the IDs that exist, in the order given."""
    ids = [i.upper() for i in comp_ids if i]
    if not ids:
        return []
    listing = "[%s]" % ", ".join('"%s"' % i for i in ids)
    data = http.request_json(GRAPHQL_URL, {"query": _FIELDS % listing}, ttl=7 * cache.DAY)
    found = {}
    for item in ((data or {}).get("data") or {}).get("chem_comps") or []:
        if item and item.get("chem_comp"):
            molecule = Molecule(item)
            found[molecule.comp_id] = molecule
    return [found[i] for i in ids if i in found]


def count_entries(comp_id):
    """How many PDB structures contain this molecule (None if unknown)."""
    payload = {"query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
        "operator": "exact_match", "value": comp_id}},
        "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 1}}}
    try:
        data = http.request_json(SEARCH_URL, payload, ttl=7 * cache.DAY)
    except http.NetworkError:
        return None
    return int(data.get("total_count", 0)) if data else 0


def lookup(name, allow_close=True, allow_id=True):
    """The molecule a name stands for, or None.

    ``molecule.exact`` is True when the name is one of the molecule's own names
    (or a well-known everyday name); otherwise it is only the closest match and
    ``molecule.alternatives`` lists the others."""
    typed = " ".join(name.split())
    key = typed.lower()
    wanted = _norm(typed)
    if not wanted:
        return None

    if key in COMMON:
        found = by_ids([COMMON[key]])
        if found:
            found[0].exact = True
            found[0].strong = True
            if _norm(found[0].name) != wanted:
                found[0].nickname = typed
            found[0].how = "'%s' is the everyday name of %s" % (typed, found[0].name)
            return found[0]

    ids = []
    for comp_id in _search("chem_comp.name", "exact_match", typed):
        ids.append(comp_id)
    for attribute in ("rcsb_chem_comp_synonyms.name", "chem_comp.name"):
        for comp_id in _search(attribute, "contains_phrase", typed):
            if comp_id not in ids:
                ids.append(comp_id)
    if allow_id and looks_like_id(typed) and typed.upper() not in ids:
        ids.append(typed.upper())
    candidates = by_ids(ids[:20])
    if not candidates:
        return None

    exact = [m for m in candidates if any(_norm(n) == wanted for n in m.all_names())]
    if exact:
        # the one whose official name it is first, then the simplest molecule
        exact.sort(key=lambda m: (_norm(m.name) != wanted, m.heavy_atoms))
        best = exact[0]
        best.exact = True
        best.strong = _norm(best.name) == wanted
        if _norm(best.name) != wanted:
            best.nickname = typed          # show the name the user knows it by
        best.how = "'%s' is a name of %s" % (typed, best.comp_id)
        best.alternatives = [m for m in exact[1:4]]
        return best
    by_id = [m for m in candidates if m.comp_id == typed.upper()]
    if allow_id and by_id and looks_like_id(typed):
        by_id[0].exact = True
        by_id[0].how = "%s is the PDB's code for %s" % (typed.upper(), by_id[0].name)
        return by_id[0]
    if not allow_close:
        return None

    def closeness(molecule):
        hits = [n for n in molecule.all_names()
                if re.search(r"(?<![a-z])%s(?![a-z])" % re.escape(key), n.lower())]
        return (0 if hits else 1, min([len(n) for n in hits] or [999]),
                molecule.heavy_atoms)
    candidates.sort(key=closeness)
    best = candidates[0]
    if closeness(best)[0] == 1:
        return None
    best.how = "no molecule is called exactly '%s'; the closest name is %s" % (
        typed, best.name)
    best.alternatives = candidates[1:4]
    return best
