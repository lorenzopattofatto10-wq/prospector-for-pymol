"""Common organism names -> scientific name + NCBI taxonomy id."""
from __future__ import annotations

import re

# (scientific name, taxonomy id, aliases).  The taxonomy ids are the ones
# used by both UniProt (organism_id) and RCSB (taxonomy_lineage.id).
_TABLE = [
    ("Homo sapiens", 9606, ["human", "homo sapiens", "h. sapiens", "hs"]),
    ("Mus musculus", 10090, ["mouse", "mus musculus", "murine"]),
    ("Rattus norvegicus", 10116, ["rat", "rattus norvegicus"]),
    ("Bos taurus", 9913, ["bovine", "cow", "bos taurus", "cattle"]),
    ("Sus scrofa", 9823, ["pig", "porcine", "sus scrofa"]),
    ("Gallus gallus", 9031, ["chicken", "hen", "gallus gallus"]),
    ("Equus caballus", 9796, ["horse", "equine", "equus caballus"]),
    ("Oryctolagus cuniculus", 9986, ["rabbit"]),
    ("Physeter macrocephalus", 9755, ["sperm whale", "whale"]),
    ("Danio rerio", 7955, ["zebrafish", "danio rerio"]),
    ("Xenopus laevis", 8355, ["xenopus", "frog", "xenopus laevis"]),
    ("Drosophila melanogaster", 7227, ["drosophila", "fruit fly", "fly",
                                       "drosophila melanogaster"]),
    ("Caenorhabditis elegans", 6239, ["c. elegans", "c.elegans", "worm", "elegans",
                                      "caenorhabditis elegans"]),
    ("Saccharomyces cerevisiae", 4932, ["yeast", "s. cerevisiae", "baker's yeast",
                                        "saccharomyces cerevisiae"]),
    ("Schizosaccharomyces pombe", 4896, ["fission yeast", "s. pombe", "pombe"]),
    ("Arabidopsis thaliana", 3702, ["arabidopsis", "arabidopsis thaliana", "thale cress"]),
    ("Escherichia coli", 562, ["e. coli", "e.coli", "ecoli", "e coli", "escherichia coli"]),
    ("Bacillus subtilis", 1423, ["b. subtilis", "bacillus subtilis"]),
    ("Mycobacterium tuberculosis", 1773, ["tuberculosis", "m. tuberculosis", "mtb"]),
    ("Staphylococcus aureus", 1280, ["s. aureus", "staph", "staphylococcus aureus"]),
    ("Pseudomonas aeruginosa", 287, ["p. aeruginosa", "pseudomonas"]),
    ("Thermus thermophilus", 274, ["t. thermophilus", "thermus thermophilus"]),
    ("Plasmodium falciparum", 5833, ["malaria", "plasmodium", "p. falciparum"]),
    ("Severe acute respiratory syndrome coronavirus 2", 2697049,
     ["sars-cov-2", "sars cov 2", "sarscov2", "sars-cov2", "covid", "covid-19", "sars2"]),
    ("Severe acute respiratory syndrome coronavirus", 694009, ["sars-cov", "sars-cov-1"]),
    ("Human immunodeficiency virus 1", 11676, ["hiv", "hiv-1", "hiv1"]),
    ("Influenza A virus", 11320, ["influenza", "influenza a", "flu"]),
    ("Aequorea victoria", 6100, ["jellyfish", "aequorea victoria"]),
]

_ALIAS = {}
for _sci, _tax, _aliases in _TABLE:
    for _alias in _aliases + [_sci.lower()]:
        _ALIAS[_alias] = (_sci, _tax)

# When the user gives no organism, results from these are still preferred.
HUMAN = 9606
MODEL_ORGANISMS = {10090, 10116, 9913, 9823, 9031, 7955, 7227, 6239, 4932, 3702, 562,
                   2697049, 11676}


class Organism(object):
    def __init__(self, name, taxid=None, typed=None):
        self.name = name        # scientific name (or the user's text if unknown)
        self.taxid = taxid      # int, or None when we do not know the organism
        self.typed = typed or name

    def __repr__(self):
        return "Organism(%r, %r)" % (self.name, self.taxid)


def lookup(text):
    """Turn user text ('human', 'e_coli', 'Homo sapiens', '9606') into an Organism."""
    clean = text.strip().strip("\"'").replace("_", " ").lower()
    if clean.isdigit():
        for sci, tax, _ in _TABLE:
            if tax == int(clean):
                return Organism(sci, tax, text)
        return Organism(clean, int(clean), text)
    if clean in _ALIAS:
        sci, tax = _ALIAS[clean]
        return Organism(sci, tax, text)
    # unknown organism: keep the text, capitalised like a species name
    return Organism(clean[:1].upper() + clean[1:], None, text)


def extract(text):
    """Find an organism mentioned inside a free-text query.

    Returns (remaining_text, Organism or None).  Only long/unambiguous aliases
    are considered, and the organism is only removed if something is left."""
    lowered = " " + text.lower() + " "
    for alias in sorted(_ALIAS, key=len, reverse=True):
        if len(alias) < 3 or alias in ("fly", "worm", "rat", "hen", "pig", "cow", "flu"):
            continue
        match = re.search(r"(?<![\w\-])%s(?![\w\-])" % re.escape(alias), lowered)
        if match:
            rest = (lowered[:match.start()] + " " + lowered[match.end():]).split()
            if rest:
                sci, tax = _ALIAS[alias]
                return " ".join(rest), Organism(sci, tax, alias)
    return text, None


def preference(taxid=None, name=None):
    """0..1: how much a structure's organism is favoured when no filter is set."""
    if taxid == HUMAN or (name or "").lower() == "homo sapiens":
        return 1.0
    if taxid in MODEL_ORGANISMS:
        return 0.7
    if name:
        for sci, tax, _ in _TABLE:
            if sci.lower() == name.lower():
                return 0.7 if tax in MODEL_ORGANISMS else 0.5
    return 0.4
