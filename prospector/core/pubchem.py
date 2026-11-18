"""Small molecules that are not in the PDB's dictionary: PubChem (NIH).

PubChem holds over a hundred million compounds; for small, not too flexible
ones it offers a computed 3D shape.  Prospector asks it only when the PDB's
own dictionary (core/chemcomp.py) does not know a name.

PUG REST, as verified against the live service (its documentation page could
not be read from the build machine):
    /compound/name/<name>/property/<list>/JSON   exact-name lookup; 404 with
                                                  PUGREST.NotFound if unknown
    /compound/cid/<cid>/description/JSON         one or more text descriptions
    /compound/cid/<cid>/SDF?record_type=3d       the 3D file; 404 if there is none
    503 PUGREST.ServerBusy                       too many requests: try later
Proteins, large lipids and metal complexes have no 3D record (ConformerCount3D
is 0); those are reported, never invented.
"""
from __future__ import annotations

import re
import urllib.parse

from . import cache, chemcomp, http

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PROPERTIES = ("Title,MolecularFormula,MolecularWeight,IUPACName,Charge,HeavyAtomCount,"
              "ConformerCount3D,XLogP,HBondDonorCount,HBondAcceptorCount")
_RISKY_NAMES = set("all none water solvent organic inorganic polymer hydro metals backbone "
                   "sidechain visible enabled center origin".split())


class Compound(chemcomp.Molecule):
    """A PubChem compound, shaped like a dictionary molecule so the rest of
    Prospector can treat both the same way."""

    def __init__(self, record):
        chemcomp.Molecule.__init__(self, {})
        self.source = "PubChem"
        self.cid = int(record.get("CID"))
        self.comp_id = ""                      # it has no PDB code
        self.name = record.get("Title") or "CID %d" % self.cid
        formula = record.get("MolecularFormula") or ""
        self.formula = " ".join(re.findall(r"[A-Z][a-z]?\d*", formula)) or formula
        try:
            self.weight = float(record.get("MolecularWeight"))
        except (TypeError, ValueError):
            self.weight = None
        self.charge = record.get("Charge") or 0
        self.heavy_atoms = record.get("HeavyAtomCount") or 0
        self.conformers = record.get("ConformerCount3D") or 0
        self.iupac = record.get("IUPACName") or ""
        self.xlogp = record.get("XLogP")
        self.donors = record.get("HBondDonorCount")
        self.acceptors = record.get("HBondAcceptorCount")
        self.identifier = "CID%d" % self.cid

    def object_name(self):
        name = re.sub(r"[^A-Za-z0-9]+", "_", self.title()).strip("_").lower()[:24]
        if not name or name[0].isdigit() or name in _RISKY_NAMES:
            name = "m_" + (name or str(self.cid))
        return name

    def sdf_url(self):
        return "%s/compound/cid/%d/SDF?record_type=3d" % (BASE, self.cid)

    def has_shape(self):
        return self.conformers > 0


def _get(path):
    """JSON for a PUG REST path, or None when PubChem says 'not found'."""
    try:
        return http.request_json("%s/%s" % (BASE, path), ttl=7 * cache.DAY)
    except http.NetworkError as error:
        if error.status == 404:
            return None
        raise


def _describe(compound):
    data = _get("compound/cid/%d/description/JSON" % compound.cid)
    for item in ((data or {}).get("InformationList") or {}).get("Information") or []:
        if item.get("Description"):
            compound.description = " ".join(item["Description"].split())
            compound.description_source = item.get("DescriptionSourceName") or ""
            break
    return compound


def lookup(name):
    """The compound PubChem files under exactly this name, or None.  When a
    name has several compounds, the first one with a 3D shape is taken."""
    typed = " ".join(name.split())
    if not typed:
        return None
    data = _get("compound/name/%s/property/%s/JSON"
                % (urllib.parse.quote(typed, safe=""), PROPERTIES))
    records = ((data or {}).get("PropertyTable") or {}).get("Properties") or []
    if not records:
        return None
    compounds = [Compound(record) for record in records[:8]]
    shaped = [c for c in compounds if c.has_shape()]
    best = (shaped or compounds)[0]
    best.exact = True
    best.how = "'%s' is filed in PubChem as %s" % (typed, best.name)
    if best.name.lower() != typed.lower():
        best.nickname = typed
    return _describe(best)


def by_cid(cid):
    data = _get("compound/cid/%d/property/%s/JSON" % (int(cid), PROPERTIES))
    records = ((data or {}).get("PropertyTable") or {}).get("Properties") or []
    return _describe(Compound(records[0])) if records else None
