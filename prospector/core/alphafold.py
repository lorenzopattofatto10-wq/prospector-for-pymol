"""AlphaFold DB client (https://alphafold.ebi.ac.uk/api).

GET /prediction/{uniprot accession} returns a LIST of models.  File URLs and
the model version are read from the answer, never hard-coded.
"""
from __future__ import annotations

from . import cache, http

API = "https://alphafold.ebi.ac.uk/api/prediction/%s"

# AlphaFold's standard confidence colours (pLDDT is stored in the B-factor column).
PLDDT_BANDS = [
    ("very high", 90, 100, "0x0053D6"),
    ("confident", 70, 90, "0x65CBF3"),
    ("low", 50, 70, "0xFFDB13"),
    ("very low", 0, 50, "0xFF7D45"),
]


class Model(object):
    def __init__(self, record):
        self.accession = record.get("uniprotAccession", "")
        self.model_id = record.get("modelEntityId") or record.get("entryId", "")
        self.name = record.get("uniprotDescription", "")
        self.gene = record.get("gene", "")
        self.organism = record.get("organismScientificName", "")
        self.cif_url = record.get("cifUrl", "")
        self.pdb_url = record.get("pdbUrl", "")
        self.version = record.get("latestVersion")
        self.mean_plddt = record.get("globalMetricValue")
        self.start = record.get("sequenceStart") or record.get("uniprotStart")
        self.end = record.get("sequenceEnd") or record.get("uniprotEnd")
        # fraction of residues in each confidence band (0..1), may be None
        self.fractions = [
            ("very high (>90)", record.get("fractionPlddtVeryHigh")),
            ("confident (70-90)", record.get("fractionPlddtConfident")),
            ("low (50-70)", record.get("fractionPlddtLow")),
            ("very low (<50)", record.get("fractionPlddtVeryLow")),
        ]

    def confidence_text(self):
        parts = ["%s %.0f%%" % (label, 100 * value) for label, value in self.fractions
                 if value is not None]
        return ", ".join(parts)

    def object_name(self):
        return self.model_id.replace("-", "_") or "AF_%s" % self.accession


def prediction(accession):
    """The canonical AlphaFold model for a UniProt accession, or None."""
    accession = accession.upper()
    try:
        data = http.request_json(API % accession, ttl=30 * cache.DAY)
    except http.NetworkError as err:
        if err.status in (400, 404, 422):
            return None
        raise
    if not data:
        return None
    records = [r for r in data if isinstance(r, dict)]
    # prefer the model of the canonical sequence (not isoforms / complexes)
    exact = [r for r in records if r.get("uniprotAccession", "").upper() == accession
             and not r.get("isComplex")]
    chosen = (exact or records)
    return Model(chosen[0]) if chosen else None
