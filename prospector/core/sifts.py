"""UniProt residue number -> PDB chain + author residue number.

Two steps, because author numbering is often irregular (offsets, gaps,
insertion codes such as 60A, 60B in thrombin):

  1. UniProt position -> position in the deposited sequence ("seq id").
     Source: SIFTS segments from PDBe
       GET https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{pdb id}
     PDBe is sometimes slow, so there is a short timeout and a fallback to the
     same SIFTS-derived alignment as served by RCSB
       (GraphQL polymer_entities.rcsb_polymer_entity_align).
  2. seq id -> the residue number printed in the file / shown by PyMOL.
     Source: RCSB GraphQL ...auth_to_entity_poly_seq_mapping (one label per
     sequence position, including insertion codes).
"""
from __future__ import annotations

from . import cache, http, rcsb

PDBE_SEGMENTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/%s"


class ChainMap(object):
    """Mapping for one chain of one UniProt protein."""

    def __init__(self, accession, chain, entity_id, numbering):
        self.accession = accession
        self.chain = chain              # author chain id, as PyMOL shows it
        self.entity_id = str(entity_id)
        self.numbering = numbering      # numbering[seq id - 1] = author label, e.g. "60A"
        self.segments = []              # [(unp_start, unp_end, seq_start)]

    def seq_id(self, position):
        for unp_start, unp_end, seq_start in self.segments:
            if unp_start <= position <= unp_end:
                return seq_start + position - unp_start
        return None

    def resi(self, position):
        """Author residue label for a UniProt position, or None if that residue
        is not part of this chain's construct."""
        seq = self.seq_id(position)
        if seq is None or not 1 <= seq <= len(self.numbering):
            return None
        return self.numbering[seq - 1]

    def uniprot_position(self, label):
        """The reverse: UniProt position of the residue labelled ``label`` in the
        file (e.g. '57' or '60A'), or None."""
        try:
            seq = self.numbering.index(str(label)) + 1
        except ValueError:
            return None
        for unp_start, unp_end, seq_start in self.segments:
            if seq_start <= seq <= seq_start + (unp_end - unp_start):
                return unp_start + seq - seq_start
        return None

    def resi_range(self, start, end):
        """All author labels for UniProt positions start..end that exist here."""
        labels = []
        for position in range(start, end + 1):
            label = self.resi(position)
            if label is not None:
                labels.append(label)
        return labels

    def covered(self):
        return sum(end - start + 1 for start, end, _ in self.segments)

    def describe_numbering(self):
        """A sentence a student can use: how PDB numbers relate to UniProt numbers."""
        offsets = set()
        example = None
        irregular = False
        for unp_start, unp_end, _ in self.segments:
            for position in range(unp_start, unp_end + 1):
                label = self.resi(position)
                if label is None:
                    continue
                if example is None:
                    example = (position, label)
                try:
                    offsets.add(int(label) - position)
                except ValueError:
                    irregular = True      # insertion code
                if len(offsets) > 3:
                    break
        if example is None:
            return "no residues mapped"
        if not irregular and offsets == {0}:
            return "same numbers as UniProt"
        if not irregular and len(offsets) == 1:
            offset = offsets.pop()
            return "PDB number = UniProt number %s %d" % ("+" if offset > 0 else "-",
                                                          abs(offset))
        return ("irregular numbering (gaps or insertion codes) - e.g. UniProt %d is "
                "residue %s" % example)


class Mapping(object):
    def __init__(self, pdb_id, source):
        self.pdb_id = pdb_id
        self.source = source      # "PDBe SIFTS" or "RCSB (SIFTS-derived)"
        self.chains = []          # [ChainMap]

    def accessions(self):
        seen = []
        for chain in self.chains:
            if chain.accession not in seen:
                seen.append(chain.accession)
        return seen

    def for_accession(self, accession):
        return [c for c in self.chains if c.accession == accession]


def _numbering(detail):
    """{author chain id: (entity id, [author label per seq position])}."""
    result = {}
    for entity in detail.get("polymer_entities") or []:
        entity_id = (entity.get("rcsb_polymer_entity_container_identifiers") or {}).get(
            "entity_id")
        for instance in entity.get("polymer_entity_instances") or []:
            ids = instance.get("rcsb_polymer_entity_instance_container_identifiers") or {}
            chain = ids.get("auth_asym_id")
            labels = ids.get("auth_to_entity_poly_seq_mapping") or []
            if chain:
                result[chain] = (str(entity_id), [str(label) for label in labels])
    return result


def _segments_from_pdbe(pdb_id):
    """{(accession, author chain): [(unp_start, unp_end, seq_start)]} or None."""
    try:
        data = http.request_json(PDBE_SEGMENTS % pdb_id.lower(), ttl=30 * cache.DAY,
                                 timeout=8, retry=False)
    except http.NetworkError:
        return None
    block = ((data or {}).get(pdb_id.lower()) or {}).get("UniProt") or {}
    segments = {}
    for accession, info in block.items():
        for item in info.get("mappings") or []:
            try:
                key = (accession, item["chain_id"])
                segment = (int(item["unp_start"]), int(item["unp_end"]),
                           int(item["start"]["residue_number"]))
            except (KeyError, TypeError, ValueError):
                continue
            segments.setdefault(key, []).append(segment)
    return segments or None


def _segments_from_rcsb(detail):
    segments = {}
    for entity in detail.get("polymer_entities") or []:
        chains = (entity.get("rcsb_polymer_entity_container_identifiers") or {}).get(
            "auth_asym_ids") or []
        for align in entity.get("rcsb_polymer_entity_align") or []:
            if (align.get("reference_database_name") or "").lower() != "uniprot":
                continue
            accession = align.get("reference_database_accession")
            for region in align.get("aligned_regions") or []:
                start, seq, length = (region.get("ref_beg_seq_id"),
                                      region.get("entity_beg_seq_id"), region.get("length"))
                if not (start and seq and length):
                    continue
                for chain in chains:
                    segments.setdefault((accession, chain), []).append(
                        (start, start + length - 1, seq))
    return segments


def mapping(pdb_id, detail=None, prefer="pdbe"):
    """Build the residue mapping for an entry.  Returns a Mapping (possibly with
    no chains, when the entry has no UniProt cross-reference) or None if the
    entry does not exist."""
    detail = detail or rcsb.entry_detail(pdb_id)
    if not detail:
        return None
    numbering = _numbering(detail)
    segments, source = None, ""
    if prefer == "pdbe":
        segments, source = _segments_from_pdbe(pdb_id), "PDBe SIFTS"
    if not segments:
        segments, source = _segments_from_rcsb(detail), "RCSB (SIFTS-derived alignment)"
    result = Mapping(pdb_id.upper(), source)
    for (accession, chain), parts in sorted(segments.items()):
        if chain not in numbering:
            continue
        entity_id, labels = numbering[chain]
        chain_map = ChainMap(accession, chain, entity_id, labels)
        chain_map.segments = sorted(set(parts))
        result.chains.append(chain_map)
    return result


def identity(accession, start, end, chain="A"):
    """Mapping for an AlphaFold model: residue numbers ARE UniProt numbers."""
    result = Mapping("AF-" + accession, "AlphaFold model (UniProt numbering)")
    labels = [str(i) for i in range(1, (end or 0) + 1)]
    chain_map = ChainMap(accession, chain, "1", labels)
    chain_map.segments = [(start or 1, end or len(labels), start or 1)]
    result.chains.append(chain_map)
    return result
