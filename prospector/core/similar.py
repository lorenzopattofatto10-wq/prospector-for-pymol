"""`similar`: related structures by sequence (pure Python)."""
from __future__ import annotations

from . import rcsb

MIN_LENGTH = 20


class Family(object):
    """All hits that are the same protein (same UniProt entry, or same name
    and organism), represented by its best structure."""

    def __init__(self, key):
        self.key = key
        self.members = []      # [(identity, evalue, entity dict)]

    @property
    def best(self):
        return min(self.members, key=lambda m: (-m[0], m[2]["resolution"] or 99,
                                                m[2]["pdb_id"]))

    @property
    def identity(self):
        return self.best[0]


class SimilarResult(object):
    def __init__(self):
        self.total = 0           # matching chains-entities reported by RCSB
        self.families = []       # [Family], most similar first
        self.query_length = 0
        self.message = ""


def find_similar(sequence, identity=30, limit=10, exclude_pdb=None):
    result = SimilarResult()
    sequence = "".join(c for c in sequence.upper() if c.isalpha())
    result.query_length = len(sequence)
    if len(sequence) < MIN_LENGTH:
        result.message = ("The chain has only %d residues; a sequence search needs at "
                          "least %d." % (len(sequence), MIN_LENGTH))
        return result
    total, hits = rcsb.search_sequence(sequence, identity / 100.0)
    result.total = total
    if not hits:
        result.message = ("No other structure shares at least %d%% sequence identity. "
                          "Try a lower threshold, e.g. identity=20." % identity)
        return result
    details = rcsb.entities([entity_id for entity_id, _, _ in hits])
    families = {}
    for entity_id, hit_identity, evalue in hits:
        info = details.get(entity_id)
        if info is None:
            continue
        if info["accessions"]:
            key = info["accessions"][0]
        else:
            key = ((info.get("description") or "").lower(),
                   (info.get("organism") or "").lower())
        family = families.setdefault(key, Family(key))
        family.members.append((hit_identity, evalue, info))
    ordered = sorted(families.values(), key=lambda f: (-f.identity, -len(f.members)))
    # put the structure being asked about out of the "best" slot where possible
    if exclude_pdb:
        for family in ordered:
            others = [m for m in family.members if m[2]["pdb_id"] != exclude_pdb.upper()]
            if others:
                family.members = others + [m for m in family.members
                                           if m[2]["pdb_id"] == exclude_pdb.upper()]
    result.families = ordered[:limit]
    result.hidden = max(0, len(ordered) - limit)
    return result


def render(result, label):
    if result.message:
        return result.message
    out = ["PROTEINS RELATED TO %s (by sequence, %d residues)" % (label, result.query_length),
           "   %d matching chains in the PDB, in %d distinct proteins. Most similar first:"
           % (result.total, len(result.families) + result.hidden), ""]
    out.append(" #  Ident.  PDB   Res(A)  Organism            Protein")
    for index, family in enumerate(result.families, 1):
        identity, evalue, info = family.best
        resolution = "%.2f" % info["resolution"] if info["resolution"] is not None else "  - "
        organism = info["organism"]
        if len(organism) > 18:
            organism = organism[:16] + ".."
        description = " ".join(info["description"].split())
        if len(description) > 34:
            description = description[:32] + ".."
        out.append("%2d  %4.0f%%   %-4s  %6s  %-18s  %s" % (
            index, 100 * identity, info["pdb_id"], resolution, organism, description))
    if result.hidden:
        out.append("    (+%d more distinct proteins; raise limit= to see them)" % result.hidden)
    out.append("")
    out.append("Each row is one protein (a group of near-identical sequences).")
    out.append("Ident. = how much of its sequence is identical to your chain.")
    return "\n".join(out)
