"""UniProt annotations mapped onto a structure (pure Python, no PyMOL)."""
from __future__ import annotations

import re

from . import alphafold, rcsb, sifts, uniprot

# key, title, UniProt feature types, PyMOL colour, how it is drawn
GROUPS = [
    ("active", "Active site", ["Active site"], "red", "sticks"),
    ("binding", "Binding site", ["Binding site"], "orange", "sticks"),
    ("site", "Site", ["Site"], "salmon", "sticks"),
    ("domain", "Domain", ["Domain", "DNA binding", "Zinc finger"], None, "cartoon"),
    ("disulfide", "Disulfide", ["Disulfide bond"], "yellow", "sticks"),
    ("ptm", "PTM", ["Modified residue", "Glycosylation", "Lipidation", "Cross-link"],
     "magenta", "sticks"),
    ("variant", "Disease variant", ["Natural variant"], "purple", "spheres"),
]
GROUP_KEYS = [g[0] for g in GROUPS]
# "site" (cleavage sites etc.) is noisy, so it is only shown when asked for.
DEFAULT_KEYS = [k for k in GROUP_KEYS if k != "site"]
_ALIASES = {"active_site": "active", "activesite": "active", "binding_site": "binding",
            "domains": "domain", "disulfides": "disulfide", "disulphide": "disulfide",
            "ss": "disulfide", "ptms": "ptm", "variants": "variant", "disease": "variant",
            "sites": "site"}
DOMAIN_COLORS = ["lightblue", "palegreen", "wheat", "lightpink", "palecyan", "lightorange",
                 "bluewhite", "limon"]
PAIR_TYPES = ("Disulfide bond", "Cross-link")   # start/end are two bonded residues

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U",
    "PYL": "O", "MSE": "M",
}


class FeatureError(ValueError):
    pass


def parse_types(text):
    """'active,binding' -> ['active', 'binding'];  '' -> all groups."""
    if not text or not text.strip():
        return list(DEFAULT_KEYS)
    if text.strip().lower() in ("all", "*"):
        return list(GROUP_KEYS)
    chosen = []
    for word in re.split(r"[,+\s]+", text.strip().lower()):
        if not word:
            continue
        word = _ALIASES.get(word, word)
        if word not in GROUP_KEYS:
            raise FeatureError("Unknown feature type '%s'. Choose from: %s"
                               % (word, ", ".join(GROUP_KEYS)))
        if word not in chosen:
            chosen.append(word)
    return chosen


class Feature(object):
    def __init__(self, group, kind, start, end, description, accession):
        self.group = group            # key from GROUPS
        self.kind = kind              # UniProt feature type
        self.start = start
        self.end = end
        self.description = description
        self.accession = accession
        self.pair = kind in PAIR_TYPES and start != end
        self.placed = {}              # {author chain: [residue labels]}

    def positions(self):
        if self.pair:
            return [self.start, self.end]
        return list(range(self.start, self.end + 1))

    def uniprot_text(self):
        if self.start == self.end:
            return str(self.start)
        return "%d%s%d" % (self.start, "<>" if self.pair else "-", self.end)

    def pdb_text(self):
        if not self.placed:
            return "not in this structure"
        if self.pair:
            both = sorted(c for c, labels in self.placed.items() if len(labels) == 2)
            if both:
                text = "%s/%s+%s" % (both[0], self.placed[both[0]][0], self.placed[both[0]][1])
                if len(both) > 1:
                    text += " (+%d)" % (len(both) - 1)
                return text
            singles = ["%s/%s" % (c, self.placed[c][0]) for c in sorted(self.placed)]
            return " + ".join(singles[:2]) if len(singles) > 1 else singles[0] + " (partner absent)"
        parts = []
        for chain in sorted(self.placed):
            labels = self.placed[chain]
            if len(labels) == 1:
                parts.append("%s/%s" % (chain, labels[0]))
            elif self.pair:
                parts.append("%s/%s" % (chain, "+".join(labels)))
            else:
                parts.append("%s/%s-%s" % (chain, labels[0], labels[-1]))
        text = parts[0]
        if len(parts) > 1:
            text += " (+%d chain%s)" % (len(parts) - 1, "s" if len(parts) > 2 else "")
        expected = len(self.positions())
        found = max(len(v) for v in self.placed.values())
        if found < expected:
            text += " [%d of %d res.]" % (found, expected)
        return text


def _is_disease(description):
    """A 'Natural variant' counts as a disease variant when UniProt names a
    condition:  'in LFS; germline mutation'  (but not just 'in dbSNP:rs123')."""
    text = re.sub(r"dbSNP:rs\d+", "", description or "")
    return bool(re.search(r"\bin [A-Za-z0-9]", text))


def extract(record, accession, wanted=None):
    """Read the features of a raw UniProt record."""
    wanted = wanted or GROUP_KEYS
    by_type = {}
    for key, _, kinds, _, _ in GROUPS:
        for kind in kinds:
            by_type[kind] = key
    found = []
    for item in record.get("features") or []:
        kind = item.get("type")
        group = by_type.get(kind)
        if group is None or group not in wanted:
            continue
        location = item.get("location") or {}
        start = (location.get("start") or {}).get("value")
        end = (location.get("end") or {}).get("value")
        if start is None and end is None:
            continue
        start = end if start is None else start
        end = start if end is None else end
        description = item.get("description") or ""
        ligand = (item.get("ligand") or {}).get("name")
        if ligand:
            description = ligand + (" - " + description if description else "")
        if kind == "Natural variant":
            if not _is_disease(description):
                continue
            alt = item.get("alternativeSequence") or {}
            change = "%s>%s" % (alt.get("originalSequence", "?"),
                                "/".join(alt.get("alternativeSequences") or ["del"]))
            description = "%s  %s" % (change, description)
        found.append(Feature(group, kind, int(start), int(end), description, accession))
    order = dict((key, i) for i, key in enumerate(GROUP_KEYS))
    found.sort(key=lambda f: (order[f.group], f.start, f.end))
    return found


class ProteinFeatures(object):
    def __init__(self, protein, record, chains, features):
        self.protein = protein        # uniprot.Protein
        self.sequence = ((record or {}).get("sequence") or {}).get("value", "")
        self.chains = chains          # [sifts.ChainMap]
        self.features = features      # [Feature], already placed

    def expected_residue(self, position):
        if 1 <= position <= len(self.sequence):
            return self.sequence[position - 1]
        return ""


class Annotation(object):
    def __init__(self, label, source):
        self.label = label            # "1PPB" or "AF-P04637-F1"
        self.source = source          # where the mapping came from
        self.proteins = []            # [ProteinFeatures]
        self.wanted = list(GROUP_KEYS)

    def selections(self):
        """What PyMOL should create: [(suffix, group key, title, {chain: [labels]},
        colour, style)].  Domains get one selection each."""
        plan = []
        domain_index = 0
        used = set()
        for key, title, _, colour, style in GROUPS:
            if key not in self.wanted:
                continue
            if key == "domain":
                for protein in self.proteins:
                    for feature in protein.features:
                        if feature.group != "domain" or not feature.placed:
                            continue
                        named = feature.description or feature.kind
                        base = re.sub(r"[^A-Za-z0-9]+", "_", named).strip("_")
                        suffix = "dom_%s" % (base[:20] or "domain")
                        number = 2
                        while suffix in used:
                            suffix = "dom_%s_%d" % (base[:20] or "domain", number)
                            number += 1
                        used.add(suffix)
                        plan.append((suffix, key, ("%s: %s" % (feature.kind,
                                                                feature.description)
                                                   ).rstrip(": "),
                                     dict(feature.placed),
                                     DOMAIN_COLORS[domain_index % len(DOMAIN_COLORS)], style))
                        domain_index += 1
                continue
            merged = {}
            for protein in self.proteins:
                for feature in protein.features:
                    if feature.group != key:
                        continue
                    for chain, labels in feature.placed.items():
                        bucket = merged.setdefault(chain, [])
                        for label in labels:
                            if label not in bucket:
                                bucket.append(label)
            if merged:
                plan.append((key, key, title, merged, colour, style))
        return plan

    def checks(self):
        """[(chain, label, expected one-letter code, UniProt position)] for every
        single mapped residue - lets the PyMOL layer verify the mapping against
        the atoms that are really there."""
        seen, out = set(), []
        for protein in self.proteins:
            for feature in protein.features:
                if feature.group == "domain":
                    continue
                for chain_map in protein.chains:
                    for position in feature.positions():
                        label = chain_map.resi(position)
                        if label is None or (chain_map.chain, label) in seen:
                            continue
                        seen.add((chain_map.chain, label))
                        out.append((chain_map.chain, label,
                                    protein.expected_residue(position), position))
        return out


def _place(features, chains):
    for feature in features:
        for chain_map in chains:
            labels = [chain_map.resi(p) for p in feature.positions()]
            labels = [label for label in labels if label is not None]
            if labels:
                feature.placed[chain_map.chain] = labels


def annotate(pdb_id, types=None):
    """Features for an experimental structure.  Raises FeatureError with a
    helpful message when it cannot be done."""
    wanted = types or list(DEFAULT_KEYS)
    detail = rcsb.entry_detail(pdb_id)
    if not detail:
        raise FeatureError("'%s' is not a PDB entry, so I cannot look up its features."
                           % pdb_id)
    mapping = sifts.mapping(pdb_id, detail)
    if not mapping.chains:
        raise FeatureError("%s has no chain linked to UniProt (it may be DNA/RNA, a "
                           "designed protein or an antibody), so there are no UniProt "
                           "features to map." % pdb_id.upper())
    annotation = Annotation(pdb_id.upper(), mapping.source)
    annotation.wanted = wanted
    for accession in mapping.accessions()[:6]:
        protein, record = uniprot.details(accession)
        if protein is None:
            continue
        chains = mapping.for_accession(accession)
        features = extract(record, accession, wanted)
        _place(features, chains)
        annotation.proteins.append(ProteinFeatures(protein, record, chains, features))
    return annotation


def annotate_model(accession, types=None):
    """Features for an AlphaFold model (its numbering is UniProt numbering)."""
    wanted = types or list(DEFAULT_KEYS)
    protein, record = uniprot.details(accession)
    if protein is None:
        raise FeatureError("UniProt has no entry %s." % accession)
    model = alphafold.prediction(accession)
    start = model.start if model else 1
    end = model.end if model else protein.length
    mapping = sifts.identity(protein.accession, start, end)
    annotation = Annotation(model.model_id if model else "AF-" + accession, mapping.source)
    annotation.wanted = wanted
    features = extract(record, protein.accession, wanted)
    _place(features, mapping.chains)
    annotation.proteins.append(ProteinFeatures(protein, record, mapping.chains, features))
    return annotation


def resi_text(labels):
    """PyMOL residue list: 57+102+195, with negative numbers escaped."""
    return "+".join(("\\" + label) if label.startswith("-") else label for label in labels)


def selection_text(obj, placed):
    """A PyMOL selection expression for {chain: [labels]} inside object ``obj``."""
    parts = []
    for chain in sorted(placed):
        if placed[chain]:
            parts.append("(chain %s and resi %s)" % (chain, resi_text(placed[chain])))
    if not parts:
        return "none"
    return "%s and (%s)" % (obj, " or ".join(parts))


ROW_LIMIT = {"active": 20, "binding": 20, "domain": 20, "disulfide": 20,
             "site": 8, "ptm": 8, "variant": 8}


def render(annotation, obj=None):
    """The table printed by `features`: what is in the structure, in detail;
    what is outside it, as a one-line summary."""
    out = []
    obj = obj or annotation.label
    titles = dict((g[0], g[1]) for g in GROUPS)
    for protein in annotation.proteins:
        info = protein.protein
        out.append("%s  <-  UniProt %s: %s%s, %s (%d aa)"
                   % (obj, info.accession, info.name,
                      " (gene %s)" % info.gene if info.gene else "", info.organism,
                      info.length))
        for chain_map in protein.chains[:4]:
            out.append("  chain %s: %s" % (chain_map.chain, chain_map.describe_numbering()))
        if len(protein.chains) > 4:
            out.append("  (+%d more chains with the same protein)" % (len(protein.chains) - 4))
        if not protein.features:
            out.append("  UniProt lists no features of the requested types for this protein.")
            out.append("")
            continue
        placed = [f for f in protein.features if f.placed]
        if placed:
            out.append("  %-16s %-10s %-22s %s" % ("Type", "UniProt", "PDB chain/residue",
                                                   "Description"))
        shown = {}
        for feature in placed:
            shown[feature.group] = shown.get(feature.group, 0) + 1
            if shown[feature.group] > ROW_LIMIT[feature.group]:
                continue
            title = titles[feature.group]
            if feature.group in ("domain", "ptm"):
                title = feature.kind if feature.kind != "Modified residue" else "Modified res."
            description = " ".join(feature.description.split())
            if len(description) > 60:
                description = description[:58] + ".."
            out.append("  %-16s %-10s %-22s %s" % (title, feature.uniprot_text(),
                                                   feature.pdb_text(), description))
        for key in GROUP_KEYS:
            extra = shown.get(key, 0) - ROW_LIMIT[key]
            if extra > 0:
                out.append("  ... +%d more %s rows (all are in the selection)"
                           % (extra, titles[key].lower()))
        missing = {}
        for feature in protein.features:
            if not feature.placed:
                missing[feature.group] = missing.get(feature.group, 0) + 1
        if missing:
            out.append("  Not in this structure (outside the crystallised/modelled part): "
                       + ", ".join("%d %s" % (missing[k], titles[k].lower())
                                   for k in GROUP_KEYS if k in missing))
        out.append("")
    return "\n".join(out).rstrip()
