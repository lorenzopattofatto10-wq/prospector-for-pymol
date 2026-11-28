"""Hands-on tools: scout, explain, measure, mutate, sequence, back, session, report."""
from __future__ import annotations

import base64
import html
import json
import os
import re
import time

from pymol import cmd

from ..core import cache, card, http, rcsb, report
from ..core import features as feat
from . import commands as base
from . import places, teach, worker

_say = base._say

# ------------------------------------------------------------ biochemistry
AA = {   # three-letter code: (name, group, what it is good for)
    "ALA": ("alanine", "small, water-avoiding", "a tiny methyl side chain; fits anywhere, "
            "common in helices"),
    "ARG": ("arginine", "positively charged", "a long arm ending in a + charge; grips "
            "phosphates, DNA and negatively charged partners"),
    "ASN": ("asparagine", "polar", "an amide that makes hydrogen bonds; the residue sugars "
            "are attached to (N-glycosylation)"),
    "ASP": ("aspartate", "negatively charged", "a short arm with a - charge; binds metal "
            "ions and often does acid/base chemistry in active sites"),
    "CYS": ("cysteine", "sulfur-containing", "a reactive -SH group; pairs up into "
            "disulfide bonds, binds metals such as zinc, and can be a catalytic nucleophile"),
    "GLN": ("glutamine", "polar", "an amide on a longer arm; makes hydrogen bonds"),
    "GLU": ("glutamate", "negatively charged", "a longer arm with a - charge; binds metal "
            "ions and does acid/base chemistry"),
    "GLY": ("glycine", "tiny, flexible", "no side chain at all, so the backbone can bend "
            "sharply here; common in tight turns and loops"),
    "HIS": ("histidine", "polar, sometimes charged", "a ring (imidazole) that can pick up "
            "and release a proton near neutral pH; a favourite in active sites and for "
            "holding metals such as zinc, iron and copper"),
    "ILE": ("isoleucine", "water-avoiding (hydrophobic)", "a branched, oily side chain; "
            "packs into the protein core"),
    "LEU": ("leucine", "water-avoiding (hydrophobic)", "an oily side chain; the most "
            "common residue in the protein core and in helices"),
    "LYS": ("lysine", "positively charged", "a long flexible arm ending in a + charge; "
            "sits on the surface, binds DNA, and is where ubiquitin and acetyl groups "
            "are attached"),
    "MET": ("methionine", "water-avoiding, sulfur-containing", "a flexible oily side "
            "chain; also the first residue of almost every newly made protein"),
    "PHE": ("phenylalanine", "aromatic, water-avoiding", "a flat benzene ring; stacks in "
            "the core and against ligands"),
    "PRO": ("proline", "rigid", "its side chain loops back onto the backbone, making a "
            "stiff kink; breaks helices and shapes turns"),
    "SER": ("serine", "polar", "a small -OH group; makes hydrogen bonds, gets "
            "phosphorylated, and is the attacking group in serine proteases"),
    "THR": ("threonine", "polar", "an -OH group plus a methyl; makes hydrogen bonds and "
            "gets phosphorylated"),
    "TRP": ("tryptophan", "aromatic", "the largest side chain, a double ring; anchors "
            "proteins in membranes and stacks against ligands and sugars"),
    "TYR": ("tyrosine", "aromatic, polar", "a benzene ring with an -OH; stacks like "
            "phenylalanine but can also hydrogen-bond and be phosphorylated"),
    "VAL": ("valine", "water-avoiding (hydrophobic)", "a small branched oily side chain; "
            "common in beta-strands and the core"),
}
MAX_AREA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225,
            "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
            "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
            "TYR": 263, "VAL": 174}
STANDARD = "+".join(sorted(AA))
SS_WORDS = {"H": "an alpha-helix", "S": "a beta-strand", "L": "a loop", "": "a loop"}


def _split_object(words):
    """(object name or None, remaining words): the object may be named at the
    start or the end, by its PDB ID or by the name used to load it."""
    for size in range(len(words), 0, -1):
        if base._known_object(" ".join(words[:size])):
            return base._known_object(" ".join(words[:size])), words[size:]
    for size in range(len(words) - 1, 0, -1):
        if base._known_object(" ".join(words[-size:])):
            return base._known_object(" ".join(words[-size:])), words[:-size]
    return None, list(words)


def _run(label, words, action, is_name=None, need_annotation=False, snapshot=True,
         allow_fetch=True, extra_job=None, protein_only=True):
    """Shared plumbing: find the structure (loading it by name if needed), fetch
    its annotations off the GUI thread, then call ``action(obj, rest, extra)``.

    ``is_name(word)`` says whether a leftover word could be part of a protein
    name (rather than an argument of the command)."""
    state = {}

    def before():
        obj, rest = _split_object(words)
        query = ""
        if obj is None:
            name_words = [w for w in rest if is_name is None or is_name(w)]
            from . import typed
            if name_words and base._molecules() and typed.not_a_name(
                    " ".join(name_words), strict=True):
                # a keyword, a typo of one, a number ...: an argument for the
                # structure on screen, not a protein to go and load
                name_words = []
            if name_words and allow_fetch:
                query = " ".join(name_words)
                rest = [w for w in rest if w not in name_words]
            else:
                obj, message = base._pick_object("", label.split()[0])
                if obj is None:
                    base._header(label)
                    _say(message)
                    _say(report.LINE)
                    return False
        state.update(obj=obj, rest=rest, query=query)
        return True

    def job():
        extra = {"annotation": None, "names": {}}
        identity = None
        if state["query"]:
            state["fetched"] = base._fetch_query(state["query"], strict=True)
            if state["fetched"]["error"]:
                return extra
            identity = state["fetched"]["identity"]
        elif state["obj"]:
            identity = base._identify(state["obj"])
        if need_annotation:
            extra["annotation"], extra["names"] = places.fetch_annotation(identity)
        extra["identity"] = identity
        if extra_job is not None:      # more network work for this command
            extra_job(identity, extra)
        return extra

    def done(extra, error):
        if error:
            return base._unexpected(error)
        base._header(label)
        fetched = state.get("fetched")
        if fetched:
            if fetched["error"]:
                _say(fetched["error"])
                _say(report.LINE)
                return
            state["obj"] = base._place_fetched(fetched)
        if protein_only and base._needs_protein(state["obj"], label.split()[0]):
            return
        if snapshot:
            base._snapshot(label)
        base._last[0] = state["obj"]
        action(state["obj"], state["rest"], extra)
        _say(report.LINE)

    worker.submit(job, done, before)


# -------------------------------------------------------- residue addressing
_SPEC = re.compile(r"^(?:([A-Za-z0-9]{1,4})/)?(?:[A-Za-z]{3})?(-?\d+[A-Za-z]?)(?:/([A-Za-z0-9']{1,4}))?$")


def _is_spec(word):
    return bool(_SPEC.match(word))


def _spec_selection(obj, word):
    """'H/57', '57', 'A/HIS87/NE2' or a ligand code -> (selection, description) or None."""
    match = _SPEC.match(word)
    if match:
        chain, number, atom = match.groups()
        label = ("\\" + number) if number.startswith("-") else number
        selection = "%s and resi %s" % (obj, label)
        if chain:
            selection += " and chain %s" % chain
        if atom:
            selection += " and name %s" % atom
        # a number means the protein's residue of that number, not a water or
        # ion that happens to carry the same number
        if cmd.count_atoms(selection + " and polymer"):
            selection += " and polymer"
        if not chain and cmd.count_atoms(selection):
            # several chains have this number: take the first, so that one
            # residue is meant and not one in every chain
            found = []
            cmd.iterate("first (%s)" % selection, "found.append(chain)",
                        space={"found": found})
            if found and found[0]:
                selection += " and chain %s" % found[0]
        if cmd.count_atoms(selection):
            return selection
        return None
    selection = "%s and resn %s and not polymer" % (obj, word.upper())
    if re.match(r"^[A-Za-z0-9]{1,3}$", word) and cmd.count_atoms(selection):
        # the first copy only, so a distance means something
        found = []
        cmd.iterate("first (%s)" % selection, "found.append((chain, resi))",
                    space={"found": found})
        chain, number = found[0]
        return "%s and resn %s and resi %s%s" % (obj, word.upper(), number,
                                                 " and chain %s" % chain if chain else "")
    return None


def _name_of(selection):
    """'HIS57 (chain H)' for a one-residue selection."""
    found = []
    cmd.iterate("first (%s)" % selection, "found.append((resn, resi, chain))",
                space={"found": found})
    if not found:
        return "?"
    resn, resi, chain = found[0]
    return "%s%s%s" % (resn, resi, " (chain %s)" % chain if chain else "")


def _show_residue(selection, colour="orange"):
    cmd.show("sticks", "byres (%s)" % selection)
    cmd.color(colour, "byres (%s) and elem C" % selection)
    cmd.util.cnc("byres (%s)" % selection)


# -------------------------------------------------------------------- scout
def _class(title, residues, blurb):
    return (title, lambda o: "%s and polymer and resn %s" % (o, residues), blurb, "sticks")


SCOUT = {}       # word -> (title, selection builder, explanation, style)
for _code, (_name, _group, _use) in AA.items():
    _entry = ("%s residues" % _name,
              (lambda code: (lambda o: "%s and polymer and resn %s" % (o, code)))(_code),
              "%s (%s, %s): %s" % (_name, _code, _group, _use), "sticks")
    for _word in (_name, _name + "s", _code.lower()):
        SCOUT[_word] = _entry
SCOUT["histidines"] = SCOUT["histidine"]
SCOUT["cysteines"] = SCOUT["cysteine"]
_GROUPS = [
    (["hydrophobic", "hydrophobics", "oily", "nonpolar", "apolar"],
     _class("water-avoiding (hydrophobic) residues", "ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO",
            "oily side chains that hide from water; they pack together in the core and "
            "hold the fold together")),
    (["polar", "polars", "hydrophilic"],
     _class("polar residues", "SER+THR+ASN+GLN+TYR+CYS+HIS",
            "side chains that like water and make hydrogen bonds; common on the surface "
            "and in binding sites")),
    (["charged", "charges", "ionic"],
     _class("charged residues", "ASP+GLU+LYS+ARG",
            "side chains with a full + or - charge; nearly always on the surface, where "
            "they attract partners and form salt bridges")),
    (["positive", "basic", "cationic", "positives"],
     _class("positively charged residues", "LYS+ARG+HIS",
            "lysine, arginine (and histidine); they bind DNA, RNA, phosphates and "
            "negatively charged partners")),
    (["negative", "acidic", "anionic", "negatives"],
     _class("negatively charged residues", "ASP+GLU",
            "aspartate and glutamate; they bind metal ions and positively charged "
            "partners")),
    (["aromatic", "aromatics", "rings"],
     _class("aromatic residues", "PHE+TYR+TRP+HIS",
            "flat ring side chains; they stack on each other, on ligands and on DNA "
            "bases, and absorb UV light at 280 nm")),
    (["small", "tiny"],
     _class("small residues", "GLY+ALA+SER",
            "the smallest side chains; they allow tight packing and sharp turns")),
    (["helix", "helices", "alpha", "alpha-helix", "alpha-helices"],
     ("alpha-helices", lambda o: "%s and polymer and ss H" % o,
      "the backbone coiled like a spring, held by hydrogen bonds between residue i and "
      "i+4", "cartoon")),
    (["strand", "strands", "sheet", "sheets", "beta", "beta-sheet", "beta-sheets",
      "beta-strand", "beta-strands"],
     ("beta-strands", lambda o: "%s and polymer and ss S" % o,
      "stretched-out backbone segments that line up side by side into sheets, held by "
      "hydrogen bonds between neighbouring strands", "cartoon")),
    (["loop", "loops", "coil", "coils", "turns"],
     ("loops", lambda o: "%s and polymer and not ss H+S" % o,
      "the irregular stretches linking helices and strands; often flexible, and often "
      "where binding sites and active sites sit", "cartoon")),
    (["metal", "metals", "metal ions"],
     ("metal ions", lambda o: "%s and metals" % o,
      "single metal atoms such as zinc, iron, magnesium or calcium; they hold folds "
      "together or take part in catalysis", "balls")),
    (["ion", "ions"],
     ("ions", lambda o: "%s and inorganic" % o,
      "single charged atoms or small inorganic groups bound to the structure", "balls")),
    (["ligand", "ligands", "small molecules", "drug", "drugs", "cofactor", "cofactors"],
     ("bound small molecules", lambda o: "%s and organic" % o,
      "molecules that are not part of the protein chain: drugs, cofactors, substrates, "
      "and also additives from the crystallisation solution", "ligand")),
    (["water", "waters", "solvent"],
     ("water molecules", lambda o: "%s and solvent" % o,
      "ordered water molecules the experiment could see; they often bridge hydrogen "
      "bonds between protein and ligand", "balls")),
    (["sugar", "sugars", "glycan", "glycans", "carbohydrate", "carbohydrates"],
     ("sugars (glycans)",
      lambda o: "%s and resn NAG+NDG+MAN+BMA+FUC+FUL+GAL+GLC+BGC+SIA+XYP+GLA+A2G" % o,
      "sugar units attached to the protein (glycosylation) or bound to it", "ligand")),
    (["dna", "rna", "nucleic", "nucleic acid", "nucleic acids", "nucleotides"],
     ("DNA / RNA", lambda o: "%s and polymer.nucleic" % o,
      "nucleic acid chains in the structure", "cartoon")),
    (["adenine", "adenines"], ("adenine bases", lambda o: "%s and resn DA+A" % o,
                               "the purine base A; pairs with T (or U) by two hydrogen "
                               "bonds", "sticks")),
    (["guanine", "guanines"], ("guanine bases", lambda o: "%s and resn DG+G" % o,
                               "the purine base G; pairs with C by three hydrogen bonds",
                               "sticks")),
    (["cytosine", "cytosines"], ("cytosine bases", lambda o: "%s and resn DC+C" % o,
                                 "the pyrimidine base C; pairs with G by three hydrogen "
                                 "bonds", "sticks")),
    (["thymine", "thymines"], ("thymine bases", lambda o: "%s and resn DT" % o,
                               "the pyrimidine base T (DNA only); pairs with A", "sticks")),
    (["uracil", "uracils"], ("uracil bases", lambda o: "%s and resn U" % o,
                             "the pyrimidine base U (RNA only); pairs with A", "sticks")),
    (["modified", "modified residues", "nonstandard", "non-standard", "unusual"],
     ("modified (non-standard) residues",
      lambda o: "%s and polymer.protein and not resn %s" % (o, STANDARD),
      "residues in the chain that are not one of the 20 standard amino acids: "
      "phosphorylated, selenium-labelled, oxidised and so on", "sticks")),
    (["phospho", "phosphorylated", "phosphorylation", "phosphosites"],
     ("phosphorylated residues", lambda o: "%s and resn SEP+TPO+PTR" % o,
      "serine, threonine or tyrosine carrying a phosphate group - the cell's commonest "
      "on/off switch", "sticks")),
]
for _words, _entry in _GROUPS:
    for _word in _words:
        SCOUT[_word] = _entry
# things that have to be computed from the coordinates
COMPUTED = {
    "disulfide": "disulfides", "disulfides": "disulfides", "disulfide bonds": "disulfides",
    "disulphides": "disulfides", "ss bonds": "disulfides", "bridges": "disulfides",
    "free cysteines": "free_cys", "free cysteine": "free_cys", "thiols": "free_cys",
    "salt bridges": "salt", "salt bridge": "salt", "ion pairs": "salt",
    "surface": "surface", "exposed": "surface", "outside": "surface",
    "buried": "buried", "core": "buried", "inside": "buried",
    "flexible": "flexible", "mobile": "flexible", "floppy": "flexible",
    "rigid": "rigid", "stable": "rigid", "ordered": "rigid",
    "interface": "interface", "interfaces": "interface", "contacts": "interface",
    "n-terminus": "nterm", "n terminus": "nterm", "nterm": "nterm", "start": "nterm",
    "c-terminus": "cterm", "c terminus": "cterm", "cterm": "cterm", "end": "cterm",
    "termini": "termini", "ends": "termini",
    "glycosylation sites": "sequon", "sequons": "sequon", "sequon": "sequon",
}
_MOTIF = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx\[\]\^]{2,30}$")


def _scout_words():
    return set(SCOUT) | set(COMPUTED) | set(["motif", "sequence", "pattern"])


def _ensure_ss(obj):
    """Some files carry no helix/strand records; let PyMOL work them out then."""
    if cmd.count_atoms("%s and polymer and name CA" % obj) and \
            not cmd.count_atoms("%s and polymer and name CA and ss H+S" % obj):
        cmd.dss(obj)
        teach.note("dss %s" % obj, "assign helices and strands from the geometry")


def _per_residue(obj):
    """[(chain, resi, resn, ss, b)] for every protein residue, in order."""
    _ensure_ss(obj)
    rows = []
    cmd.iterate("%s and polymer and name CA" % obj,
                "rows.append((chain, resi, resn, ss, b))", space={"rows": rows})
    return rows


def _residue_areas(obj):
    """{(chain, resi): relative exposure 0..1} from solvent-accessible surface."""
    scratch = "_prospector_area"
    cmd.delete(scratch)
    cmd.create(scratch, "%s and polymer" % obj, zoom=0)   # zoom=0: keep the camera
    saved = cmd.get("dot_solvent")
    cmd.set("dot_solvent", 1)
    cmd.get_area(scratch, load_b=1)
    cmd.set("dot_solvent", saved)
    totals, names = {}, {}
    cmd.iterate(scratch, "totals[(chain, resi)] = totals.get((chain, resi), 0.0) + b; "
                         "names[(chain, resi)] = resn",
                space={"totals": totals, "names": names})
    cmd.delete(scratch)
    return dict((key, min(1.0, area / float(MAX_AREA.get(names[key].upper(), 200))))
                for key, area in totals.items())


def _keys_selection(obj, keys):
    """Selection for a list of (chain, resi)."""
    by_chain = {}
    for chain, resi in keys:
        by_chain.setdefault(chain, []).append(resi)
    if not by_chain:
        return "none"
    parts = []
    for chain, numbers in by_chain.items():
        text = feat.resi_text(numbers)
        parts.append("(chain %s and resi %s)" % (chain, text) if chain else "(resi %s)" % text)
    return "%s and polymer and (%s)" % (obj, " or ".join(parts))


def _motif_regex(pattern):
    """'CxxC' or 'N[^P][ST]' -> compiled regex over one-letter sequence."""
    text = pattern.upper().replace("X", ".")
    try:
        return re.compile(text)
    except re.error:
        return None


def _computed(obj, kind, extra_word=""):
    """(title, selection, explanation, style) for the coordinate-based scouts."""
    rows = _per_residue(obj)
    if kind in ("disulfides", "free_cys"):
        sulfurs = cmd.get_model("%s and polymer and resn CYS and name SG" % obj).atom
        bonded = set()
        for i, a in enumerate(sulfurs):
            for b in sulfurs[i + 1:]:
                d2 = sum((a.coord[k] - b.coord[k]) ** 2 for k in range(3))
                if d2 <= 2.5 ** 2:
                    bonded.add((a.chain, a.resi))
                    bonded.add((b.chain, b.resi))
        if kind == "disulfides":
            return ("cysteines in disulfide bonds (%d bonds)" % (len(bonded) // 2),
                    _keys_selection(obj, bonded),
                    "two cysteines whose sulfur atoms are bonded (S-S, about 2.0 A apart); "
                    "these bonds staple the fold together, mostly in proteins that work "
                    "outside the cell", "sticks")
        free = [(a.chain, a.resi) for a in sulfurs if (a.chain, a.resi) not in bonded]
        return ("free cysteines (not in a disulfide)", _keys_selection(obj, free),
                "cysteines with a free -SH group; they may bind a metal, be a catalytic "
                "residue, or simply sit in the core", "sticks")
    if kind == "salt":
        acid = "(%s and polymer and resn ASP+GLU and name OD1+OD2+OE1+OE2)" % obj
        basic = "(%s and polymer and resn LYS+ARG+HIS and name NZ+NH1+NH2+NE+ND1+NE2)" % obj
        selection = "byres ((%s within 4.0 of %s) or (%s within 4.0 of %s))" % (
            acid, basic, basic, acid)
        return ("residues in salt bridges", selection,
                "a negatively charged side chain (Asp, Glu) within 4 A of a positively "
                "charged one (Lys, Arg, His); the opposite charges attract and steady the "
                "fold", "sticks")
    if kind in ("surface", "buried"):
        areas = _residue_areas(obj)
        if kind == "surface":
            keys = [key for key, value in areas.items() if value >= 0.5]
            return ("exposed residues (over half of the side chain touches water)",
                    _keys_selection(obj, keys),
                    "residues on the outside of the protein; this is what other molecules "
                    "meet, and it is mostly polar and charged", "cartoon")
        keys = [key for key, value in areas.items() if value <= 0.05]
        return ("buried residues (hidden from water)", _keys_selection(obj, keys),
                "residues packed in the interior; mostly oily side chains - this is the "
                "hydrophobic core that drives folding", "cartoon")
    if kind in ("flexible", "rigid"):
        values = sorted(row[4] for row in rows)
        if not values:
            return None
        cut = max(1, len(values) // 7)
        if kind == "flexible":
            limit = values[-cut]
            keys = [(r[0], r[1]) for r in rows if r[4] >= limit]
            return ("most mobile residues (highest B-factors)", _keys_selection(obj, keys),
                    "the 15% of residues whose atoms were most smeared out in the "
                    "experiment: loops, chain ends and anything that moves", "cartoon")
        limit = values[cut - 1]
        keys = [(r[0], r[1]) for r in rows if r[4] <= limit]
        return ("most rigid residues (lowest B-factors)", _keys_selection(obj, keys),
                "the 15% of residues that were pinned down best in the experiment, "
                "usually the tightly packed core", "cartoon")
    if kind == "interface":
        chains = [c for c in cmd.get_chains("%s and polymer" % obj) if c]
        if len(chains) < 2:
            return ("chain interfaces", "none",
                    "this structure has only one chain, so there is no interface "
                    "between chains", "sticks")
        parts = ["(chain %s and byres ((%s and polymer and chain %s) within 4 of "
                 "(%s and polymer and not chain %s)))" % (c, obj, c, obj, c) for c in chains]
        return ("residues where chains touch each other",
                "%s and polymer and (%s)" % (obj, " or ".join(parts)),
                "residues within 4 A of another chain; these contacts hold the "
                "assembly together", "sticks")
    if kind in ("nterm", "cterm", "termini"):
        firsts, lasts = {}, {}
        for chain, resi, resn, ss, b in rows:
            firsts.setdefault(chain, (chain, resi))
            lasts[chain] = (chain, resi)
        keys = []
        if kind in ("nterm", "termini"):
            keys += list(firsts.values())
        if kind in ("cterm", "termini"):
            keys += list(lasts.values())
        what = {"nterm": "N-terminus (start) of each chain",
                "cterm": "C-terminus (end) of each chain",
                "termini": "both ends of each chain"}[kind]
        return (what, _keys_selection(obj, keys),
                "a protein chain is made from the N-terminus (free amino group) to the "
                "C-terminus (free carboxyl group); only the residues that were modelled "
                "are shown, so the true ends may be missing", "balls")
    if kind in ("sequon", "motif"):
        pattern = "N[^P][ST]" if kind == "sequon" else extra_word
        regex = _motif_regex(pattern)
        if regex is None:
            return None
        by_chain = {}
        for chain, resi, resn, ss, b in rows:
            by_chain.setdefault(chain, []).append((resi, feat.THREE_TO_ONE.get(resn.upper(), "X")))
        keys = []
        for chain, items in by_chain.items():
            sequence = "".join(letter for _, letter in items)
            for match in re.finditer("(?=(%s))" % regex.pattern, sequence):
                for offset in range(len(match.group(1))):
                    key = (chain, items[match.start() + offset][0])
                    if key not in keys:
                        keys.append(key)
        if kind == "sequon":
            return ("possible N-glycosylation sites (Asn-X-Ser/Thr)",
                    _keys_selection(obj, keys),
                    "the sequence pattern N-x-S/T (x is not proline); cells can attach a "
                    "sugar chain to the asparagine. A match only means it is possible",
                    "sticks")
        return ("the sequence motif %s" % pattern.upper(), _keys_selection(obj, keys),
                "every place in the amino-acid sequence that matches this pattern "
                "(x = any residue)", "sticks")
    return None


def _scout_menu():
    _say("")
    _say("WHAT SCOUT CAN SPOT   (use:  scout <protein> <what>)")
    _say("   Amino acids   any of the 20 by name: histidine, cysteine, tryptophan ...")
    _say("   By property   hydrophobic, polar, charged, positive, negative, aromatic,")
    _say("                 small")
    _say("   Shape         helices, strands, loops, termini, n-terminus, c-terminus")
    _say("   Bonds         disulfides, free cysteines, salt bridges")
    _say("   Location      surface, buried, interface (where chains touch),")
    _say("                 flexible, rigid")
    _say("   Not protein   ligands, metals, ions, sugars, water, dna, rna,")
    _say("                 adenine, guanine, cytosine, thymine, uracil")
    _say("   Chemistry     modified, phospho, glycosylation sites")
    _say("   Sequence      motif <pattern>, e.g.  scout p53 motif CxxC   (x = any)")
    _say("")
    _say("EXAMPLES     scout p53 histidine   |   scout haemoglobin salt bridges")
    _say("             scout thrombin disulfides   |   scout p53 surface")


def scout(*args, **kwargs):
    """
DESCRIPTION

    Spot something in a structure: highlight every histidine, every
    disulfide bond, the hydrophobic core, the metals, a sequence motif ...
    Type "scout" alone to see everything it can look for.

USAGE

    scout                         list what can be spotted
    scout <protein> <what>
    scout <what>                  in the loaded structure

EXAMPLES

    scout p53 histidine
    scout thrombin disulfides
    scout haemoglobin salt bridges
    scout p53 motif CxxC
    """
    words = base._raw(args, kwargs).replace(",", " ").split()
    label = "scout %s" % " ".join(words)
    if not words:
        worker.submit(lambda: None, lambda r, e: (base._header(label), _scout_menu(),
                                                  _say(report.LINE)))
        return
    vocabulary = _scout_words()
    lowered = [w.lower() for w in words]

    after_motif = set(words[i + 1] for i, w in enumerate(lowered[:-1])
                      if w in ("motif", "sequence", "pattern"))

    def is_name(word):
        low = word.lower()
        if low in vocabulary or word in after_motif:
            return False
        # part of a two-word target such as "salt bridges" or "free cysteines"
        return not any(low in phrase.split() for phrase in vocabulary if " " in phrase
                       and all(p in lowered for p in phrase.split()))

    def action(obj, rest, extra):
        what = " ".join(rest).lower().strip()
        found = None
        _ensure_ss(obj)
        if what in SCOUT:
            title, builder, blurb, style = SCOUT[what]
            found = (title, builder(obj), blurb, style)
        elif what in COMPUTED:
            found = _computed(obj, COMPUTED[what])
        elif rest and rest[0].lower() in ("motif", "sequence", "pattern") and \
                len(rest) == 2 and _MOTIF.match(rest[1]):
            found = _computed(obj, "motif", rest[1])
        if found is None:
            import difflib
            close = difflib.get_close_matches(what, sorted(vocabulary), n=1, cutoff=0.7)
            if not what:
                _say("Say what to look for, e.g.  scout %s histidine" % obj)
            else:
                codes = set()
                cmd.iterate("%s and not polymer and not solvent" % obj, "codes.add(resn)",
                            space={"codes": codes})
                bound = [c for c in codes if what == c.lower()
                         or what in places.LIGAND_WORDS.get(c, [])]
                if bound:
                    _say("%s is a molecule bound in %s, not something to spot. To see it:"
                         "  focus %s   (and what holds it:  pocket %s %s)"
                         % (what, obj, bound[0], obj, bound[0]))
                    return
                parts = what.split()
                if len(parts) > 1 and all(p in vocabulary for p in parts):
                    _say("scout spots one kind of thing at a time. Run them one after "
                         "the other:")
                    for part in parts[:4]:
                        _say("   scout %s %s" % (obj, part))
                    return
                _say("I don't know how to spot '%s'.%s" % (
                    what, " Did you mean '%s'?" % close[0] if close else ""))
            _scout_menu()
            return
        _show_scout(obj, found)

    _run(label, words, action, is_name=is_name)


def _show_scout(obj, found):
    title, selection, blurb, style = found
    name = "%s_scout" % obj
    cmd.delete(name)
    cmd.select(name, selection, enable=0)
    atoms = cmd.count_atoms(name)
    seen = set()
    cmd.iterate(name, "seen.add((chain, resn, resi))", space={"seen": seen})

    def order(item):
        digits = re.sub(r"[^0-9-]", "", item[2]) or "0"
        return (item[0], int(digits) if digits not in ("-",) else 0, item[2])
    listed = sorted(seen, key=order)
    count = len(listed)
    _say("")
    if not count:
        cmd.delete(name)
        _say("SPOTTED      none: %s has no %s." % (obj, title))
        _say("WHAT IT IS   " + blurb)
        return
    total = cmd.count_atoms("%s and polymer and name CA" % obj) or 1
    protein_like = style in ("sticks", "cartoon") and "polymer" in selection
    share = "  (%.0f%% of its %d residues)" % (100.0 * count / total, total) \
        if protein_like and "residues" in title else ""
    _say("SPOTTED      %d %s in %s%s" % (count, title, obj, share)
         if "bonds)" not in title else "SPOTTED      %s in %s" % (title, obj))
    by_chain = {}
    for chain, resn, resi in listed:
        by_chain.setdefault(chain, []).append("%s%s" % (resn, resi))
    if style == "balls" and "metals" in selection:
        # say which metal it is when it sits inside a larger group (Fe in a heme)
        atoms_found, by_chain = [], {}
        cmd.iterate(name, "atoms_found.append((chain, elem, resn, resi))",
                    space={"atoms_found": atoms_found})
        for chain, element, resn, resi in atoms_found:
            text = "%s%s" % (resn, resi) if element.upper() == resn.upper() \
                else "%s (in %s%s)" % (element, resn, resi)
            by_chain.setdefault(chain, []).append(text)
    shown_chains = 0
    for chain in sorted(by_chain):
        items = by_chain[chain]
        if shown_chains >= 4:
            _say("             (+%d more chains)" % (len(by_chain) - shown_chains))
            break
        text = " ".join(items[:14]) + (" ... +%d more" % (len(items) - 14)
                                      if len(items) > 14 else "")
        _say("             %s%s" % ("chain %s: " % chain if chain else "", text))
        shown_chains += 1

    # draw it
    cmd.enable(obj)
    base._clean(obj)
    if cmd.count_atoms("%s and rep surface" % obj):
        cmd.hide("surface", obj)
    cmd.show("cartoon", "%s and polymer" % obj)
    cmd.color("grey80", "%s and polymer and elem C" % obj)
    teach.note("select %s, <the residues listed above>" % name, "a named selection")
    teach.note("color grey80, %s and polymer and elem C" % obj, "everything else grey")
    if style == "cartoon" or count > 150:
        cmd.set("cartoon_transparency", 0, obj)
        cmd.color("orange", "%s and elem C" % name)
        teach.note("color orange, %s and elem C" % name, "paint that part of the ribbon")
        how = "those parts of the ribbon painted orange; the rest is grey"
    elif style == "balls":
        cmd.set("cartoon_transparency", 0.5, obj)
        cmd.show("spheres", "%s and (name CA or not polymer)" % name)
        cmd.set("sphere_scale", 0.8, name)
        cmd.color("orange", "%s and (name CA or not polymer)" % name)
        teach.note("show spheres, %s" % name)
        how = "an orange ball at each one, inside a faint grey ribbon"
    else:
        cmd.set("cartoon_transparency", 0.5, obj)
        target = name if style == "ligand" else "%s and (sidechain or name CA) or " \
                                                "(%s and not polymer)" % (name, name)
        cmd.show("sticks", target)
        cmd.color("orange", "%s and elem C" % name)
        cmd.util.cnc(name)
        teach.note("show sticks, %s and (sidechain or name CA)" % name)
        teach.note("color orange, %s and elem C" % name, "carbons orange")
        how = ("each one as orange sticks (blue = nitrogen, red = oxygen, yellow = "
               "sulfur) inside a faint grey ribbon")
        if count <= 25:
            cmd.label("%s and name CA" % name, '"%s%s" % (resn, resi)')
            how += ", with a name tag on each"
    if count <= 12 and atoms:
        cmd.orient(name)
        cmd.zoom(name, 10)
        teach.note("zoom %s, 10" % name)
    else:
        cmd.orient(obj)
    cmd.deselect()
    _say("ON SCREEN    " + how)
    _say("WHAT IT IS   " + blurb)
    first = listed[0]
    example = "%s%s" % ("%s/" % first[0] if first[0] else "", first[2])
    _say("GO CLOSER    focus %s  (one of them)  |  zoom %s  |  explain %s"
         % (example, name, example))
    _say("             back  undoes this view")


# ------------------------------------------------------------------ explain
def explain(*args, **kwargs):
    """
DESCRIPTION

    Explain one residue: which amino acid it is and what that means, where it
    sits (helix/strand/loop, buried/exposed), what it touches, and anything
    UniProt says about that position.

USAGE

    explain [protein] <chain/residue>

EXAMPLES

    explain H/57
    explain p53 A/248
    explain 1PPB 195
    """
    words = base._raw(args, kwargs).replace(",", " ").split()
    label = "explain %s" % " ".join(words)

    def action(obj, rest, extra):
        specs = [w for w in rest if _is_spec(w)]
        if len(specs) != 1:
            _say("Say which residue, e.g.  explain %s %s"
                 % (obj, places._example_residue(obj)))
            return
        selection = _spec_selection(obj, specs[0])
        if selection is None or not cmd.count_atoms("%s and polymer" % selection):
            _say("%s has no residue '%s'. Chains: %s. Example:  explain %s %s"
                 % (obj, specs[0], ", ".join(cmd.get_chains(obj)), obj,
                    places._example_residue(obj)))
            return
        residue = "(byres (first (%s and polymer)))" % selection
        _ensure_ss(obj)
        info = []
        cmd.iterate("first (%s and name CA)" % residue,
                    "info.append((chain, resn, resi, ss, b))", space={"info": info})
        if not info:
            _say("That residue has no backbone atoms in this model.")
            return
        chain, resn, resi, ss, b_factor = info[0]
        name, group, use = AA.get(resn.upper(), (resn, "non-standard", "a modified or "
                                                 "unusual residue"))
        # zoom there first, then describe
        place = places._residue_place(obj, "%s%s" % ("%s/" % chain if chain else "", resi))
        if place is not None:
            places._apply_quiet(obj, place)
        exposure = _residue_areas(obj).get((chain, resi))
        neighbours = set()
        cmd.iterate("byres ((%s and polymer) within 4 of (%s)) and not (%s)"
                    % (obj, residue, residue), "neighbours.add((chain, resn, resi))",
                    space={"neighbours": neighbours})
        nearby = set()
        cmd.iterate("(%s and not polymer and not solvent) within 4 of (%s)"
                    % (obj, residue), "nearby.add(resn)", space={"nearby": nearby})
        rows = _per_residue(obj)
        b_values = sorted(r[4] for r in rows)
        rank_b = sum(1 for v in b_values if v <= b_factor) / float(len(b_values) or 1)

        _say("")
        _say("RESIDUE      %s%s%s of %s" % (resn, resi, " in chain %s" % chain if chain
                                             else "", obj))
        _say("AMINO ACID   %s - %s" % (name, group))
        _say("             " + use)
        _say("")
        _say("WHERE        in %s" % SS_WORDS.get(ss, "a loop"))
        if exposure is not None:
            word = ("buried inside the protein" if exposure < 0.1 else
                    "partly exposed" if exposure < 0.5 else "on the surface")
            _say("             %s (%.0f%% of it can touch water)" % (word, 100 * exposure))
        move = ("one of the most rigid parts" if rank_b < 0.2 else
                "one of the most mobile parts" if rank_b > 0.85 else "average mobility")
        _say("             %s (B-factor %.0f)" % (move, b_factor))
        def order(item):
            digits = re.sub(r"[^0-9-]", "", item[2]) or "0"
            return (item[0], int(digits) if digits != "-" else 0)
        close = ["%s%s" % (n, i) for c, n, i in sorted(neighbours, key=order)]
        _say("TOUCHES      %d residues within 4 A: %s%s"
             % (len(close), " ".join(close[:10]), " ..." if len(close) > 10 else ""))
        if nearby:
            _say("             and the bound molecule%s %s"
                 % ("" if len(nearby) == 1 else "s", ", ".join(sorted(nearby))))
        # what UniProt knows about this position
        notes, variants, position = [], [], None
        annotation = extra.get("annotation")
        if annotation:
            for protein in annotation.proteins:
                for chain_map in protein.chains:
                    if chain_map.chain != chain:
                        continue
                    position = chain_map.uniprot_position(resi)
                    if position is None:
                        continue
                    for feature in protein.features:
                        if feature.group == "variant" and position in feature.positions():
                            variants.append(" ".join(feature.description.split()))
                        elif position in feature.positions() and feature.group != "domain":
                            notes.append("%s: %s" % (feature.kind, " ".join(
                                feature.description.split())[:110] or "yes"))
                        elif feature.group == "domain" and \
                                feature.start <= position <= feature.end:
                            notes.append(("inside the %s %s" % (
                                feature.description, feature.kind.lower())).replace(
                                    "  ", " ").replace("the dna binding", "the DNA-binding "
                                                                          "region"))
                    expected = protein.expected_residue(position)
                    if expected and feat.THREE_TO_ONE.get(resn.upper()) not in (None, expected):
                        notes.append("UniProt has %s here - this structure carries a "
                                     "mutation" % expected)
        _say("")
        if position is not None:
            _say("IN UNIPROT   position %d%s" % (position, "" if str(position) == resi
                                                 else "  (the file numbers it %s)" % resi))
        if notes:
            seen = []
            for note in notes:
                if note not in seen:
                    seen.append(note)
            for note in seen[:8]:
                _say("ANNOTATED    " + note)
            if len(seen) > 8:
                _say("             ... +%d more (often many disease variants)"
                     % (len(seen) - 8))
        elif annotation and not variants:
            _say("ANNOTATED    nothing special is recorded for this position")
        if variants:
            changes = []
            for text in variants:
                change = text.split(" ", 1)[0]
                if change not in changes:
                    changes.append(change)
            _say("DISEASE      %d known disease variant%s at this position: %s"
                 % (len(changes), "" if len(changes) == 1 else "s", ", ".join(changes[:10])))
            _say("             e.g. %s" % variants[-1][:150])
            # the suggestion uses the FILE's numbering (UniProt's can differ by one
            # or more: haemoglobin's Glu6 is UniProt position 7)
            first = changes[-1].replace(">", "%s" % resi) if ">" in changes[-1] else ""
            if re.match(r"^[A-Z]-?\d+[A-Z]?[A-Z]$", first):
                _say("             model one with:  whatif %s %s%s"
                     % (obj, "%s/" % chain if chain else "", first))
        _say("")
        _say("ON SCREEN    the residue as orange sticks with its name tag; the residues "
             "around it are the ribbon")
        _say("TRY          measure %s%s <another residue>   |   scout %s %s"
             % ("%s/" % chain if chain else "", resi, obj, name))

    _run(label, words, action, is_name=lambda w: not _is_spec(w), need_annotation=True)


# ------------------------------------------------------------------ measure
_measures = [0]


_METAL_WORDS = {"iron": "Fe", "fe": "Fe", "zinc": "Zn", "zn": "Zn", "copper": "Cu",
                "cu": "Cu", "magnesium": "Mg", "mg": "Mg", "calcium": "Ca",
                "manganese": "Mn", "mn": "Mn", "cobalt": "Co", "nickel": "Ni",
                "sodium": "Na", "potassium": "K", "metal": "Fe+Zn+Cu+Mg+Ca+Mn+Co+Ni"}


def measure(*args, **kwargs):
    """
DESCRIPTION

    Measure the distance between two residues, atoms or ligands, draw it,
    and say what that distance means.

USAGE

    measure [protein] <A> <B>

    A and B can be a residue (H/57), an atom (H/57/NE2) or a ligand code (HEM).

EXAMPLES

    measure H/57 H/195
    measure p53 A/176 ZN
    measure A/87/NE2 HEM
    """
    words = base._raw(args, kwargs).replace(",", " ").split()
    label = "measure %s" % " ".join(words)

    def action(obj, rest, extra):
        if len(rest) != 2:
            _say("measure needs two things, e.g.  measure %s %s <another residue>"
                 % (obj, places._example_residue(obj)))
            _say("Each can be a residue (H/57), an atom (H/57/NE2) or a ligand code (HEM).")
            return
        selections = [_spec_selection(obj, word) for word in rest]
        metal_words = []
        for index, word in enumerate(rest):
            element = _METAL_WORDS.get(word.lower())
            if selections[index] is None and element and \
                    cmd.count_atoms("%s and elem %s" % (obj, element)):
                # "the iron": a metal atom, also one inside a molecule (the Fe of heme)
                selections[index] = "(%s and elem %s)" % (obj, element)
                metal_words.append(index)
        for word, selection in zip(rest, selections):
            if selection is None:
                _say("I can't find '%s' in %s. Use chain/number (e.g. %s) or a ligand "
                     "code." % (word, obj, places._example_residue(obj)))
                return
        first = cmd.get_model(selections[0] + " and not elem H").atom
        second = cmd.get_model(selections[1] + " and not elem H").atom
        best = None
        for a in first:
            for b in second:
                d = sum((a.coord[k] - b.coord[k]) ** 2 for k in range(3)) ** 0.5
                if best is None or d < best[0]:
                    best = (d, a, b)
        distance, a, b = best
        for index in metal_words:          # of all the irons, the one that was measured
            atom = a if index == 0 else b
            selections[index] = "(%s and chain %s and resi %s and name %s)" % (
                obj, atom.chain or "''", atom.resi, atom.name)
        existing = cmd.get_names("all")
        _measures[0] = 1
        while "measure_%d" % _measures[0] in existing:      # first free name
            _measures[0] += 1
        name = "measure_%d" % _measures[0]

        def atom_sel(atom):
            chain = " and chain %s" % atom.chain if atom.chain else ""
            number = ("\\" + atom.resi) if atom.resi.startswith("-") else atom.resi
            return "%s%s and resi %s and name %s" % (obj, chain, number, atom.name)
        cmd.distance(name, atom_sel(a), atom_sel(b), zoom=0)
        teach.note("distance %s, %s, %s" % (name, atom_sel(a), atom_sel(b)),
                   "draws a dashed line with the length on it")
        base._clean(obj)
        for selection in selections:
            _show_residue(selection)
        cmd.zoom("(%s) or (%s)" % tuple(selections), 6)
        cmd.deselect()
        polar = a.symbol in ("N", "O") and b.symbol in ("N", "O")
        metal = cmd.count_atoms("(%s or %s) and metals" % (atom_sel(a), atom_sel(b))) > 0
        sulfur = a.symbol == "S" and b.symbol == "S"
        if sulfur and distance < 2.5:
            meaning = "a disulfide bond (two sulfur atoms bonded, about 2.0 A)"
        elif metal and distance < 2.8:
            meaning = "close enough to be bonded to the metal ion (typically 1.9-2.5 A)"
        elif distance < 1.8:
            meaning = "a covalent bond (atoms this close are chemically bonded)"
        elif polar and distance <= 3.5:
            meaning = "the right distance for a hydrogen bond (N or O atoms, 2.6-3.5 A)"
        elif distance <= 4.5:
            meaning = "in contact (atoms touch at about 3.5-4.5 A) but not bonded"
        elif distance <= 8:
            meaning = "near each other, but not touching"
        else:
            meaning = "far apart - they do not interact directly"
        _say("")
        _say("DISTANCE     %.2f A  between %s and %s" % (distance, _name_of(selections[0]),
                                                         _name_of(selections[1])))
        _say("             closest atoms: %s %s  -  %s %s" % (a.resn + a.resi, a.name,
                                                             b.resn + b.resi, b.name))
        _say("WHAT IT MEANS  " + meaning)
        ca_one = cmd.get_model("(%s) and name CA" % ("byres " + selections[0])).atom
        ca_two = cmd.get_model("(%s) and name CA" % ("byres " + selections[1])).atom
        if ca_one and ca_two:
            ca = sum((ca_one[0].coord[k] - ca_two[0].coord[k]) ** 2 for k in range(3)) ** 0.5
            _say("             (backbone to backbone, CA-CA: %.1f A)" % ca)
        _say("ON SCREEN    both shown as sticks, joined by a yellow dashed line with the "
             "distance written on it")
        _say("             1 A = 0.1 nanometre; a carbon-carbon bond is 1.5 A")
        _say("REMOVE IT    delete %s" % name)

    _run(label, words, action, protein_only=False,
         is_name=lambda w: not _is_spec(w) and not re.match(r"^[A-Z0-9]{1,3}$", w))


# ------------------------------------------------------------------- mutate
_MUTATION = re.compile(r"^(?:([A-Za-z0-9]{1,4})/)?([A-Za-z])(-?\d+[A-Za-z]?)([A-Za-z])$")
ONE_TO_THREE = dict((one, three) for three, one in feat.THREE_TO_ONE.items()
                    if three not in ("MSE", "SEC", "PYL"))


def mutate(*args, **kwargs):
    """
DESCRIPTION

    Model a point mutation and show what changes: the new side chain next to
    the original one, and the contacts that are lost or gained.  The original
    structure is not altered - the mutant is a separate object.

USAGE

    mutate [protein] <mutation>

    Write the mutation the usual way: original letter, position, new letter.

EXAMPLES

    mutate p53 R248Q
    mutate A/R248Q
    mutate haemoglobin B/E6V      (sickle-cell haemoglobin)
    """
    words = base._raw(args, kwargs).replace(",", " ").split()
    label = "mutate %s" % " ".join(words)

    def action(obj, rest, extra):
        specs = [w for w in rest if _MUTATION.match(w)]
        if len(specs) != 1:
            _say("Write one mutation as original letter + position + new letter, e.g.  "
                 "mutate %s R248Q" % obj)
            return
        chain, old, number, new = _MUTATION.match(specs[0]).groups()
        old, new = old.upper(), new.upper()
        if old not in ONE_TO_THREE or new not in ONE_TO_THREE:
            _say("'%s' and '%s' must be one-letter amino-acid codes (e.g. R = arginine, "
                 "Q = glutamine)." % (old, new))
            return
        # where is that residue?  Try the file's numbering, then UniProt numbering.
        rows = _per_residue(obj)
        matches = [(c, i) for c, i, n, s, b in rows
                   if i == number and (not chain or c == chain)
                   and feat.THREE_TO_ONE.get(n.upper()) == old]
        note = ""
        annotation = extra.get("annotation")
        if not matches and annotation and re.match(r"^\d+$", number):
            for protein in annotation.proteins:
                for chain_map in protein.chains:
                    if chain and chain_map.chain != chain:
                        continue
                    mapped = chain_map.resi(int(number))
                    if mapped and protein.expected_residue(int(number)) == old:
                        if any(c == chain_map.chain and i == mapped and
                               feat.THREE_TO_ONE.get(n.upper()) == old
                               for c, i, n, s, b in rows):
                            matches.append((chain_map.chain, mapped))
                            note = ("(%s%s is UniProt numbering; in this file that "
                                    "residue is number %s)" % (old, number, mapped))
        if not matches:
            actual = [(c, n) for c, i, n, s, b in rows if i == number
                      and (not chain or c == chain)]
            if actual:
                _say("Residue %s is %s in %s, not %s (%s). Check the position or the "
                     "letter." % (number, actual[0][1], obj, ONE_TO_THREE[old], old))
            else:
                _say("%s has no residue number %s%s. It may not be part of this "
                     "structure." % (obj, number, " in chain %s" % chain if chain else ""))
            return
        target_chain, target_resi = matches[0]
        label_resi = ("\\" + target_resi) if target_resi.startswith("-") else target_resi
        where = "%s and resi %s" % ("chain %s" % target_chain if target_chain
                                    else "polymer", label_resi)
        mutant = re.sub(r"[^A-Za-z0-9_]", "_", "%s_%s%s%s" % (obj, old, target_resi, new))
        _copy_for_change(obj, mutant)
        teach.note("create %s, %s" % (mutant, obj), "work on a copy, keep the original")
        try:
            changed = _swap(mutant, target_chain, target_resi, ONE_TO_THREE[new])
            teach.note("wizard mutagenesis", "PyMOL's built-in mutation tool "
                                             "(menu: Wizard > Mutagenesis)")
        except Exception as err:
            cmd.set_wizard()
            changed = False
            _say("PyMOL's mutagenesis tool failed here: %s" % str(err).strip()[:120])
        if not changed:
            cmd.delete(mutant)
            _say("PyMOL could not build %s at position %s, so nothing is shown. (This "
                 "can happen at chain ends or where atoms are missing.)"
                 % (AA[ONE_TO_THREE[new]][0], target_resi))
            return

        def contacts(source):
            found = set()
            cmd.iterate("byres ((%s and polymer) within 4 of (%s and %s and sidechain)) "
                        "and not (%s and %s)" % (source, source, where, source, where),
                        "found.add(resn + resi)", space={"found": found})
            ligands = set()
            cmd.iterate("(%s and not polymer and not solvent) within 4 of "
                        "(%s and %s and sidechain)" % (source, source, where),
                        "ligands.add(resn)", space={"ligands": ligands})
            return found, ligands
        before, before_lig = contacts(obj)
        after, after_lig = contacts(mutant)
        # draw: original green-grey, mutant magenta, same place
        for name in base._molecules(enabled_only=True):
            if name not in (obj, mutant):
                cmd.disable(name)
        cmd.hide("everything", mutant)
        base._clean(obj)
        cmd.set("cartoon_transparency", 0.6, obj)
        cmd.show("cartoon", "%s and polymer" % obj)
        cmd.show("sticks", "%s and %s" % (obj, where))
        cmd.color("grey50", "%s and %s and elem C" % (obj, where))
        cmd.show("sticks", "%s and %s" % (mutant, where))
        cmd.color("magenta", "%s and %s and elem C" % (mutant, where))
        cmd.util.cnc("(%s or %s) and %s" % (obj, mutant, where))
        nearby = "byres ((%s and polymer) within 4 of (%s and %s))" % (obj, obj, where)
        cmd.show("lines", nearby)
        cmd.label(obj, "")
        cmd.label("%s and %s and name CA" % (obj, where), '"%s%s" % (resn, resi)')
        cmd.orient("%s and %s" % (obj, where))
        cmd.zoom("%s and %s" % (obj, where), 9)
        cmd.deselect()

        old_name, new_name = AA[ONE_TO_THREE[old]], AA[ONE_TO_THREE[new]]
        _say("")
        _say("MUTATION     %s%s%s in %s%s:  %s  ->  %s"
             % (old, target_resi, new, obj, " chain %s" % target_chain if target_chain
                else "", old_name[0], new_name[0]))
        if note:
            _say("             " + note)
        _say("CHANGE       from %s to %s" % (old_name[1], new_name[1]))
        _say("             was: " + old_name[2])
        _say("             now: " + new_name[2])
        lost = sorted(before - after)
        gained = sorted(after - before)
        _say("CONTACTS     the original side chain touched %d residue%s; the new one "
             "touches %d" % (len(before), "" if len(before) == 1 else "s", len(after)))
        if lost:
            _say("             lost:   %s" % " ".join(lost[:10]))
        if gained:
            _say("             gained: %s" % " ".join(gained[:10]))
        if before_lig - after_lig:
            _say("             no longer touches the bound molecule%s %s"
                 % ("" if len(before_lig - after_lig) == 1 else "s",
                    ", ".join(sorted(before_lig - after_lig))))
        # is this a known disease variant?
        if annotation:
            for protein in annotation.proteins:
                for chain_map in protein.chains:
                    if chain_map.chain != target_chain:
                        continue
                    position = chain_map.uniprot_position(target_resi)
                    for feature in protein.features:
                        if feature.group == "variant" and feature.start == position and \
                                feature.description.startswith("%s>%s" % (old, new)):
                            _say("KNOWN        UniProt lists this exact variant: %s"
                                 % " ".join(feature.description.split("  ", 1)[-1]
                                            .split())[:120])
        _say("")
        _say("ON SCREEN    grey sticks = the original %s; magenta sticks = the new %s in "
             "the same place; thin lines = the residues around it"
             % (old_name[0], new_name[0]))
        _say("CAUTION      this is a quick model: only this one side chain was swapped, in "
             "its most likely shape. The rest of the protein was not allowed to adjust, "
             "so treat it as a sketch.")
        _say("OBJECTS      original: %s   |   mutant: %s  (hidden except the new side "
             "chain;  enable/show it to see more)" % (obj, mutant))

    _run(label, words, action, is_name=lambda w: not _MUTATION.match(w),
         need_annotation=True)


# ------------------------------------------------------------------- whatif
CHARGE = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1}
SIDE_ATOMS = {"GLY": 0, "ALA": 1, "SER": 2, "CYS": 2, "VAL": 3, "THR": 3, "PRO": 3,
              "ILE": 4, "LEU": 4, "ASP": 4, "ASN": 4, "MET": 4, "GLU": 5, "GLN": 5,
              "LYS": 5, "HIS": 6, "PHE": 7, "ARG": 7, "TYR": 8, "TRP": 10}
HYDROPATHY = {"ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5, "GLN": -3.5,
              "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5, "LEU": 3.8, "LYS": -3.9,
              "MET": 1.9, "PHE": 2.8, "PRO": -1.6, "SER": -0.8, "THR": -0.7, "TRP": -0.9,
              "TYR": -1.3, "VAL": 4.2}


def _copy_for_change(obj, mutant):
    """A copy of ``obj`` that PyMOL's mutagenesis tool can work on.  Where the
    authors modelled a residue in two alternative shapes, only the first is
    kept: the tool silently does nothing on residues that have two."""
    cmd.delete(mutant)
    cmd.create(mutant, obj, zoom=0)
    cmd.remove("%s and not alt ''+A" % mutant)
    cmd.alter(mutant, "alt=''")
    cmd.sort(mutant)


def _swap(mutant, chain, resi, new3):
    """Change one residue of ``mutant`` to ``new3`` and CHECK that it happened."""
    cmd.wizard("mutagenesis")
    cmd.refresh_wizard()
    wizard = cmd.get_wizard()
    wizard.set_mode(new3)
    wizard.do_select("/%s//%s/%s/" % (mutant, chain, resi))
    wizard.apply()
    cmd.set_wizard()
    cmd.delete("mutation")
    now = set()
    cmd.iterate("%s and %s" % (mutant, _where(chain, resi)), "now.add(resn)",
                space={"now": now})
    return now == set([new3])


def _locate(obj, chain, old, number, annotation, rows):
    """Where is residue <old><number>?  Tries the file's numbering, then UniProt
    numbering.  Returns ((chain, resi), note) or (None, reason)."""
    matches = [(c, i) for c, i, n, s, b in rows
               if i == number and (not chain or c == chain)
               and feat.THREE_TO_ONE.get(n.upper()) == old]
    if matches:
        return matches[0], ""
    if annotation and re.match(r"^\d+$", number):
        for protein in annotation.proteins:
            for chain_map in protein.chains:
                if chain and chain_map.chain != chain:
                    continue
                mapped = chain_map.resi(int(number))
                if mapped and protein.expected_residue(int(number)) == old and any(
                        c == chain_map.chain and i == mapped
                        and feat.THREE_TO_ONE.get(n.upper()) == old
                        for c, i, n, s, b in rows):
                    return (chain_map.chain, mapped), (
                        "%s%s is UniProt numbering; in this file that residue is "
                        "number %s" % (old, number, mapped))
    actual = [(c, n) for c, i, n, s, b in rows if i == number and (not chain or c == chain)]
    if actual:
        return None, ("Residue %s is %s in %s, not %s (%s). Check the position or the "
                      "letter." % (number, actual[0][1], obj, ONE_TO_THREE[old], old))
    return None, ("%s has no residue number %s%s. It may not be part of this structure."
                  % (obj, number, " in chain %s" % chain if chain else ""))


def _where(chain, resi):
    label = ("\\" + resi) if resi.startswith("-") else resi
    return "%s and resi %s" % ("chain %s" % chain if chain else "polymer", label)


def _side_chain_facts(source, where):
    """What one residue's side chain does in ``source``: neighbours, hydrogen-bond
    partners, bound molecules it touches, and atoms it bumps into."""
    side = "(%s and %s and (sidechain or (resn GLY and name CA)))" % (source, where)
    own = "(%s and %s)" % (source, where)
    neighbours, ligands = set(), set()
    cmd.iterate("byres ((%s and polymer) within 4 of %s) and not %s" % (source, side, own),
                "neighbours.add(resn + resi)", space={"neighbours": neighbours})
    cmd.iterate("(%s and not polymer and not solvent) within 4 of %s" % (source, side),
                "ligands.add(resn)", space={"ligands": ligands})
    ligands -= rcsb.ADDITIVES          # buffer components are not worth a mention
    mine = cmd.get_model("%s and (elem N or elem O) and not name N+O" % side).atom
    others = cmd.get_model("(%s and (elem N or elem O) and not solvent) and not %s "
                           "within 3.6 of %s" % (source, own, side)).atom
    bonds = set()
    for a in mine:
        for b in others:
            d = sum((a.coord[k] - b.coord[k]) ** 2 for k in range(3)) ** 0.5
            if d <= 3.5:
                bonds.add("%s%s" % (b.resn, b.resi))
    heavy = cmd.get_model("%s and not elem H and not name CA+CB" % side).atom
    around = cmd.get_model("(%s and not elem H and not solvent) and not %s within 3 of %s"
                           % (source, own, side)).atom
    clashes = set()
    for a in heavy:
        for b in around:
            d = sum((a.coord[k] - b.coord[k]) ** 2 for k in range(3)) ** 0.5
            if d < 2.5:
                clashes.add("%s%s" % (b.resn, b.resi))
    nucleic = set()
    cmd.iterate("byres ((%s and polymer.nucleic) within 4.5 of %s)" % (source, side),
                "nucleic.add(resn + resi)", space={"nucleic": nucleic})
    return {"neighbours": neighbours, "ligands": ligands, "bonds": bonds,
            "clashes": clashes, "nucleic": nucleic}


def _whatif_off(quiet=False):
    base._leave_panels()
    _whatif_restore()
    if not quiet:
        _say("")
        _say("Back to one panel. The changed version stays loaded (hidden); "
             " enable <name>  shows it,  delete <name>  removes it.")


_whatif_before = {}     # what whatif changed, so that  whatif off  can put it back


def _whatif_remember(obj):
    """Store how things look before whatif changes them: colours, what is drawn
    and which objects are shown (a scene without the camera), plus the
    see-through ribbon and stick thickness that a scene does not keep."""
    _whatif_before.clear()
    try:
        cmd.scene(base.WHATIF_SCENE, "store", view=0)
    except Exception:
        return
    _whatif_before.update(obj=obj, settings=dict(
        (name, (cmd.get("cartoon_transparency", name), cmd.get("stick_radius", name)))
        for name in base._molecules()))


def _whatif_restore():
    """Undo what whatif did to the original: its colours, ribbon and sticks, and
    the objects it hid (other structures, measurements, pocket dashes)."""
    if not _whatif_before:
        return
    names = cmd.get_names("all")
    if base.WHATIF_SCENE in (cmd.get_scene_list() or []):
        cmd.scene(base.WHATIF_SCENE, "recall", animate=0)
        cmd.scene(base.WHATIF_SCENE, "clear")
    for name, (transparency, radius) in _whatif_before.get("settings", {}).items():
        if name in names:
            cmd.set("cartoon_transparency", transparency, name)
            cmd.set("stick_radius", radius, name)
    obj = _whatif_before.get("obj")
    if obj in names:
        cmd.label(obj, "")
    mutant = _whatif_before.get("mutant")
    if mutant in names:             # the changed copy stays loaded, but hidden
        cmd.disable(mutant)
    _whatif_before.clear()


_OPEN_CHANGE = re.compile(r"^(?:([A-Za-z0-9]{1,4})/)?(first|last|-?\d+)([A-Za-z])?$", re.I)


def _fill_in_changes(obj, words, annotation):
    """Complete changes that leave something out: '248' (becomes what?), '248Q'
    (no original letter), 'first', 'last'.  The original letter is read from the
    structure; a missing new residue becomes alanine, the standard test.
    Returns (words, notes).  '1V 1A' (two halves of one change) is left alone
    so the caller can suggest V1A."""
    loose = [w for w in words if _OPEN_CHANGE.match(w) and not _MUTATION.match(w)]
    if not loose:
        return words, []
    numbers = [_OPEN_CHANGE.match(w).group(2) for w in loose]
    if len(loose) == 2 and numbers[0] == numbers[1] and len(words) == 2:
        return words, []
    rows = [r for r in _per_residue(obj) if r[2].upper() in feat.THREE_TO_ONE]
    out, notes = [], []
    for word in words:
        match = _OPEN_CHANGE.match(word)
        if _MUTATION.match(word) or not match or not rows:
            out.append(word)
            continue
        chain, position, new = match.groups()
        pool = [r for r in rows if not chain or r[0] == chain] or rows
        position = position.lower()
        row = None
        if position in ("first", "last"):
            row = pool[0] if position == "first" else pool[-1]
            notes.append("NOTE         the %s residue IN THIS STRUCTURE is %s%s (chain %s)."
                         % (position, row[2], row[1], row[0]))
            if position == "first" and row[1] != "1":
                notes.append("             Residues before it were not captured in the "
                             "experiment, so they cannot be shown.")
        else:
            found = [r for r in pool if r[1] == position]
            if not found and annotation:
                for protein in annotation.proteins:
                    for chain_map in protein.chains:
                        mapped = chain_map.resi(int(position))
                        found = found or [r for r in pool if r[0] == chain_map.chain
                                          and r[1] == mapped]
            row = found[0] if found else None
        if row is None:
            out.append(word)
            continue
        old = feat.THREE_TO_ONE[row[2].upper()]
        if not new:
            new = "G" if old == "A" else "A"
            notes.append("NOTE         you did not say what %s%s becomes, so I used %s: "
                         "the standard test, it cuts the side chain away. For another "
                         "one type e.g.  whatif %s %s%s%s"
                         % (row[2], row[1], "glycine" if new == "G" else "alanine", obj,
                            old, row[1], "D" if old != "D" else "K"))
        out.append("%s/%s%s%s" % (row[0], old, row[1], new.upper()) if row[0]
                   else "%s%s%s" % (old, row[1], new.upper()))
    return out, notes


def whatif(*args, **kwargs):
    """
DESCRIPTION

    "What if this residue were different?"  Changes one residue (or several)
    in a copy of the protein, shows the original and the changed version SIDE
    BY SIDE in two panels with the same camera, and explains the differences:
    charge, size, hydrogen bonds, contacts, clashes, and whether the change is
    a known disease variant.

USAGE

    whatif [protein] <change> [<change> ...]
    whatif off                       back to one panel

    Write a change as original letter + position + new letter (R248Q).
    Add a chain if needed (B/E6V).
    Shorter forms are completed for you: 248W (the original is read from
    the structure), 248 / first / last (becomes alanine, the standard test).

EXAMPLES

    whatif p53 R248Q
    whatif haemoglobin B/E6V
    whatif p53 R248Q R273H
    whatif off
    """
    words = base._raw(args, kwargs).replace(",", " ").replace("+", " ").split()
    label = "whatif %s" % " ".join(words)
    if [w.lower() for w in words] in (["off"], ["end"], ["stop"]):
        worker.submit(lambda: None, lambda r, e: (base._header(label), _whatif_off(),
                                                  _say(report.LINE)))
        return

    def action(obj, rest, extra):
        rest, open_notes = _fill_in_changes(obj, rest, extra.get("annotation"))
        for note in open_notes:
            _say(note)
        specs = [w for w in rest if _MUTATION.match(w)]
        if not specs or len(specs) != len(rest) or len(specs) > 5:
            wrong = [w for w in rest if not _MUTATION.match(w)]
            first = _per_residue(obj)[:1]
            example = "%s%sA" % (feat.THREE_TO_ONE.get(first[0][2].upper(), "G"),
                                 first[0][1]) if first else "R248Q"
            if wrong:
                _say("I can't read '%s' as a change." % " ".join(wrong))
                # "1V 1A" -> the user probably means V1A
                halves = [re.match(r"^(\d+)([A-Za-z])$", w) for w in wrong]
                if len(wrong) == 2 and all(halves) and \
                        halves[0].group(1) == halves[1].group(1):
                    guess = "%s%s%s" % (halves[0].group(2).upper(), halves[0].group(1),
                                        halves[1].group(2).upper())
                    _say("Did you mean  whatif %s %s  ?  (%s at position %s becomes %s)"
                         % (obj, guess, guess[0], halves[0].group(1), guess[-1]))
                    example = guess
            _say("Write a change as ONE word: original letter + position + new letter.")
            _say("   whatif %s %s     = at position %s, change %s into %s"
                 % (obj, example, example[1:-1], example[0], example[-1]))
            _say("   Up to five changes, each its own word:  whatif p53 R248Q R273H")
            _say("   See what is at each position with:  sequence %s" % obj)
            return
        annotation = extra.get("annotation")
        rows = _per_residue(obj)
        changes = []
        for spec in specs:
            chain, old, number, new = _MUTATION.match(spec).groups()
            old, new = old.upper(), new.upper()
            if old not in ONE_TO_THREE or new not in ONE_TO_THREE:
                _say("In '%s', '%s' and '%s' must be one-letter amino-acid codes (R = "
                     "arginine, Q = glutamine ...)." % (spec, old, new))
                return
            if old == new:
                _say("'%s' changes nothing: both letters are %s." % (spec, old))
                return
            place, note = _locate(obj, chain, old, number, annotation, rows)
            if place is None:
                _say(note)
                return
            changes.append((place[0], place[1], old, new, note))

        # build the changed version on a copy
        _whatif_off(quiet=True)
        _whatif_remember(obj)
        cmd.enable(obj)
        tag = "_".join("%s%s%s" % (old, resi, new) for _, resi, old, new, _ in changes)
        mutant = re.sub(r"[^A-Za-z0-9_]", "_", "%s_%s" % (obj, tag))
        _copy_for_change(obj, mutant)
        _whatif_before["mutant"] = mutant
        teach.note("create %s, %s" % (mutant, obj), "change a copy, keep the original")
        failed = ""
        try:
            for chain, resi, old, new, _ in changes:
                if not _swap(mutant, chain, resi, ONE_TO_THREE[new]):
                    failed = "%s%s%s" % (old, resi, new)
                    break
            teach.note("wizard mutagenesis", "PyMOL's mutation tool (Wizard > Mutagenesis)")
        except Exception as err:
            cmd.set_wizard()
            failed = str(err).strip()[:120]
        if failed:
            cmd.delete(mutant)
            _say("PyMOL could not build the change %s, so nothing is shown. (This can "
                 "happen at chain ends or where atoms are missing.)" % failed)
            return

        # draw both the same way, then put them in two panels.  Everything else
        # is hidden - not only other structures but also dashed lines left by
        # pocket or measure, which would each get a panel of their own
        for name in cmd.get_names("objects", enabled_only=1) or []:
            if name not in (obj, mutant):
                cmd.disable(name)
        cmd.enable(mutant)
        every = " or ".join("(%s)" % _where(c, r) for c, r, _, _, _ in changes)
        for name, colour in ((obj, "green"), (mutant, "magenta")):
            base._clean(name)
            cmd.hide("everything", "%s and solvent" % name)
            cmd.show("cartoon", "%s and polymer" % name)
            cmd.set("cartoon_transparency", 0.6, name)
            cmd.color("grey80", "%s and polymer and elem C" % name)
            around = "byres ((%s and not solvent) within 4 of (%s and (%s)))" % (
                name, name, every)
            cmd.show("sticks", around)
            cmd.set("stick_radius", 0.12, name)
            cmd.set("stick_radius", 0.3, "%s and (%s)" % (name, every))
            cmd.color(colour, "%s and (%s) and elem C" % (name, every))
            cmd.util.cnc(around)
            cmd.label("%s and (%s) and name CA" % (name, every), '"%s%s" % (resn, resi)')
        for name, source in (("whatif_bonds_before", obj), ("whatif_bonds_after", mutant)):
            cmd.delete(name)
            side = "%s and (%s) and sidechain" % (source, every)
            if cmd.count_atoms(side):
                cmd.distance(name, side, "(%s and not solvent) and not (%s and (%s))"
                             % (source, source, every), mode=2, zoom=0)
                cmd.hide("labels", name)
        # a group shares one panel, so each structure keeps its own dashes
        for group, members in (("whatif_original", [obj, "whatif_bonds_before"]),
                               ("whatif_changed", [mutant, "whatif_bonds_after"])):
            present = [m for m in members if m in cmd.get_names("all")]
            cmd.group(group, " ".join(present))
            cmd.group(group, action="open")
        cmd.order("whatif_original whatif_changed", location="top")
        cmd.set("grid_mode", 1)
        teach.note("group whatif_original, %s" % obj, "a group shares one panel")
        teach.note("group whatif_changed, %s" % mutant)
        teach.note("set grid_mode, 1", "one panel per object or group, same camera")
        first = "%s and (%s)" % (obj, _where(changes[0][0], changes[0][1]))
        focus_sel = "%s and (%s)" % (obj, every)
        cmd.orient(first)
        cmd.zoom(focus_sel, 9)
        cmd.deselect()

        # ---- explain
        names = ", ".join("%s%s%s" % (old, resi, new) for _, resi, old, new, _ in changes)
        _say("")
        _say("WHAT IF      %s in %s" % (names, obj))
        _say("ON SCREEN    two panels with the same camera - turn one and both turn")
        _say("             LEFT   the original; the residue in question is green")
        _say("             RIGHT  the changed version; the new residue is magenta")
        _say("             thin sticks = everything within 4 A; yellow dashes = hydrogen "
             "bonds of that side chain; grey ribbon = the rest")
        exposure = _residue_areas(obj)
        verdicts = []
        for chain, resi, old, new, note in changes:
            where = _where(chain, resi)
            old3, new3 = ONE_TO_THREE[old], ONE_TO_THREE[new]
            before = _side_chain_facts(obj, where)
            after = _side_chain_facts(mutant, where)
            _say("")
            _say("%s%s%s%s:  %s  ->  %s" % (old, resi, new, " (chain %s)" % chain if chain
                                           else "", AA[old3][0], AA[new3][0]))
            if note:
                _say("   (%s)" % note)
            reasons = []
            # kind of side chain
            _say("   %-11s %s  ->  %s" % ("Kind", AA[old3][1], AA[new3][1]))
            # charge
            sign = {1: "positive (+)", -1: "negative (-)", 0: "none"}
            c_old, c_new = CHARGE.get(old3, 0), CHARGE.get(new3, 0)
            if c_old != c_new:
                if c_old and c_new:
                    text, why = "the charge is reversed", "reverses a charge"
                elif c_old:
                    text, why = "a charge is lost", "removes a %s charge" % sign[c_old][:8]
                else:
                    text, why = "a charge appears", "adds a %s charge" % sign[c_new][:8]
                reasons.append(why)
            else:
                text = "unchanged"
            _say("   %-11s %s  ->  %s   (%s)" % ("Charge", sign[c_old], sign[c_new], text))
            # size
            s_old, s_new = SIDE_ATOMS[old3], SIDE_ATOMS[new3]
            if s_new - s_old >= 3:
                text = "much bigger - it may not fit"
                reasons.append("is much bulkier")
            elif s_old - s_new >= 3:
                text = "much smaller - it may leave an empty space"
                reasons.append("is much smaller")
            else:
                text = "about the same size" if abs(s_new - s_old) <= 1 else (
                    "bigger" if s_new > s_old else "smaller")
            _say("   %-11s %d  ->  %d side-chain atoms   (%s)" % ("Size", s_old, s_new, text))
            # water-loving or not
            h_old, h_new = HYDROPATHY[old3], HYDROPATHY[new3]
            def word(value):
                return "oily" if value > 1.5 else "water-loving" if value < -1.5 \
                    else "in between"
            exposed = exposure.get((chain, resi))
            place = "" if exposed is None else (
                "buried" if exposed < 0.1 else "partly exposed" if exposed < 0.5
                else "on the surface")
            text = "%s  ->  %s" % (word(h_old), word(h_new))
            if word(h_old) != word(h_new):
                if word(h_new) == "oily" and place == "on the surface":
                    text += "   (an oily patch on the surface: can make proteins stick " \
                            "together)"
                    reasons.append("puts an oily side chain on the surface")
                elif word(h_new) == "water-loving" and place == "buried":
                    text += "   (a water-loving group in the oily core: destabilising)"
                    reasons.append("puts a water-loving side chain in the core")
            _say("   %-11s %s" % ("Water", text))
            if place:
                _say("   %-11s %s (%.0f%% of the original side chain can touch water)"
                     % ("Location", place, 100 * exposed))
            # hydrogen bonds
            lost = sorted(before["bonds"] - after["bonds"])
            gained = sorted(after["bonds"] - before["bonds"])
            _say("   %-11s %d  ->  %d%s%s" % (
                "H-bonds", len(before["bonds"]), len(after["bonds"]),
                "   lost: %s" % " ".join(lost[:6]) if lost else "",
                "   new: %s" % " ".join(gained[:6]) if gained else ""))
            if len(lost) > len(gained):
                reasons.append("loses %d hydrogen bond%s" % (len(lost) - len(gained),
                               "" if len(lost) - len(gained) == 1 else "s"))
            # neighbours
            lost_n = sorted(before["neighbours"] - after["neighbours"])
            gained_n = sorted(after["neighbours"] - before["neighbours"])
            _say("   %-11s touches %d  ->  %d residues%s%s" % (
                "Contacts", len(before["neighbours"]), len(after["neighbours"]),
                "   lost: %s" % " ".join(lost_n[:6]) if lost_n else "",
                "   new: %s" % " ".join(gained_n[:6]) if gained_n else ""))
            gone = sorted(before["ligands"] - after["ligands"])
            if gone:
                _say("   %-11s no longer touches the bound molecule%s %s"
                     % ("", "" if len(gone) == 1 else "s", ", ".join(gone)))
                reasons.append("stops touching %s" % ", ".join(gone))
            if before["ligands"] & after["ligands"]:
                _say("   %-11s still touches %s" % ("", ", ".join(sorted(
                    before["ligands"] & after["ligands"]))))
            if before["nucleic"]:
                still = "still reaches it" if after["nucleic"] else "no longer reaches it"
                _say("   %-11s the original side chain touches DNA/RNA (%s); the new one "
                     "%s" % ("DNA or RNA", " ".join(sorted(before["nucleic"])[:4]), still))
                if c_old == 1 and c_new <= 0:
                    reasons.append("removes a positive charge that holds the negatively "
                                   "charged DNA/RNA")
                elif not after["nucleic"]:
                    reasons.append("loses its contact with DNA/RNA")
            # clashes
            if after["clashes"]:
                _say("   %-11s the new side chain bumps into %s (atoms closer than "
                     "2.5 A) - the protein would have to shift to make room"
                     % ("Clashes", " ".join(sorted(after["clashes"])[:5])))
                reasons.append("does not fit without the neighbours moving")
            else:
                _say("   %-11s none - the new side chain fits where the old one was"
                     % "Clashes")
            # known variant?
            if annotation:
                for protein in annotation.proteins:
                    for chain_map in protein.chains:
                        if chain_map.chain != chain:
                            continue
                        position = chain_map.uniprot_position(resi)
                        for feature in protein.features:
                            if feature.group == "variant" and feature.start == position \
                                    and feature.description.startswith("%s>%s" % (old, new)):
                                _say("   %-11s UniProt lists this exact change: %s" % (
                                    "Known", " ".join(feature.description.split("  ", 1)[-1]
                                                      .split())[:110]))
                            elif feature.group in ("active", "binding") and \
                                    position in feature.positions():
                                _say("   %-11s this residue is annotated as %s%s" % (
                                    "Role", feature.kind.lower(),
                                    " (%s)" % feature.description[:40]
                                    if feature.description else ""))
                                reasons.append("changes a %s residue" % feature.kind.lower())
            verdicts.append((old, resi, new, reasons))
        _say("")
        for old, resi, new, reasons in verdicts:
            if reasons:
                _say("IN SHORT     %s%s%s %s." % (old, resi, new, "; ".join(
                    dict.fromkeys(reasons))))
            else:
                _say("IN SHORT     %s%s%s is a mild change by these measures: similar "
                     "size and chemistry, no contacts lost." % (old, resi, new))
        _say("CAUTION      a sketch, not a prediction: only the side chain%s swapped, in "
             "%s most likely shape; the rest of the protein was not allowed to adjust."
             % (" was" if len(changes) == 1 else "s were",
                "its" if len(changes) == 1 else "their"))
        _say("WHEN DONE    whatif off  (one panel again)  |  back  (the view from before)")
        _say("             original = %s, changed version = %s" % (obj, mutant))

    def is_name(word):
        # short words with a digit ("1V", "V1", "248") are mistyped changes,
        # not part of a protein name - they must never be sent to the search
        return not _MUTATION.match(word) and not _OPEN_CHANGE.match(word)

    _run(label, words, action, is_name=is_name, need_annotation=True)



# ----------------------------------------------------------------- sequence
MARKS = [("active", "A"), ("binding", "B"), ("disulfide", "C"), ("ptm", "P"),
         ("variant", "v")]


def sequence(*args, **kwargs):
    """
DESCRIPTION

    Print the amino-acid sequence of a chain with its secondary structure
    and UniProt features underneath, so you can link the 1D sequence to the
    3D structure.

USAGE

    sequence [protein] [chain=A]

EXAMPLES

    sequence p53
    sequence haemoglobin chain=B
    """
    words, options = base._split_options(base._raw(args, kwargs))
    label = "sequence %s" % base._raw(args, kwargs)

    def action(obj, rest, extra):
        rows = _per_residue(obj)
        chains = []
        for row in rows:
            if row[0] not in chains:
                chains.append(row[0])
        if not rows:
            _say("%s has no protein chain to print." % obj)
            return
        chain = options.get("chain", chains[0])
        if chain not in chains:
            _say("%s has no protein chain '%s'. Chains: %s" % (obj, chain,
                                                               ", ".join(chains)))
            return
        mine = [row for row in rows if row[0] == chain]
        letters = "".join(feat.THREE_TO_ONE.get(r[2].upper(), "X") for r in mine)
        structure = "".join({"H": "H", "S": "E"}.get(r[3], "-") for r in mine)
        marks = [" "] * len(mine)
        index = dict((r[1], i) for i, r in enumerate(mine))
        used = []
        annotation = extra.get("annotation")
        name = ""
        if annotation:
            for protein in annotation.proteins:
                if not any(c.chain == chain for c in protein.chains):
                    continue
                name = protein.protein.name
                for group, mark in MARKS:
                    for feature in protein.features:
                        if feature.group != group:
                            continue
                        for label_resi in feature.placed.get(chain, []):
                            if label_resi in index and marks[index[label_resi]] in (" ", "v"):
                                marks[index[label_resi]] = mark
                                if mark not in used:
                                    used.append(mark)
        helix = structure.count("H") * 100.0 / len(mine)
        strand = structure.count("E") * 100.0 / len(mine)
        _say("")
        _say("SEQUENCE     %s chain %s%s: %d residues modelled (%.0f%% helix, %.0f%% strand)"
             % (obj, chain, " - %s" % name if name else "", len(mine), helix, strand))
        _say("")
        width = 60
        for start in range(0, len(mine), width):
            block = slice(start, start + width)
            first, last = mine[start][1], mine[min(start + width, len(mine)) - 1][1]
            _say("%6s  %s  %s" % (first, letters[block], last))
            _say("        %s" % structure[block])
            feature_line = "".join(marks[block]).rstrip()
            if feature_line:
                _say("        %s" % feature_line)
            _say("")
        _say("HOW TO READ  line 1: one letter per amino acid, with the residue numbers of "
             "the row's first and last residue")
        _say("             line 2: H = alpha-helix, E = beta-strand, - = loop")
        if used:
            legend = {"A": "A = active site", "B": "B = binding site",
                      "C": "C = disulfide cysteine", "P": "P = modified residue (PTM)",
                      "v": "v = disease variant position"}
            _say("             line 3: " + ", ".join(legend[m] for m in used))
        if len(chains) > 1:
            _say("OTHER CHAINS %s   e.g.  sequence %s chain=%s"
                 % (", ".join(c for c in chains if c != chain), obj,
                    [c for c in chains if c != chain][0]))
        _say("SEE IT IN 3D focus %s%s  (any residue number)  |  scout %s helices"
             % ("%s/" % chain if chain else "", mine[len(mine) // 2][1], obj))
        _say("             Only residues present in the model are printed; disordered "
             "ends and loops are missing.")

    _run(label, words, action, need_annotation=True, snapshot=False)


# --------------------------------------------------------------------- back
def back(*args, **kwargs):
    """
DESCRIPTION

    Undo the last Prospector command that changed the picture (look, viewas,
    focus, scout, features, pocket, compare, tour, mutate, measure ...).
    Type it again to go further back (up to 10 steps).

USAGE

    back             one step
    back 3           three steps
    """
    words = base._raw(args, kwargs).split()
    steps = max(1, min(10, int(words[0]))) if words and words[0].isdigit() else 1

    def action():
        base._header("back")
        _say("")
        if _cleared[0] and not base._molecules():      # the last thing was `clear`
            path, _cleared[0] = _cleared[0], None
            cmd.load(path)
            try:
                with open(path[:-4] + ".json", "r", encoding="utf-8") as handle:
                    data = json.load(handle)
                base._objects.update(dict((k, tuple(v)) for k, v in data["objects"].items()))
                base._aliases.update(data["aliases"])
                base._last[0] = data.get("last")
            except (OSError, ValueError, KeyError):
                pass
            names = base._molecules()
            _say("WENT BACK    to before  clear : %d structure%s on screen again%s"
                 % (len(names), "" if len(names) == 1 else "s",
                    " (%s)" % ", ".join(names[:6]) if names else ""))
            _say(report.LINE)
            return
        for step in range(steps):
            if step and not base._undo:
                break
            _say(base._restore())
        _say(report.LINE)

    worker.submit(lambda: None, lambda result, error: action())


# -------------------------------------------------------------------- clear
_cleared = [None]     # the rescue file written by the last `clear`, for `back`


def clear(*args, **kwargs):
    """
DESCRIPTION

    Empty the screen: removes every loaded structure, selection, measurement
    and highlight so you can start fresh.  A copy is kept first, so  back
    brings everything back.

USAGE

    clear            remove everything from the 3D window
    clear text       only wipe the text in the output area
    back             (right after clear) put everything back
    """
    words = base._raw(args, kwargs).lower().split()

    def action():
        if words and words[0] in ("text", "log", "output", "terminal"):
            cmd.do("cls", echo=0)
            return
        base._header("clear")
        _say("")
        if words:
            _say("Use:  clear  (empty the 3D window)  |  clear text  (wipe this text)")
            _say("To remove one structure only:  delete <its name>")
            _say(report.LINE)
            return
        names = cmd.get_names("all")
        if not names:
            _say("The 3D window is already empty.")
            _say(report.LINE)
            return
        molecules = base._molecules()
        rescue = os.path.join(_sessions_dir(), "before_clear.pse")
        try:
            cmd.save(rescue)
            with open(rescue[:-4] + ".json", "w", encoding="utf-8") as handle:
                json.dump({"objects": base._objects, "aliases": base._aliases,
                           "last": base._last[0]}, handle)
            _cleared[0] = rescue
        except Exception:
            _cleared[0] = None
        _whatif_off(quiet=True)
        cmd.delete("all")
        cmd.set("grid_mode", 0)
        del base._undo[:]
        base._objects.clear()
        base._aliases.clear()
        base._last[0] = None
        teach.note("delete all", "removes every object and selection")
        _say("CLEARED      %d structure%s removed%s; the 3D window is empty."
             % (len(molecules), "" if len(molecules) == 1 else "s",
                " (%s)" % ", ".join(molecules[:6]) if molecules else ""))
        if _cleared[0]:
            _say("CHANGED YOUR MIND?   back   puts everything back (until you load "
                 "something new; after that:  session open before_clear)")
        _say("WHAT NEXT    grab <protein>   |   surprise   |   clear text  (wipe this "
             "text too)")
        _say(report.LINE)

    worker.submit(lambda: None, lambda result, error: action())


# ------------------------------------------------------------------ session
def _sessions_dir():
    path = os.path.join(cache.home(), "sessions")
    os.makedirs(path, exist_ok=True)
    return path


def _session_name(text):
    return re.sub(r"[^A-Za-z0-9_\-]", "_", text.strip())[:60]


def session(*args, **kwargs):
    """
DESCRIPTION

    Keep your work between PyMOL runs: every loaded structure, its colours,
    selections and view, plus the names you gave things.

USAGE

    session                  list saved sessions
    session save <name>
    session open <name>
    session delete <name>

EXAMPLES

    session save p53_homework
    session open p53_homework
    """
    words = base._raw(args, kwargs).split()

    def action():
        base._header("session %s" % " ".join(words))
        _say("")
        verb = words[0].lower() if words else "list"
        folder = _sessions_dir()
        saved = sorted(f[:-4] for f in os.listdir(folder)
                       if f.endswith(".pse") and not f.endswith(".bak.pse"))
        if verb in ("list", "ls"):
            if not saved:
                _say("No saved sessions yet. Save your work with:  session save my_name")
            else:
                _say("SAVED SESSIONS   (open one with:  session open <name>)")
                for name in saved:
                    path = os.path.join(folder, name + ".pse")
                    _say("   %-28s %s   %.1f MB" % (name, time.strftime(
                        "%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(path))),
                        os.path.getsize(path) / 1e6))
            _say(report.LINE)
            return
        name = _session_name(" ".join(words[1:])) if len(words) > 1 else ""
        if verb == "save" and not name and base._molecules():
            name = _session_name(base._last[0] or base._molecules()[-1])   # a default
        if verb not in ("save", "open", "load", "delete", "rm") or not name:
            _say("Use:  session save <name>  |  session open <name>  |  session  (list)  "
                 "|  session delete <name>")
            _say(report.LINE)
            return
        path = os.path.join(folder, name + ".pse")
        side = os.path.join(folder, name + ".json")
        if verb == "save":
            if not base._molecules():
                _say("Nothing is loaded, so there is nothing to save.")
                _say(report.LINE)
                return
            replaced = os.path.exists(path)
            if replaced:                      # keep the previous version as a backup
                os.replace(path, os.path.join(folder, name + ".bak.pse"))
            cmd.save(path)
            with open(side, "w", encoding="utf-8") as handle:
                json.dump({"objects": base._objects, "aliases": base._aliases}, handle)
            teach.note("save %s" % path, "a .pse file holds the whole PyMOL session")
            _say("SAVED        '%s' - %d structure%s with colours, selections and view"
                 % (name, len(base._molecules()), "" if len(base._molecules()) == 1
                    else "s"))
            if replaced:
                _say("             (the earlier '%s' was kept as a backup: %s.bak.pse)"
                     % (name, name))
            _say("OPEN LATER   session open %s" % name)
            _say("FILE         %s" % path)
        elif verb in ("open", "load"):
            if not os.path.exists(path):
                import difflib
                close = difflib.get_close_matches(name, saved, n=1, cutoff=0.5)
                _say("There is no saved session called '%s'.%s" % (
                    name, " Did you mean '%s'?" % close[0] if close else ""))
                _say("Saved sessions: %s" % (", ".join(saved) or "none"))
                _say(report.LINE)
                return
            if base._molecules():             # opening replaces everything: keep a copy
                rescue = os.path.join(folder, "before_last_open.pse")
                cmd.save(rescue)
                with open(rescue[:-4] + ".json", "w", encoding="utf-8") as handle:
                    json.dump({"objects": base._objects, "aliases": base._aliases}, handle)
                _say("(What was on screen was saved first as 'before_last_open'.)")
            cmd.load(path)
            teach.note("load %s" % path, "loading a .pse replaces the current session")
            try:
                with open(side, "r", encoding="utf-8") as handle:
                    data = json.load(handle)
                base._objects.clear()
                base._objects.update(dict((k, tuple(v)) for k, v in data["objects"].items()))
                base._aliases.clear()
                base._aliases.update(data["aliases"])
            except (OSError, ValueError, KeyError):
                pass
            names = base._molecules()
            _say("OPENED       '%s' - %d structure%s: %s"
                 % (name, len(names), "" if len(names) == 1 else "s", ", ".join(names[:8])))
            known = sorted(set(a for a, o in base._aliases.items() if o in names))
            if known:
                _say("             names you can use again: %s" % ", ".join(known[:8]))
        else:
            if not os.path.exists(path):
                _say("There is no saved session called '%s'." % name)
            else:
                os.replace(path, os.path.join(folder, name + ".bak.pse"))
                if os.path.exists(side):
                    os.remove(side)
                _say("Removed '%s' from the list (a backup copy stays as %s.bak.pse in "
                     "%s)." % (name, name, folder))
        _say(report.LINE)

    worker.submit(lambda: None, lambda result, error: action())


# ------------------------------------------------------------------- report
def report_command(*args, **kwargs):
    """
DESCRIPTION

    Write a one-page report about a structure: its info card, key places and
    three pictures, in a single HTML file you can open in a browser, print or
    hand in.

USAGE

    report [protein]

EXAMPLES

    report p53
    report haemoglobin
    """
    words = base._raw(args, kwargs).split()
    label = "report %s" % " ".join(words)
    def action(obj, rest, extra):
        identity = extra.get("identity")
        text = extra.get("card")
        places.build(obj, extra.get("annotation"), extra.get("names"))
        folder = os.getcwd() if os.access(os.getcwd(), os.W_OK) else \
            os.path.join(cache.home(), "reports")
        os.makedirs(folder, exist_ok=True)
        stem, number = os.path.join(folder, "%s_report" % obj), 2
        path = stem + ".html"
        while os.path.exists(path):
            path = "%s_%d.html" % (stem, number)
            number += 1
        _say("")
        _say("Drawing three pictures - PyMOL is busy for a few seconds...")
        pictures = []
        view_w, view_h = cmd.get_viewport()
        width = 900
        height = int(round(width * float(view_h) / max(1, view_w)))
        saved_bg = cmd.get("bg_rgb")
        cmd.bg_color("white")

        def shot(caption):
            temp = os.path.join(cache.home(), "_report_shot.png")
            if os.path.exists(temp):
                os.remove(temp)
            cmd.png(temp, width=width, height=height, dpi=150, ray=1, quiet=1)
            if os.path.exists(temp):
                with open(temp, "rb") as handle:
                    pictures.append((caption, base64.b64encode(handle.read()).decode("ascii")))
                os.remove(temp)

        for name in base._molecules(enabled_only=True):
            if name != obj:
                cmd.disable(name)
        cmd.set("cartoon_transparency", 0, obj)
        cmd.label(obj, "")
        cmd.hide("everything", obj)
        cmd.show("cartoon", "%s and polymer" % obj)
        cmd.util.cbc(obj) if len(cmd.get_chains("%s and polymer" % obj)) > 1 \
            else cmd.util.chainbow(obj)
        cmd.show("sticks", "%s and organic" % obj)
        cmd.orient(obj)
        shot("Overall fold: ribbons follow the backbone (spirals = alpha-helices, arrows = "
             "beta-strands); sticks are bound molecules.")
        listed = places._places.get(obj, [])
        special = [p for p in listed if p.kind in ("active", "ligand", "binding")]
        if special:
            places._apply_quiet(obj, special[0])
            shot("%s: %s. %s." % (special[0].title, special[0].detail,
                                   places.MEANING.get(special[0].kind, "").capitalize()))
        cmd.set("cartoon_transparency", 0, obj)
        cmd.label(obj, "")
        cmd.hide("everything", obj)
        cmd.show("surface", "%s and polymer" % obj)
        for colour_name, rgb, residues, _ in base_hydro():
            cmd.set_color(colour_name, rgb)
            cmd.color(colour_name, "%s and polymer and resn %s" % (obj, residues))
        cmd.orient(obj)
        shot("Surface coloured by hydrophobicity: orange = water-avoiding patches, blue = "
             "charged, white = in between.")
        cmd.bg_color("black")
        try:
            cmd.set("bg_rgb", saved_bg)
        except Exception:
            pass
        base._restore(silent=True)       # put the screen back as it was

        rows = "".join("<tr><td>%d</td><td>%s</td><td>%s</td><td><code>focus %d</code></td>"
                       "</tr>" % (i, html.escape(p.title), html.escape(p.detail), i)
                       for i, p in enumerate(listed, 1))
        figures = "".join('<figure><img src="data:image/png;base64,%s" alt="%s">'
                          '<figcaption>%s</figcaption></figure>'
                          % (data, html.escape(caption), html.escape(caption))
                          for caption, data in pictures)
        title = "%s - structure report" % obj
        page = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>%s</title>
<style>
 body{font-family:Georgia,serif;max-width:900px;margin:2em auto;padding:0 1em;color:#222}
 h1{font-size:1.6em;margin-bottom:.2em} h2{font-size:1.15em;margin-top:1.6em}
 pre{background:#f6f6f3;padding:1em;overflow-x:auto;font-size:.85em;line-height:1.35}
 table{border-collapse:collapse;width:100%%;font-size:.95em}
 td,th{border-bottom:1px solid #ddd;padding:.35em .6em;text-align:left}
 figure{margin:1.2em 0} img{max-width:100%%;border:1px solid #ddd}
 figcaption{font-size:.9em;color:#555;margin-top:.3em}
 .note{color:#666;font-size:.9em}
</style></head><body>
<h1>%s</h1>
<p class="note">Made with Prospector for PyMOL on %s. Data: RCSB PDB, UniProt%s.</p>
<h2>What it is</h2>
<pre>%s</pre>
<h2>Key places</h2>
<table><tr><th>#</th><th>Place</th><th>Details</th><th>Zoom there in PyMOL</th></tr>%s</table>
<h2>Pictures</h2>
%s
<h2>Reproduce this in PyMOL</h2>
<pre>grab %s
about %s
tour %s</pre>
</body></html>""" % (html.escape(title), html.escape(title), time.strftime("%Y-%m-%d"),
                     ", AlphaFold DB" if identity and identity[0] == "af" else "",
                     html.escape(text or "(the info card could not be fetched)"),
                     rows, figures, obj, obj, obj)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(page)
        _say("SAVED        %s" % path)
        _say("             %d picture%s, the info card and %d key places, in one file "
             "(%.1f MB)" % (len(pictures), "" if len(pictures) == 1 else "s", len(listed),
                            os.path.getsize(path) / 1e6))
        _say("OPEN IT      double-click the file: it opens in your web browser, where you "
             "can read or print it")
        _say("ON SCREEN    unchanged - the view was put back as it was")

    def fetch_card(identity, extra):
        try:
            if identity and identity[0] == "pdb":
                extra["card"] = card.entry_card(identity[1])
            elif identity:
                from ..core import alphafold as af_api
                model = af_api.prediction(identity[1])
                extra["card"] = card.model_card(model) if model else None
        except http.NetworkError:
            extra["card"] = None

    _run(label, words, action, need_annotation=True, extra_job=fetch_card,
         protein_only=False)


def base_hydro():
    from . import tier3
    return tier3._HYDRO


def register():
    table = {"scout": scout, "explain": explain, "measure": measure, "mutate": mutate,
             "sequence": sequence, "back": back, "session": session, "clear": clear,
             "whatif": whatif,
             "report": report_command}
    for name, function in table.items():
        cmd.extend(name, function)
        cmd.extend("p_" + name, function)
    try:
        cmd.auto_arg[0]["session"] = [lambda: cmd.Shortcut(["save", "open", "delete"]),
                                      "action", " "]
    except Exception:
        pass
