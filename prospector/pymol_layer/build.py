"""`build`: make biological structures that are not simply downloaded.

    build dna ATGCGTACGT [form=A]      a DNA double helix of any sequence
    build rna AUGGCUACG                a single RNA strand (A-form helix)
    build helix AKLEELLKKL             a peptide as an alpha-helix
    build strand KVTVKV                a peptide as an extended beta-strand
    build sheet KVTVKV TVKVTV          two or more strands as a beta-sheet
    build peptide ACDEFGH              a peptide stretched out (no folding)
    build membrane [popc] [10] [cholesterol=30]   a lipid bilayer patch
    build micelle [sds] [60]           a ball of detergent molecules
    build assembly ferritin            the complete biological unit from the PDB

Two kinds of result, always said plainly:
  IDEALISED MODEL      built from geometry rules; nothing was measured
  FROM EXPERIMENT      a real structure, completed with the symmetry stored
                       with it in the PDB (the "biological assembly")
"""
from __future__ import annotations

import contextlib
import io
import json
import math
import os
import random
import re

from pymol import cmd

from ..core import cache, http, rcsb, report, search
from . import commands as base
from . import teach, worker

_say = base._say
_DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     "core", "data", "shapes")

LIPIDS = {"popc": "POV", "dopc": "PCW", "dope": "PEE"}
DETERGENTS = {"sds": "SDS", "ldao": "LDA"}
LIPID_FACTS = {
    "POV": ("POPC", "the most common lipid in animal cell membranes: a phosphocholine head, "
            "one straight saturated tail and one kinked unsaturated tail"),
    "PCW": ("DOPC", "a phosphocholine lipid with two kinked (unsaturated) tails; it makes "
            "membranes more fluid"),
    "PEE": ("DOPE", "a lipid with a small ethanolamine head; common on the inner side of cell "
            "membranes and in bacteria"),
    "SDS": ("SDS", "sodium dodecyl sulfate: a negatively charged detergent - the one in "
            "SDS-PAGE, where it unfolds proteins and coats them with charge"),
    "LDA": ("LDAO", "lauryldimethylamine oxide: a gentle detergent used to take membrane "
            "proteins out of their membrane without unfolding them"),
}
MAX_ASSEMBLY_ATOMS = 300000      # bigger assemblies are shown coarse, never downloaded whole
MAX_POINTS = 120000              # balls in a coarse view
MAX_DNA, MAX_PEPTIDE, MAX_MEMBRANE_NM = 120, 80, 20
_built = {}                      # object -> its card (text lines), for `about`
_shapes = [None]


# ------------------------------------------------------------------- helpers
def _quiet(function, *args, **kwargs):
    """Run a PyMOL builder without the internal noise it prints."""
    with contextlib.redirect_stderr(io.StringIO()):
        return function(*args, **kwargs)


def _new_name(stem):
    existing = set(cmd.get_names("all"))
    number = 1
    while "%s_%d" % (stem, number) in existing:
        number += 1
    return "%s_%d" % (stem, number)


def _keep(name, kind, card, aliases=()):
    """Register a built object: its kind, its card for `about`, its names."""
    base._objects[name] = ("built", kind)
    _built[name] = card
    base._last[0] = name
    for alias in aliases:
        base._remember(alias, name)


def card(name):
    """The card of a built object (for `about`), or None."""
    return _built.get(name)


def _kind_line(experimental=False):
    if experimental:
        return ("FROM EXPERIMENT  a real structure from the PDB, completed with the symmetry "
                "deposited with it")
    return ("IDEALISED MODEL  built from geometry rules - a picture of the textbook shape, "
            "not a measured structure")


def _shape_info():
    if _shapes[0] is None:
        with open(os.path.join(_DATA, "shapes.json"), "r", encoding="utf-8") as handle:
            _shapes[0] = json.load(handle)
    return _shapes[0]


def _template(code):
    """(atoms, bonds) of a bundled straight shape: head at the origin, tails -z.
    atoms = [(element, [x, y, z])], bonds = [(i, j, order)]."""
    with open(os.path.join(_DATA, code + ".mol"), "r", encoding="utf-8") as handle:
        lines = handle.read().split("\n")
    count, nbonds = int(lines[3][:3]), int(lines[3][3:6])
    atoms = [(l[31:34].strip(), [float(l[0:10]), float(l[10:20]), float(l[20:30])])
             for l in lines[4:4 + count]]
    bonds = [(int(l[0:3]) - 1, int(l[3:6]) - 1, int(l[6:9]))
             for l in lines[4 + count:4 + count + nbonds]]
    return atoms, bonds


def _rotation(axis, degrees):
    x, y, z = axis
    length = math.sqrt(x * x + y * y + z * z) or 1.0
    x, y, z = x / length, y / length, z / length
    a = math.radians(degrees)
    c, s, t = math.cos(a), math.sin(a), 1 - math.cos(a)
    return [[t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c]]


def _apply(matrix, point):
    return [sum(matrix[r][k] * point[k] for k in range(3)) for r in range(3)]


def _align(u, v):
    """Rotation matrix taking direction u onto direction v."""
    def unit(a):
        n = math.sqrt(sum(x * x for x in a))
        return [x / n for x in a]
    u, v = unit(u), unit(v)
    c = sum(u[k] * v[k] for k in range(3))
    axis = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
    s = math.sqrt(sum(x * x for x in axis))
    if s < 1e-9:
        return [[1, 0, 0], [0, 1, 0], [0, 0, 1]] if c > 0 else _rotation([1, 0, 0], 180)
    return _rotation(axis, math.degrees(math.atan2(s, c)))


def _model_from_copies(code, placements):
    """One PyMOL object made of many copies of a bundled shape.
    placements = [(matrix, shift, chain, resi, code)]: each copy's atoms are
    rotated by matrix, moved by shift, and named as residue code/resi."""
    from chempy import Atom, Bond, models
    model = models.Indexed()
    templates = {}
    for matrix, shift, chain, resi, which in placements:
        if which not in templates:
            atoms, bonds = _template(which)
            counts, names = {}, []
            totals = {}
            for element, _ in atoms:
                totals[element] = totals.get(element, 0) + 1
            for element, _ in atoms:
                counts[element] = counts.get(element, 0) + 1
                # an element that occurs once keeps its plain name (P, N, S)
                names.append(element if totals[element] == 1
                             else "%s%d" % (element, counts[element]))
            templates[which] = (atoms, bonds, names)
        atoms, bonds, names = templates[which]
        start = len(model.atom)
        for (element, xyz), name in zip(atoms, names):
            atom = Atom()
            atom.symbol, atom.name = element, name
            atom.resn, atom.resi, atom.resi_number = which, str(resi), resi
            atom.chain, atom.segi, atom.hetatm = chain, chain, 1
            moved = _apply(matrix, xyz)
            atom.coord = [moved[k] + shift[k] for k in range(3)]
            model.atom.append(atom)
        for i, j, order in bonds:
            bond = Bond()
            bond.index, bond.order = [start + i, start + j], order
            model.bond.append(bond)
    return model


def _options(words):
    """('key=value' options, the other words)."""
    options, rest = {}, []
    for word in words:
        match = re.match(r"^([a-z]+)=(.+)$", word.lower())
        if match:
            options[match.group(1)] = match.group(2)
        else:
            rest.append(word)
    return options, rest


def _finish(label):
    _say(report.LINE)


# ----------------------------------------------------------- nucleic acids
_BASE_COLOURS = (("DA+A", "green", "A"), ("DT+U", "red", "T/U"), ("DG+G", "orange", "G"),
                 ("DC+C", "skyblue", "C"))


def _nucleic(kind, words):
    options, rest = _options(words)
    sequence = "".join(rest).upper().replace("5'", "").replace("3'", "")
    sequence = re.sub(r"[^A-Z]", "", sequence)
    letters = "ACGT" if kind == "dna" else "ACGU"
    if not sequence:
        _say("Give the sequence, 5' to 3', e.g.  build %s %s" % (
            kind, "ATGCGTACGTTA" if kind == "dna" else "GGCUAGCUAC"))
        return
    if kind == "rna" and "T" in sequence and "U" not in sequence:
        sequence = sequence.replace("T", "U")
        _say("NOTE         RNA has U instead of T, so every T was read as U.")
    if kind == "dna" and "U" in sequence and "T" not in sequence:
        sequence = sequence.replace("U", "T")
        _say("NOTE         DNA has T instead of U, so every U was read as T.")
    wrong = sorted(set(sequence) - set(letters))
    if wrong:
        _say("%s only uses the letters %s; '%s' is not one of them."
             % (kind.upper(), ", ".join(letters), "".join(wrong)))
        return
    if len(sequence) > MAX_DNA:
        _say("That is %d bases; I build up to %d (longer helices are slow and hard to look at)."
             % (len(sequence), MAX_DNA))
        return
    form = (options.get("form") or "B").upper() if kind == "dna" else "A"
    if form not in ("A", "B"):
        _say("form must be A or B (Z-DNA is not available in PyMOL's builder).")
        return
    single = kind == "rna" or "single" in [w.lower() for w in rest] or \
        options.get("strands") == "1"
    name = _new_name(kind)
    base._snapshot("build %s" % kind, allow_empty=True)
    _quiet(cmd.fnab, sequence, name=name, mode=kind.upper(), form=form,
           dbl_helix=0 if single else 1)
    if not cmd.count_atoms(name):
        _say("PyMOL could not build that sequence.")
        return
    cmd.hide("everything", name)
    cmd.show("cartoon", name)
    cmd.set("cartoon_ring_mode", 3, name)
    cmd.set("cartoon_ladder_mode", 1, name)
    cmd.color("grey70", name)
    for residues, colour, _ in _BASE_COLOURS:
        cmd.color(colour, "%s and resn %s" % (name, residues))
    # the backbone of each strand in one plain colour, so the bases stand out
    for chain, colour in (("A", "grey50"), ("B", "slate")):
        cmd.set("cartoon_color", colour, "%s and chain %s" % (name, chain))
    cmd.orient(name)
    teach.note("fnab %s, name=%s, mode=%s, form=%s, dbl_helix=%d"
               % (sequence, name, kind.upper(), form, 0 if single else 1),
               "PyMOL's own nucleic-acid builder")
    pairs = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}
    length = len(sequence)
    rise, per_turn, width = (3.4, 10.5, 20) if form == "B" else (2.6, 11.0, 23)
    lines = ["BUILT        %s %s, %d %s: %s" % (
        "a single" if single else "a double-stranded", kind.upper(), length,
        "bases" if single else "base pairs", name)]
    lines.append(_kind_line())
    if not single:
        lines.append("SEQUENCE     5' %s 3'" % sequence)
        lines.append("             %s" % ("   " + "|" * length))
        lines.append("             3' %s 5'   (the partner strand, paired base by base)"
                     % "".join(pairs[b] for b in sequence))
    else:
        lines.append("SEQUENCE     5' %s 3'" % sequence)
    lines.append("ON SCREEN    the sugar-phosphate backbone as a ribbon%s; each base as a "
                 "flat ring: green = A, red = %s, orange = G, blue = C"
                 % (" (grey and blue: one per strand)" if not single else "",
                    "U" if kind == "rna" else "T"))
    if form == "B" and not single:
        lines.append("WHAT IT SHOWS  B-DNA, the form in cells: a right-handed helix about "
                     "%d A wide, %.1f base pairs per turn, %.1f A between base pairs - here "
                     "about %.0f A long and %.1f turns. Look for the two grooves: the wide "
                     "MAJOR groove is where most proteins read the sequence."
                     % (width, per_turn, rise, length * rise, length / per_turn))
    elif form == "A" and not single:
        lines.append("WHAT IT SHOWS  A-form DNA: wider and shorter than B-DNA, about 11 "
                     "base pairs per turn, bases tilted. Double-stranded RNA takes this "
                     "shape, and DNA does too when it dries out.")
    else:
        lines.append("WHAT IT SHOWS  one strand twisted as it would be inside an A-form "
                     "double helix. A real single strand folds back on itself wherever "
                     "bases can pair (hairpins, loops); that folding is not modelled here.")
        if kind == "rna":
            lines.append("NOTE         PyMOL cannot build double-stranded RNA. Its shape is "
                         "close to  build dna <sequence> form=A .")
    lines.append("TRY          scout %s dna  |  measure %s A/1/P A/%d/P  |  viewas %s sticks"
                 % (name, name, length, name))
    for line in lines:
        _say(line)
    _keep(name, "%s %s" % ("single-stranded" if single else "double-stranded", kind.upper()),
          lines, aliases=(kind, "%s helix" % kind))


# --------------------------------------------------------------- peptides
_AMINO = set("ACDEFGHIKLMNPQRSTVWY")


def _sequence(words, what):
    sequence = re.sub(r"[^A-Za-z]", "", "".join(words)).upper()
    if not sequence:
        _say("Give the sequence in one-letter code, e.g.  build %s AKLEELLKKLA" % what)
        return None
    wrong = sorted(set(sequence) - _AMINO)
    if wrong:
        _say("'%s' is not one of the 20 amino-acid letters (ACDEFGHIKLMNPQRSTVWY)."
             % "".join(wrong))
        return None
    if len(sequence) > MAX_PEPTIDE:
        _say("That is %d residues; I build up to %d." % (len(sequence), MAX_PEPTIDE))
        return None
    return sequence


def _peptide(kind, words):
    options, rest = _options(words)
    if kind == "sheet":
        return _sheet(rest)
    sequence = _sequence(rest, kind)
    if sequence is None:
        return
    name = _new_name(kind)
    base._snapshot("build %s" % kind, allow_empty=True)
    ss = {"helix": 1, "strand": 2, "peptide": 0}[kind]
    _quiet(cmd.fab, sequence, name, ss=ss)
    cmd.hide("everything", name)
    cmd.show("cartoon", name)
    cmd.show("sticks", "%s and not hydro and (sidechain or name CA)" % name)
    cmd.set("cartoon_transparency", 0.2, name)
    cmd.util.cbag(name)
    cmd.color("grey80", "%s and name CA" % name)
    lines = ["BUILT        a %d-residue peptide as %s: %s" % (
        len(sequence), {"helix": "an alpha-helix", "strand": "an extended beta-strand",
                        "peptide": "a stretched-out chain"}[kind], name), _kind_line(),
        "SEQUENCE     N-%s-C" % sequence]
    if kind == "helix":
        cmd.distance(name + "_hbonds", "%s and name O" % name, "%s and name N" % name,
                     3.5, mode=2)
        cmd.color("yellow", name + "_hbonds")
        cmd.hide("labels", name + "_hbonds")
        lines += ["ON SCREEN    the helix as a ribbon, side chains as sticks pointing outward, "
                  "yellow dashes = the hydrogen bonds that hold the helix together",
                  "WHAT IT SHOWS  3.6 residues per turn and 1.5 A of rise per residue: "
                  "%d residues make %.1f turns and about %.0f A. Each C=O bonds to the N-H "
                  "four residues further on. Whether a real sequence forms a helix depends "
                  "on its amino acids - here the shape is imposed."
                  % (len(sequence), len(sequence) / 3.6, len(sequence) * 1.5)]
    elif kind == "strand":
        lines += ["ON SCREEN    the chain as a flat arrow, side chains as sticks",
                  "WHAT IT SHOWS  a beta-strand: almost fully stretched, about 3.3 A per "
                  "residue, side chains alternating up and down. A strand alone is not "
                  "stable - it needs a partner strand to bond with. Try  build sheet "
                  "%s %s" % (sequence, sequence[::-1])]
    else:
        lines += ["ON SCREEN    the chain as a thin ribbon, side chains as sticks",
                  "WHAT IT SHOWS  the primary structure only: the chain fully stretched, "
                  "nothing folded. Every peptide bond is flat; the chain bends only at the "
                  "alpha-carbons. See it fold into a shape with  build helix %s" % sequence]
    cmd.orient(name)
    teach.note("fab %s, %s, ss=%d" % (sequence, name, ss), "PyMOL's peptide builder: "
               "ss=1 helix, ss=2 strand")
    lines.append("TRY          scout %s hydrophobic  |  sequence %s  |  viewas %s sticks"
                 % (name, name, name))
    for line in lines:
        _say(line)
    _keep(name, "peptide " + kind, lines, aliases=(kind, "peptide"))


def _ca(obj):
    found = []
    cmd.iterate_state(1, "%s and name CA" % obj, "found.append((int(resi), list((x,y,z))))",
                      space={"found": found})
    return [xyz for _, xyz in sorted(found)]


def _sheet(words):
    strands = [w for w in words if re.sub(r"[^A-Za-z]", "", w)]
    sequences = []
    for word in strands:
        sequence = _sequence([word], "sheet")
        if sequence is None:
            return
        sequences.append(sequence)
    if not sequences:
        _say("Give two or more strand sequences, e.g.  build sheet KVTVKVE TVKVTVE")
        return
    note = ""
    if len(sequences) == 1:
        sequences = sequences * 2
        note = "One sequence was given, so the sheet has two copies of it."
    if len(sequences) > 6:
        _say("I build sheets of up to 6 strands.")
        return
    name = _new_name("sheet")
    base._snapshot("build sheet", allow_empty=True)
    parts = []
    for index, sequence in enumerate(sequences):
        part = "%s_part%d" % (name, index)
        _quiet(cmd.fab, sequence, part, ss=2, chain=chr(65 + index))
        parts.append(part)
    # lay each strand antiparallel next to the previous one, sliding it until
    # the most backbone hydrogen bonds (N...O 2.6-3.3 A) form
    for index in range(1, len(parts)):
        previous, part = parts[index - 1], parts[index]
        ca = _ca(previous)
        u = [ca[-1][k] - ca[0][k] for k in range(3)]
        centre = [sum(p[k] for p in ca) / len(ca) for k in range(3)]
        cb = []
        cmd.iterate_state(1, "%s and name CB" % previous, "cb.append((int(resi),(x,y,z)))",
                          space={"cb": cb})
        cas = dict((i + 1, p) for i, p in enumerate(ca))
        normal = [0.0, 0.0, 0.0]
        for resi, xyz in cb:
            sign = 1 if resi % 2 else -1
            for k in range(3):
                normal[k] += sign * (xyz[k] - cas.get(resi, xyz)[k])
        un = math.sqrt(sum(x * x for x in u))
        u = [x / un for x in u]
        along = sum(normal[k] * u[k] for k in range(3))
        normal = [normal[k] - along * u[k] for k in range(3)]
        nn = math.sqrt(sum(x * x for x in normal)) or 1.0
        normal = [x / nn for x in normal]
        w = [normal[1] * u[2] - normal[2] * u[1], normal[2] * u[0] - normal[0] * u[2],
             normal[0] * u[1] - normal[1] * u[0]]
        # start from the previous strand's position, turned round (antiparallel)
        mine = _ca(part)
        mc = [sum(p[k] for p in mine) / len(mine) for k in range(3)]
        cmd.translate([-mc[k] for k in range(3)], part, camera=0)
        cmd.fit("%s and name CA" % part, "%s and name CA" % previous,
                matchmaker=-1) if len(mine) == len(ca) else None
        mine = _ca(part)
        mc = [sum(p[k] for p in mine) / len(mine) for k in range(3)]
        cmd.translate([centre[k] - mc[k] for k in range(3)], part, camera=0)
        cmd.rotate(normal, 180, part, camera=0, origin=centre)
        best = None
        for side in (1, -1):
            for slide in [x * 0.5 for x in range(-8, 9)]:
                for gap in (4.4, 4.8, 5.2):
                    shift = [side * gap * w[k] + slide * u[k] for k in range(3)]
                    cmd.translate(shift, part, camera=0)
                    bonds = cmd.count_atoms("(%s and name N) within 3.3 of (%s and name O)"
                                            % (part, previous)) + \
                        cmd.count_atoms("(%s and name O) within 3.3 of (%s and name N)"
                                        % (part, previous))
                    placed = " or ".join(parts[:index])      # every strand already laid
                    bumps = cmd.count_atoms("(%s and not hydro) within 3.0 of ((%s) and not "
                                            "hydro)" % (part, placed))
                    score = bonds - 3 * bumps
                    if index >= 2:            # keep the sheet flat: away from the strand
                        far_side = cmd.count_atoms("(%s and name CA) within 7 of (%s and name "
                                                   "CA)" % (part, parts[index - 2]))
                        score -= 5 * far_side
                    if best is None or score > best[0]:
                        best = (score, shift)
                    cmd.translate([-x for x in shift], part, camera=0)
        cmd.translate(best[1], part, camera=0)
    cmd.create(name, " or ".join(parts))
    for part in parts:
        cmd.delete(part)
    cmd.hide("everything", name)
    cmd.show("cartoon", name)
    cmd.show("sticks", "%s and not hydro and (sidechain or name CA)" % name)
    cmd.util.cbc(name)
    cmd.util.cnc(name)
    cmd.distance(name + "_hbonds", "%s and name N" % name, "%s and name O" % name, 3.3, mode=2)
    cmd.color("yellow", name + "_hbonds")
    cmd.hide("labels", name + "_hbonds")
    cmd.dss(name)
    cmd.orient(name)
    hbonds = 0
    for index in range(1, len(sequences)):
        hbonds += cmd.count_atoms("(%s and chain %s and name N) within 3.3 of (%s and chain %s "
                                  "and name O)" % (name, chr(65 + index), name, chr(64 + index)))
    lines = ["BUILT        a %d-stranded beta-sheet: %s" % (len(sequences), name), _kind_line(),
             "STRANDS      " + "   ".join("%s: %s" % (chr(65 + i), s)
                                          for i, s in enumerate(sequences)),
             "ON SCREEN    each strand as a flat arrow in its own colour, side chains as sticks, "
             "yellow dashes = hydrogen bonds between neighbouring strands (%d placed)" % hbonds,
             "WHAT IT SHOWS  an antiparallel sheet: neighbouring strands run in opposite "
             "directions and are held together by backbone N-H...O=C hydrogen bonds; side "
             "chains alternate above and below the sheet. The strands were slid along each "
             "other to the position with the most hydrogen bonds - a simple rule, not a "
             "prediction of how these sequences would fold."]
    if note:
        lines.append("NOTE         " + note)
    lines.append("TRY          scout %s hydrophobic  (one face of a sheet is often oily)"
                 % name)
    for line in lines:
        _say(line)
    _keep(name, "beta-sheet", lines, aliases=("sheet", "beta sheet"))


# --------------------------------------------------------------- membranes
def _membrane(words):
    options, rest = _options(words)
    lipid, size = "popc", 8.0
    for word in rest:
        low = word.lower().rstrip(",")
        sized = re.match(r"^(\d+(?:\.\d+)?)(?:\s*nm)?(?:x(\d+(?:\.\d+)?))?$", low)
        if low in LIPIDS or low.upper() in LIPIDS.values():
            lipid = low if low in LIPIDS else [k for k, v in LIPIDS.items()
                                                if v == low.upper()][0]
        elif sized:
            size = float(sized.group(1))
        elif low in ("membrane", "bilayer", "lipid", "lipids", "nm", "of", "a", "with", "and"):
            continue
        elif low in ("cholesterol", "chol", "clr"):
            options.setdefault("cholesterol", "30")
        else:
            _say("I can build membranes of %s. '%s' is not one of them."
                 % (", ".join(sorted(LIPIDS)), word))
            return
    size = float(options.get("size", size))
    if not 2 <= size <= MAX_MEMBRANE_NM:
        _say("The patch can be 2 to %d nm wide (you asked for %g)." % (MAX_MEMBRANE_NM, size))
        return
    try:
        chol = float(str(options.get("cholesterol", "0")).rstrip("%"))
    except ValueError:
        _say("cholesterol= takes a percentage, e.g. cholesterol=30")
        return
    if not 0 <= chol <= 50:
        _say("cholesterol= can be 0 to 50 (%).")
        return
    code = LIPIDS[lipid]
    shapes = _shape_info()
    info = shapes[code]
    sx, sy = info["lattice"]
    length = info["length"]
    nx = max(1, int(round(size * 10 / sx)))
    ny = max(1, int(round(size * 10 / sy)))
    gap = 3.0                             # between the tail ends of the two leaflets
    flip = _rotation([1, 0, 0], 180)      # turns a lipid upside down
    upright = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    chooser = random.Random(7)            # the same membrane every time
    placements, counts = [], {"U": 0, "L": 0}
    sterols = 0
    for chain, matrix, z, offset in (("U", upright, length + gap / 2, 0.0),
                                     ("L", flip, -(length + gap / 2), sy / 2)):
        resi = 0
        for j in range(ny):
            for i in range(nx):
                x = i * sx + (sx / 2 if j % 2 else 0.0) - nx * sx / 2
                y = j * sy + offset - ny * sy / 2
                resi += 1
                if chol and chooser.random() * 100 < chol:
                    # cholesterol: its OH sits at the level of the lipids' ester groups
                    drop = 7.0 if chain == "U" else -7.0
                    placements.append((matrix, [x, y, z - drop], chain, resi, "CLR"))
                    sterols += 1
                else:
                    placements.append((matrix, [x, y, z], chain, resi, code))
                counts[chain] += 1
    name = _new_name("membrane")
    base._snapshot("build membrane", allow_empty=True)
    cmd.load_model(_model_from_copies(code, placements), name)
    cmd.hide("everything", name)
    cmd.show("sticks", name)
    cmd.set("stick_radius", 0.18, name)
    cmd.show("spheres", "%s and elem P+N" % name)
    cmd.set("sphere_scale", 0.9, name)
    cmd.color("grey70", "%s and elem C" % name)
    cmd.color("red", "%s and elem O" % name)
    cmd.color("blue", "%s and elem N" % name)
    cmd.color("orange", "%s and elem P" % name)
    cmd.color("yellow", "%s and resn CLR and elem C" % name)
    cmd.orient(name)
    cmd.turn("x", -70)
    total = counts["U"] + counts["L"]
    lipid_name, about = LIPID_FACTS[code]
    span = 2 * length + gap
    lines = ["BUILT        a %s bilayer patch, %.0f x %.0f nm: %d molecules (%d per leaflet)%s"
             % (lipid_name, nx * sx / 10, ny * sy / 10, total, total // 2,
                ", %d of them cholesterol" % sterols if sterols else ""),
             _kind_line(),
             "ON SCREEN    two layers of lipids, tails facing each other: orange balls = "
             "phosphate, blue balls = nitrogen of the head group, red = oxygen, grey sticks "
             "= fatty-acid tails%s" % (", yellow sticks = cholesterol" if sterols else ""),
             "WHAT IT SHOWS  why membranes form: the heads like water and face out, the "
             "oily tails hide from water in the middle. %s." % (about[0].upper() + about[1:]),
             "HOW REAL    the tails are drawn fully stretched, as in a cold, gel-like "
             "membrane, so this bilayer is about %.0f A thick between the phosphates. In a "
             "living cell's fluid membrane the tails wriggle and it is thinner, about 38-40 A."
             % (span + 0)]
    if sterols:
        lines.append("CHOLESTEROL  its small OH group sits among the head groups and its "
                     "rigid rings lie along the tails: it stiffens the membrane and makes "
                     "it less leaky.")
    lines.append("TRY          measure %s U/1/P L/1/P  (thickness)  |  viewas %s spheres  "
                 "(the packed volume)  |  choose phosphorus atoms in %s" % (name, name, name))
    for line in lines:
        _say(line)
    teach.note("load_model <copies of %s>, %s" % (code, name),
               "the lipid is copied onto a grid; the lower layer is a turned-over copy")
    _keep(name, "membrane", lines, aliases=("membrane", "bilayer", lipid + " membrane"))


# ----------------------------------------------------------------- micelles
def _micelle(words):
    options, rest = _options(words)
    detergent, count = "sds", None
    for word in rest:
        low = word.lower()
        if low in DETERGENTS or low.upper() in DETERGENTS.values():
            detergent = low if low in DETERGENTS else [k for k, v in DETERGENTS.items()
                                                       if v == low.upper()][0]
        elif low.isdigit():
            count = int(low)
        elif low in ("micelle", "of", "a", "with", "molecules"):
            continue
        else:
            _say("I can build micelles of %s. '%s' is not one of them."
                 % (", ".join(sorted(DETERGENTS)), word))
            return
    count = int(options.get("count", count or (60 if detergent == "sds" else 75)))
    if not 10 <= count <= 200:
        _say("A micelle can have 10 to 200 molecules (you asked for %d)." % count)
        return
    code = DETERGENTS[detergent]
    length = _shape_info()[code]["length"]
    # straight tails meet near the centre: keep their ends at least ~4 A apart
    inner = max(3.0, math.sqrt(count * 14.0 / (4 * math.pi)))
    radius = inner + length
    golden = math.pi * (3 - math.sqrt(5))
    chooser = random.Random(11)
    placements = []
    for i in range(count):
        z = 1 - 2 * (i + 0.5) / count
        r = math.sqrt(1 - z * z)
        out = [r * math.cos(golden * i), r * math.sin(golden * i), z]
        spin = _rotation([0, 0, 1], chooser.uniform(0, 360))
        turn = _align([0, 0, 1], out)              # head (+z of the shape) points outward
        matrix = [[sum(turn[a][k] * spin[k][b] for k in range(3)) for b in range(3)]
                  for a in range(3)]
        placements.append((matrix, [radius * x for x in out], "M", i + 1, code))
    name = _new_name("micelle")
    base._snapshot("build micelle", allow_empty=True)
    cmd.load_model(_model_from_copies(code, placements), name)
    cmd.hide("everything", name)
    cmd.show("spheres", name)
    cmd.color("grey70", "%s and elem C" % name)
    cmd.color("red", "%s and elem O" % name)
    cmd.color("yellow", "%s and elem S" % name)
    cmd.color("blue", "%s and elem N" % name)
    cmd.orient(name)
    cmd.clip("slab", 200)
    name_text, about = LIPID_FACTS[code]
    lines = ["BUILT        a micelle of %d %s molecules, about %.0f A across: %s"
             % (count, name_text, 2 * radius + 4, name), _kind_line(),
             "ON SCREEN    every atom as a ball: grey = the carbon tails, %s = the head groups "
             "on the outside" % ("yellow and red (sulfate)" if code == "SDS"
                                 else "blue and red (amine oxide)"),
             "WHAT IT SHOWS  detergent molecules in water hide their oily tails in the "
             "middle and keep their heads outside: a ball with a dry core. %s."
             % (about[0].upper() + about[1:]),
             "HOW REAL    the tails are drawn straight, so the very centre is emptier than "
             "in a real micelle, whose tails are bent and fill the core like a liquid. Real "
             "micelles also change shape and swap molecules all the time.",
             "TRY          viewas %s sticks  then  clip slab, 10  (a slice through the middle "
             "shows the core)" % name]
    for line in lines:
        _say(line)
    _keep(name, "micelle", lines, aliases=("micelle", detergent + " micelle"))


# --------------------------------------------------------------- assemblies
_ASSEMBLY_QUERY = """{ entry(entry_id: "%s") {
  rcsb_id struct { title }
  rcsb_entry_info { deposited_atom_count deposited_polymer_entity_instance_count }
  assemblies { rcsb_id
    rcsb_assembly_info { atom_count polymer_entity_instance_count }
    pdbx_struct_assembly { oligomeric_details oligomeric_count details method_details }
  } } }"""
_OPERATOR_QUERY = """{ assembly(entry_id: "%s", assembly_id: "%s") {
  pdbx_struct_assembly_gen { asym_id_list oper_expression }
  pdbx_struct_oper_list { id matrix_1_1 matrix_1_2 matrix_1_3 matrix_2_1 matrix_2_2
    matrix_2_3 matrix_3_1 matrix_3_2 matrix_3_3 vector_1 vector_2 vector_3 } } }"""


def _operators(expression):
    """'(1-60)(61-88)' -> [[ids of group 1], [ids of group 2]] (a product)."""
    groups = re.findall(r"\(([^)]*)\)", expression) or [expression]
    out = []
    for group in groups:
        ids = []
        for part in group.split(","):
            part = part.strip()
            if re.match(r"^\d+-\d+$", part):
                a, b = part.split("-")
                ids.extend(str(i) for i in range(int(a), int(b) + 1))
            elif part:
                ids.append(part)
        out.append(ids)
    return out


def assembly_size(pdb_id):
    """(chains in the biological assembly, chains in the file, oligomer words), or
    None.  One small cached request - used by grab to say "this is only part"."""
    try:
        data = http.request_json(rcsb.GRAPHQL_URL, {"query": _ASSEMBLY_QUERY % pdb_id},
                                 ttl=7 * cache.DAY)
    except http.NetworkError:
        return None
    entry = ((data or {}).get("data") or {}).get("entry") or {}
    assemblies = sorted(entry.get("assemblies") or [], key=lambda a: a["rcsb_id"])
    if not assemblies:
        return None
    info = assemblies[0].get("rcsb_assembly_info") or {}
    details = assemblies[0].get("pdbx_struct_assembly") or {}
    deposited = (entry.get("rcsb_entry_info") or {}).get(
        "deposited_polymer_entity_instance_count") or 0
    return (info.get("polymer_entity_instance_count") or 0, deposited,
            details.get("oligomeric_details") or "", info.get("atom_count") or 0,
            assemblies[0]["rcsb_id"].split("-")[-1])


def whole_molecule(pdb_id):
    """Worker half, used whenever a PDB entry is loaded: when the file holds only
    part of the working molecule (one chain of a tetramer), download the PDB's
    assembly file instead.  Returns {"assembly": (...) or None, "path": path of the
    assembly file or None}."""
    size = assembly_size(pdb_id)
    out = {"assembly": size, "path": None}
    # load the assembly when the file holds less (one chain of a tetramer) or more
    # (two copies packed side by side in the crystal) than the working molecule
    if size and size[0] > 0 and size[1] > 0 and size[0] != size[1] and \
            0 < size[3] <= MAX_ASSEMBLY_ATOMS and size[0] <= 62:
        target = os.path.join(cache.structures_dir(), "%s-assembly%s.cif" % (pdb_id, size[4]))
        try:
            out["path"] = http.download("https://files.rcsb.org/download/%s-assembly%s.cif"
                                        % (pdb_id, size[4]), target)
        except http.NetworkError:
            out["path"] = None           # fall back to the file as deposited
    return out


def _id_variants(pdb_id):
    """8Y80 -> 8Y8O and back: the letter O and the digit 0 look alike."""
    out = []
    for index in range(1, len(pdb_id)):
        swap = {"0": "O", "O": "0"}.get(pdb_id[index])
        if swap:
            out.append(pdb_id[:index] + swap + pdb_id[index + 1:])
    return out


def _assembly_plan(pdb_id):
    """Worker half: what the biological assembly of an entry is - trying the
    letter O for the digit 0 (and back) if the entry does not exist."""
    plan = _assembly_plan_exact(pdb_id)
    if (plan.get("error") or "").startswith("The PDB has no entry"):
        for other in _id_variants(pdb_id):
            alternative = _assembly_plan_exact(other)
            if not alternative.get("error"):
                alternative["note"] = ("There is no entry %s; I used %s (the letter O and "
                                       "the digit 0 look alike)." % (pdb_id, other))
                return alternative
    return plan


def _assembly_plan_exact(pdb_id):
    """Worker half: what the biological assembly of an entry is, and how to show it."""
    data = http.request_json(rcsb.GRAPHQL_URL, {"query": _ASSEMBLY_QUERY % pdb_id},
                             ttl=7 * cache.DAY)
    entry = ((data or {}).get("data") or {}).get("entry")
    if not entry:
        return {"error": "The PDB has no entry %s." % pdb_id}
    assemblies = sorted(entry.get("assemblies") or [], key=lambda a: a["rcsb_id"])
    if not assemblies:
        return {"error": "%s has no biological assembly on record." % pdb_id}
    chosen = assemblies[0]
    info = chosen["rcsb_assembly_info"] or {}
    details = chosen["pdbx_struct_assembly"] or {}
    plan = {"pdb": pdb_id, "title": (entry.get("struct") or {}).get("title") or "",
            "assembly": chosen["rcsb_id"].split("-")[-1],
            "atoms": info.get("atom_count") or 0,
            "chains": info.get("polymer_entity_instance_count") or 0,
            "deposited_atoms": (entry.get("rcsb_entry_info") or {}).get(
                "deposited_atom_count") or 0,
            "deposited_chains": (entry.get("rcsb_entry_info") or {}).get(
                "deposited_polymer_entity_instance_count") or 0,
            "oligomer": details.get("oligomeric_details") or "",
            "how": details.get("details") or "", "others": len(assemblies) - 1,
            "error": None}
    if plan["chains"] <= plan["deposited_chains"] and plan["atoms"] <= plan["deposited_atoms"] * 1.05:
        plan["mode"] = "same"
        plan["path"] = base._download_pdb(pdb_id)
    elif plan["atoms"] <= MAX_ASSEMBLY_ATOMS:
        plan["mode"] = "full"
        path = os.path.join(cache.structures_dir(), "%s-assembly%s.cif"
                            % (pdb_id, plan["assembly"]))
        plan["path"] = http.download("https://files.rcsb.org/download/%s-assembly%s.cif"
                                     % (pdb_id, plan["assembly"]), path)
    else:
        plan["mode"] = "coarse"
        plan["path"] = base._download_pdb(pdb_id)
        ops = http.request_json(rcsb.GRAPHQL_URL, {"query": _OPERATOR_QUERY
                                                   % (pdb_id, plan["assembly"])},
                                ttl=7 * cache.DAY)
        plan["operators"] = ((ops or {}).get("data") or {}).get("assembly")
    return plan


def _coarse(plan, name):
    """A huge assembly as balls: every k-th amino acid of every copy."""
    temp = "_prospector_asym"
    cmd.load(plan["path"], temp)
    data = plan["operators"] or {}
    matrices = {}
    for op in data.get("pdbx_struct_oper_list") or []:
        m = [[op["matrix_%d_%d" % (r, c)] for c in (1, 2, 3)] for r in (1, 2, 3)]
        matrices[op["id"]] = (m, [op["vector_1"], op["vector_2"], op["vector_3"]])
    jobs = []
    for gen in data.get("pdbx_struct_assembly_gen") or []:
        groups = _operators(gen["oper_expression"])
        combos = [[]]
        for group in groups:
            combos = [c + [g] for c in combos for g in group]
        listed = gen["asym_id_list"]
        if isinstance(listed, str):
            listed = listed.split(",")
        asyms = [a.strip() for a in listed if a and a.strip()]
        for combo in combos:
            jobs.append((combo, asyms))
    points = {}
    for asym in set(a for _, asyms in jobs for a in asyms):
        found = []
        cmd.iterate_state(1, "%s and polymer and name CA and segi %s" % (temp, asym),
                          "found.append((resn, int(resi) if resi.lstrip('-').isdigit() else 0, "
                          "(x, y, z)))", space={"found": found})
        points[asym] = found
    total = sum(len(points.get(a, [])) for _, asyms in jobs for a in asyms)
    step = max(1, int(math.ceil(total / float(MAX_POINTS))))
    from chempy import Atom, models
    model = models.Indexed()
    for copy, (combo, asyms) in enumerate(jobs):
        # mmCIF: "(A)(B)" means apply B first, then A
        matrix, shift = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0]
        for op_id in reversed(combo):
            m, v = matrices.get(op_id, ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0]))
            shift = [sum(m[r][k] * shift[k] for k in range(3)) + v[r] for r in range(3)]
            matrix = [[sum(m[r][k] * matrix[k][c] for k in range(3)) for c in range(3)]
                      for r in range(3)]
        for asym in asyms:
            for n, (resn, resi, xyz) in enumerate(points.get(asym, [])):
                if n % step:
                    continue
                atom = Atom()
                atom.symbol, atom.name, atom.resn = "C", "CA", resn
                atom.resi, atom.resi_number = str(resi), resi
                atom.chain = atom.segi = "%s%d" % (asym, copy + 1)
                atom.b = float(copy)
                atom.coord = [sum(matrix[r][k] * xyz[k] for k in range(3)) + shift[r]
                              for r in range(3)]
                model.atom.append(atom)
    cmd.delete(temp)
    cmd.load_model(model, name)
    cmd.hide("everything", name)
    cmd.show("spheres", name)
    cmd.set("sphere_scale", 1.2 + 0.4 * min(step, 10), name)
    cmd.spectrum("b", "rainbow", name)
    return len(jobs) * len(set(a for _, asyms in jobs for a in asyms)), step, len(model.atom)


def _plain_chain_names(obj):
    """Assembly files call the copies A, B, A-2, B-2 ...  The first copy keeps its
    chain names (so UniProt annotations still match); the other copies get the next
    free letters (C, D ...), so chains look - and compare - like a normal entry.
    Returns [(old, new)], empty when nothing was renamed."""
    order = []
    cmd.iterate(obj, "order.append(chain) if not order or order[-1] != chain else None",
                space={"order": order})
    chains = []
    for chain in order:
        if chain not in chains:
            chains.append(chain)
    if not any("-" in c for c in chains) or len(chains) > 62:
        return []
    bases = []
    for chain in chains:
        if chain.split("-")[0] not in bases:
            bases.append(chain.split("-")[0])

    def copy_order(chain):            # copy 1 (A, B), then copy 2 (A-2, B-2) ...
        base_name, _, number = chain.partition("-")
        return (int(number) if number.isdigit() else 1, bases.index(base_name))
    chains.sort(key=copy_order)
    pool = [chr(c) for c in range(ord("A"), ord("Z") + 1)] + \
        [chr(c) for c in range(ord("a"), ord("z") + 1)] + [str(d) for d in range(10)]
    keep = [c for c in chains if "-" not in c]
    free = [c for c in pool if c not in keep]
    mapping = dict((c, c) for c in keep)
    for chain in chains:
        if chain not in mapping:
            if not free:
                return []
            mapping[chain] = free.pop(0)
    cmd.alter(obj, "chain = segi = mapping.get(chain, chain)", space={"mapping": mapping})
    cmd.sort(obj)
    return [(old, mapping[old]) for old in chains if old != mapping[old]]


def _assembly(words, label):
    text = " ".join(w for w in words if w.lower() not in ("of", "the", "a", "whole", "complete",
                                                           "full", "assembly"))
    state = {}

    def before():
        obj = base._known_object(text) if text else (base._last[0]
                                                    if base._last[0] in base._molecules()
                                                    else None)
        identity = base._identify(obj) if obj else None
        if identity and identity[0] == "pdb":
            state["pdb"] = identity[1]
        elif text and rcsb.looks_like_pdb_id(text):
            state["pdb"] = text.upper()
            # 8Y80 typed while 8Y8O is on screen: the letter O and the digit 0 look
            # alike, and both can be real entries - use the one on screen
            loaded = set(n.upper() for n in cmd.get_names("objects"))
            if state["pdb"] not in loaded:
                for other in _id_variants(state["pdb"]):
                    if other in loaded:
                        state["note"] = ("%s is on screen, so I used it (you typed %s: "
                                         "the letter O and the digit 0 look alike)."
                                         % (other, state["pdb"]))
                        state["pdb"] = other
                        break
        elif text:
            state["query"] = text
        else:
            base._header(label)
            _say("Name a protein or a PDB ID, e.g.  build assembly ferritin  or  build "
                 "assembly 1AON")
            _say(report.LINE)
            return False
        return True

    def job():
        try:
            if "query" in state:
                result = search.find(state["query"] + " limit=3")
                if result.kind == "ambiguous":
                    choices = result.resolution.choices
                    lines = ["'%s' could mean several proteins, so I did not guess. Pick one:"
                             % state["query"]]
                    for number, group in enumerate(choices[:5], 1):
                        lines.append("  %d. build assembly %-12s %s" % (
                            number, group.suggested_query(choices),
                            report._cut(group.label(), 40)))
                    lines.append("   (or just type its number: 1, 2 ...)")
                    base._choices[:] = ["build assembly %s" % g.suggested_query(choices)
                                        for g in choices[:5]]
                    base._choice_labels[:] = [g.label(1) for g in choices[:5]]
                    return {"error": "\n".join(lines)}
                if result.kind != "results":
                    return {"error": "'%s' did not lead to a PDB structure. %s"
                                     % (state["query"], result.message or
                                        "Try its PDB ID, e.g.  build assembly 2FHA")}
                state["pdb"] = result.top.pdb_id
            return _assembly_plan(state["pdb"])
        except http.NetworkError as error:
            return {"error": base._network_message(error)}

    def done(plan, error):
        if error:
            return base._unexpected(error)
        base._header(label)
        _say("")
        if plan.get("error"):
            _say(plan["error"])
            _say(report.LINE)
            return
        for note in (state.get("note"), plan.get("note")):
            if note:
                _say("Note: " + note)
        pdb_id = plan["pdb"]
        name = "%s_assembly" % pdb_id
        if pdb_id in base._assembled and pdb_id in cmd.get_names("objects"):
            _say("ALREADY WHOLE  %s on screen already is the complete molecule (%s, %d "
                 "chains): it was loaded from the PDB's assembly data."
                 % (pdb_id, plan["oligomer"] or "assembly", plan["chains"]))
            cmd.zoom(pdb_id)
            _say(report.LINE)
            return
        base._snapshot(label, allow_empty=True)
        if name in cmd.get_names("objects"):
            cmd.delete(name)
        lines = []
        if plan["mode"] == "same":
            if pdb_id not in cmd.get_names("objects"):
                cmd.load(plan["path"], pdb_id)
            base._objects[pdb_id] = ("pdb", pdb_id)
            name = pdb_id
            lines.append("ALREADY WHOLE  %s - %s" % (pdb_id, report._cut(plan["title"], 50)))
            lines.append("             The deposited structure already is the whole biological "
                         "unit (%s, %d chain%s): nothing had to be added."
                         % (plan["oligomer"] or "as deposited", plan["chains"],
                            "" if plan["chains"] == 1 else "s"))
        elif plan["mode"] == "full":
            cmd.load(plan["path"], name)
            renamed = _plain_chain_names(name)
            cmd.hide("everything", name)
            cmd.show("cartoon", name)
            cmd.show("sticks", "%s and organic" % name)
            cmd.show("spheres", "%s and inorganic and not solvent" % name)
            cmd.util.cbc(name)
            base._objects[name] = ("pdb", pdb_id)
            lines.append("BUILT        the biological assembly of %s: %s"
                         % (pdb_id, report._cut(plan["title"], 46)))
            lines.append(_kind_line(experimental=True))
            lines.append("             As deposited: %d chain%s. Complete: %d chains (%s), "
                         "%d atoms." % (plan["deposited_chains"],
                                        "" if plan["deposited_chains"] == 1 else "s",
                                        plan["chains"], plan["oligomer"] or "?", plan["atoms"]))
            if renamed:
                lines.append("CHAINS       named %s, like a normal PDB entry (the file calls the "
                             "copies %s)" % (", ".join(new for _, new in renamed),
                                             ", ".join(old for old, _ in renamed)))
            lines.append("ON SCREEN    every chain as a ribbon in its own colour; bound "
                         "molecules as sticks, metal ions as balls")
        else:
            copies, step, balls = _coarse(plan, name)
            base._objects[name] = ("pdb", pdb_id)
            lines.append("BUILT        a COARSE view of the biological assembly of %s: %s"
                         % (pdb_id, report._cut(plan["title"], 36)))
            lines.append(_kind_line(experimental=True))
            lines.append("             The complete assembly has %s atoms in %d chains (%s) - "
                         "far too many to load. So every %s amino acid of every copy is shown "
                         "as one ball: %d balls."
                         % ("{:,}".format(plan["atoms"]), plan["chains"],
                            plan["oligomer"] or "?",
                            "%d%s" % (step, {1: "st", 2: "nd", 3: "rd"}.get(
                                step if step < 20 else step % 10, "th")) if step > 1 else "",
                            balls))
            lines.append("ON SCREEN    one colour per copy of the deposited unit, blue to red, "
                         "so you can see how the copies tile the whole")
        lines.append("WHAT IT SHOWS  many proteins only work as a team of chains. A PDB entry "
                     "often stores just one repeating piece; the assembly is the whole object "
                     "as it exists in the cell, rebuilt from that piece and its symmetry.")
        if plan["others"]:
            lines.append("NOTE         this entry lists %d other possible assembl%s; this is the "
                         "first (usually the one the authors describe)."
                         % (plan["others"], "y" if plan["others"] == 1 else "ies"))
        if plan["mode"] == "coarse":
            lines.append("TRY          turn it with the mouse  |  snap %s_picture  |  about %s"
                         % (pdb_id, pdb_id))
        else:
            lines.append("TRY          viewas %s surface  |  look rainbow %s  |  about %s"
                         % (name, name, pdb_id))
        cmd.orient(name)
        for line in lines:
            _say(line)
        teach.note("load https://files.rcsb.org/download/%s-assembly%s.cif"
                   % (pdb_id, plan["assembly"]), "the PDB's file of the whole assembly")
        base._last[0] = name
        base._remember(text or pdb_id, name)
        base._remember("%s assembly" % (text or pdb_id), name)
        _say(report.LINE)

    worker.submit(job, done, before)


# ------------------------------------------------------------------ command
_MENU = [
    ("build dna ATGCGTACGTTA", "a DNA double helix of your sequence (form=A for A-DNA)"),
    ("build rna GGCUAGCUAC", "a single RNA strand"),
    ("build helix AKLEELLKKLAE", "a peptide folded as an alpha-helix"),
    ("build strand KVTVKVE", "a peptide as a stretched beta-strand"),
    ("build sheet KVTVKVE TVKVTVE", "strands side by side as a beta-sheet"),
    ("build peptide ACDEFGHIK", "a peptide stretched out, nothing folded"),
    ("build membrane popc 10", "a lipid bilayer, 10 x 10 nm (popc, dopc, dope)"),
    ("build membrane popc cholesterol=30", "the same with 30% cholesterol"),
    ("build micelle sds", "a ball of detergent (sds, ldao)"),
    ("build assembly ferritin", "the complete biological unit of a PDB structure"),
]
_KINDS = {"dna": "dna", "rna": "rna", "helix": "helix", "alpha": "helix", "strand": "strand",
          "sheet": "sheet", "peptide": "peptide", "membrane": "membrane", "bilayer": "membrane",
          "micelle": "micelle", "assembly": "assembly", "complex": "assembly",
          "oligomer": "assembly", "biological": "assembly"}


def build(*args, **kwargs):
    """
DESCRIPTION

    Build a structure: DNA or RNA of any sequence, a peptide as a helix,
    strand or sheet, a lipid membrane, a detergent micelle, or the complete
    biological assembly of a PDB entry.  Every result says whether it is an
    idealised model (geometry rules) or experimental data.

USAGE

    build                                 list what can be built
    build dna <sequence> [form=A|B] [single]
    build rna <sequence>
    build helix | strand | peptide <sequence>
    build sheet <sequence> <sequence> [...]
    build membrane [popc|dopc|dope] [size in nm] [cholesterol=<percent>]
    build micelle [sds|ldao] [number of molecules]
    build assembly <protein name | PDB ID>

EXAMPLES

    build dna ATGCGTACGTTA
    build helix AKLEELLKKLAELLKK
    build sheet KVTVKVE TVKVTVE
    build membrane popc 10 cholesterol=30
    build micelle sds 60
    build assembly ferritin
    """
    words = base._raw(args, kwargs).split()
    label = "build %s" % " ".join(words)
    kind = _KINDS.get(words[0].lower()) if words else None
    rest = words[1:]
    if words and kind is None:
        # "build a double helix of ATGC": find the kind anywhere in the words
        for index, word in enumerate(words):
            if word.lower() in _KINDS:
                kind, rest = _KINDS[word.lower()], words[:index] + words[index + 1:]
                break
    if words and kind is None:
        # a typo of a kind: "build asembly 1J4N", "build membrain popc"
        import difflib
        close = difflib.get_close_matches(words[0].lower(), sorted(_KINDS), n=1, cutoff=0.75)
        if close:
            kind, rest = _KINDS[close[0]], words[1:]
            _say("Reading '%s' as '%s'." % (words[0], close[0]))
    if kind == "assembly":
        _assembly([w for w in rest if w.lower() not in ("unit", "assembly")], label)
        return

    def action():
        base._header(label)
        _say("")
        if kind is None:
            if words:
                _say("I don't know how to build '%s'." % " ".join(words))
            _say("WHAT I CAN BUILD")
            import textwrap
            for example, what in _MENU:
                if len(example) > 30:        # long example: description underneath
                    _say("   " + example)
                    _say(textwrap.fill(what, 78, initial_indent=" " * 34,
                                       subsequent_indent=" " * 34))
                else:
                    _say(textwrap.fill(what, 78, initial_indent="   %-30s " % example,
                                       subsequent_indent=" " * 34))
            _say("Models built from rules are labelled IDEALISED MODEL; assemblies come "
                 "from real PDB data.")
            _say(report.LINE)
            return
        try:
            if kind in ("dna", "rna"):
                _nucleic(kind, rest)
            elif kind in ("helix", "strand", "sheet", "peptide"):
                _peptide(kind, [w for w in rest if w.lower() not in ("of", "a", "an")])
            elif kind == "membrane":
                _membrane(rest)
            elif kind == "micelle":
                _micelle(rest)
        except Exception as error:      # never leave PyMOL in a half-built state
            _say("Building failed: %s" % error)
        _say(report.LINE)

    worker.submit(lambda: None, lambda result, error: action() if not error
                  else base._unexpected(error))


def register():
    for name in ("build", "p_build"):
        cmd.extend(name, build)
    try:
        cmd.auto_arg[0]["build"] = [lambda: cmd.Shortcut(sorted(set(_KINDS.values()))),
                                    "what to build", " "]
    except Exception:
        pass
