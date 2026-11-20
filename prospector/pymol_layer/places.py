"""Key places of a structure, and the `focus` command that zooms to them.

`about` lists the places (active site, ligands, domains, ...); `focus 1` or
`focus active site` then goes there, shows it clearly and says what is on
screen.
"""
from __future__ import annotations

import re

from pymol import cmd

from ..core import http, rcsb, report
from ..core import features as feat
from . import commands as base
from . import teach, worker

_say = base._say
_places = {}        # object name -> [Place]
_current = [None]   # the object whose places were listed last

MEANING = {
    "active": "the residues that carry out the enzyme's chemical reaction",
    "binding": "the residues that hold a partner molecule or metal ion",
    "site": "residues with a special role, such as where the chain gets cut",
    "domain": "a stretch of the chain that folds, and often works, as a unit",
    "disulfide": "pairs of cysteines joined by a sulfur-sulfur bond; they staple the "
                 "fold together",
    "ptm": "residues the cell chemically modifies after the protein is made",
    "variant": "positions where known mutations cause disease",
    "ligand": "a small molecule (or metal ion) bound to the protein",
    "chain": "one of the separate chains that make up this structure",
    "whole": "the complete structure",
}
TITLES = {"active": "Active site", "binding": "Binding sites", "variant": "Disease variants",
          "disulfide": "Disulfide bonds", "ptm": "Modified residues (PTMs)",
          "site": "Other sites"}
# everyday names for common ligand codes, so "focus heme" works
LIGAND_WORDS = {"HEM": ["heme", "haem"], "HEC": ["heme", "haem"], "ZN": ["zinc"],
                "MG": ["magnesium"], "CA": ["calcium"], "FE": ["iron"], "MN": ["manganese"],
                "CU": ["copper"], "NAG": ["sugar", "glycan"], "STI": ["imatinib", "gleevec"],
                "BTN": ["biotin"], "ATP": ["atp"], "ADP": ["adp"], "GTP": ["gtp"],
                "NAD": ["nad"], "FAD": ["fad"], "RET": ["retinal"], "OXY": ["oxygen"]}
AMINO_ACIDS = {
    "ALA": "alanine", "ARG": "arginine", "ASN": "asparagine", "ASP": "aspartate",
    "CYS": "cysteine", "GLN": "glutamine", "GLU": "glutamate", "GLY": "glycine",
    "HIS": "histidine", "ILE": "isoleucine", "LEU": "leucine", "LYS": "lysine",
    "MET": "methionine", "PHE": "phenylalanine", "PRO": "proline", "SER": "serine",
    "THR": "threonine", "TRP": "tryptophan", "TYR": "tyrosine", "VAL": "valine"}
_tipped = [False]
# what people type for each kind of place (used to explain when one is missing)
KIND_WORDS = {
    "active": ("active site", "active", "catalytic site", "active sites"),
    "binding": ("binding site", "binding sites", "binding"),
    "ligand": ("ligand", "ligands", "drug", "cofactor", "inhibitor", "substrate"),
    "variant": ("variants", "variant", "disease variants", "mutations"),
    "domain": ("domain", "domains"),
    "disulfide": ("disulfides", "disulfide", "disulfide bonds"),
    "ptm": ("ptm", "ptms", "modifications", "modified residues"),
}
# words that describe a place rather than name a protein
PLACE_WORDS = set("active site sites catalytic ligand ligands binding domain domains "
                  "variant variants disease mutations disulfide disulfides bonds ptm ptms "
                  "modifications modified residues residue chain all whole everything "
                  "structure back out the of in to".split())


class Place(object):
    def __init__(self, kind, title, selection, colour, style, words, detail=""):
        self.kind = kind            # key of MEANING
        self.title = title          # "Active site", "Ligand HEM", "Domain: Kringle 1"
        self.selection = selection  # PyMOL selection expression
        self.colour = colour
        self.style = style          # sticks | spheres | cartoon | ligand | chain | whole
        self.words = words          # what the user may type to get here
        self.detail = detail        # short text for the list (residues, name ...)
        self.residues = 0


def _residue_labels(selection, limit=6):
    """(count, 'H/HIS57, H/ASP102, ...') for the residues really in the object."""
    found = []
    cmd.iterate("(%s) and name CA" % selection, "found.append((chain, resn, resi))",
                space={"found": found})
    text = ", ".join("%s%s%s" % ("%s/" % c if c else "", n, i) for c, n, i in found[:limit])
    if len(found) > limit:
        text += " +%d more" % (len(found) - limit)
    return len(found), text


def build(obj, annotation=None, names=None):
    """Work out the key places of a loaded object.  ``annotation`` (UniProt
    features) and ``names`` (ligand full names) make the list richer but are
    optional."""
    names = names or {}
    places = []
    plan = annotation.selections() if annotation else []

    def add_feature(key, word_list):
        for suffix, group, title, placed, colour, style in plan:
            if group != key:
                continue
            selection = feat.selection_text(obj, placed)
            count, labels = _residue_labels(selection)
            if not count:
                continue
            words = list(word_list)
            chain_text = "chain%s %s" % ("" if len(placed) == 1 else "s",
                                         ", ".join(sorted(placed)))
            if key == "domain":
                short = title.split(": ", 1)[-1]
                words = ["domain " + short.lower(), short.lower()] + words
                detail = "%d residues, %s" % (count, chain_text)
            else:
                title = TITLES.get(key, title)
                detail = _residue_labels(selection, 3)[1] if count <= 6 \
                    else "%d residues, %s" % (count, chain_text)
            place = Place(key, title, selection, colour or "orange", style, words, detail)
            place.residues = count
            places.append(place)

    add_feature("active", ["active site", "active", "catalytic site"])

    ligands = base._ligands_in(obj)
    seen = set()
    for chain, resn, resi, atoms in ligands:
        if resn in seen or resn in rcsb.ADDITIVES:
            continue
        seen.add(resn)
        full = names.get(resn, ("", ""))[0]
        copies = len([1 for item in ligands if item[1] == resn])
        detail = full[:38] if full else "%d atoms" % atoms
        if resn in LIGAND_WORDS:
            detail = LIGAND_WORDS[resn][0]
        if copies > 1:
            detail += " (%d copies)" % copies
        words = [resn.lower(), "ligand " + resn.lower()] + LIGAND_WORDS.get(resn, [])
        words += [w for w in re.findall(r"[a-z]{5,}", full.lower()) if w != "containing"]
        if not any(p.kind == "ligand" for p in places):
            words.append("ligand")
        places.append(Place("ligand", "Ligand %s" % resn,
                            base._residue_selection(obj, chain, resn, resi),
                            "yellow", "ligand", words, detail))

    add_feature("binding", ["binding site", "binding", "binding sites"])
    add_feature("domain", ["domain", "domains"])
    add_feature("variant", ["disease variants", "variants", "variant", "mutations"])
    add_feature("disulfide", ["disulfides", "disulfide", "disulfide bonds", "ss"])
    add_feature("ptm", ["ptm", "ptms", "modifications", "modified residues"])

    chain_names = {}
    if annotation:
        for protein in annotation.proteins:
            for chain_map in protein.chains:
                chain_names[chain_map.chain] = protein.protein.name
    chains = [c for c in cmd.get_chains("%s and polymer" % obj) if c]
    if 1 < len(chains) <= 8:
        for chain in chains:
            detail = chain_names.get(chain, "")
            if not detail and cmd.count_atoms("%s and chain %s and polymer.nucleic"
                                              % (obj, chain)):
                detail = "DNA or RNA"
            places.append(Place("chain", "Chain %s" % chain,
                                "%s and chain %s" % (obj, chain), None, "chain",
                                ["chain " + chain.lower()], detail))
    places.append(Place("whole", "Whole structure", obj, None, "whole",
                        ["all", "whole", "everything", "whole structure", "back", "out"],
                        "%d chain%s" % (len(chains) or 1, "" if len(chains) <= 1 else "s")))
    _places[obj] = places
    _current[0] = obj
    return places


def listing(obj):
    """The KEY PLACES block printed by `about` and by `focus` without arguments."""
    places = _places.get(obj) or []
    out = [" KEY PLACES in %s  (type the command on the right to zoom there)" % obj,
           " " + "-" * 76]
    for number, place in enumerate(places, 1):
        out.append("   %2d  %-26s %-34s focus %d" % (number, place.title[:26],
                                                    place.detail[:34], number))
    out.append("   You can also use the name, e.g.   focus %s"
               % (places[0].words[0] if places else "all"))
    out.append("   or go to one residue, e.g.        focus %s"
               % _example_residue(obj))
    return "\n".join(out)


def _example_residue(obj):
    found = []
    cmd.iterate("first (%s and polymer and name CA)" % obj, "found.append((chain, resi))",
                space={"found": found})
    if found and found[0][0]:
        return "%s/%s" % found[0]
    return found[0][1] if found else "A/10"


def _match(places, text):
    """The place the user means: a number, a name, or None."""
    wanted = " ".join(text.lower().split())
    if len(wanted) == 1 and wanted.isalpha():      # "focus B" = chain B, if there is one
        for place in places:
            if place.title.lower().startswith("chain %s" % wanted):
                return place
    if wanted.isdigit():
        index = int(wanted) - 1
        return places[index] if 0 <= index < len(places) else None
    for place in places:
        if wanted in place.words or wanted == place.title.lower():
            return place
    for place in places:      # looser: "kringle" finds "Domain: Kringle 1"
        if any(wanted in word for word in place.words) or wanted in place.title.lower():
            return place
    import difflib            # typos: "activ site", "lignd"
    every = dict((word, place) for place in reversed(places) for word in place.words)
    close = difflib.get_close_matches(wanted, list(every), n=1, cutoff=0.75)
    return every[close[0]] if close else None


_RESIDUE = re.compile(r"^(?:(?:chain\s+)?([A-Za-z0-9]{1,4})\s*/\s*)?"
                      r"(?:resi(?:due)?\s+)?(?:[A-Za-z]{3})?(-?\d+[A-Za-z]?)$")


def _residue_place(obj, text):
    """'H/57', '57', 'A/HIS87' or 'residue 57' -> a one-residue Place, or None."""
    match = _RESIDUE.match(text.strip())
    if not match:
        return None
    chain, number = match.group(1), match.group(2)
    label = ("\\" + number) if number.startswith("-") else number
    selection = "%s and polymer and resi %s" % (obj, label)
    if chain:
        selection += " and chain %s" % chain
    count, labels = _residue_labels(selection)
    if not count:
        return None
    names = []
    cmd.iterate("(%s) and name CA" % selection, "names.append(resn)", space={"names": names})
    what = AMINO_ACIDS.get(names[0].upper(), "") if names else ""
    place = Place("residue", "Residue %s" % labels, selection, "orange", "sticks", [],
                  "%s = %s, residue number %s%s" % (names[0], what, number,
                                                   " of chain %s" % chain if chain else "")
                  if what else "")
    place.residues = count
    return place


def _apply(obj, place):
    """Go to a place: make it visible, zoom, and describe what is on screen."""
    base._snapshot("focus %s" % place.title.lower())
    seen, hidden, uncovered = _apply_quiet(obj, place)
    _describe(obj, place, seen, hidden, uncovered)


def _apply_quiet(obj, place):
    """Draw and zoom to a place without printing.  Returns what is on screen."""
    cmd.enable(obj)
    cmd.label(obj, "")                       # remove labels of the previous focus
    sel = place.selection
    seen = []
    uncovered = False
    # one structure at a time: others would be drawn straight through this one
    hidden = [name for name in base._molecules(enabled_only=True) if name != obj]
    for name in hidden:
        cmd.disable(name)
        teach.note("disable %s" % name, "hide the other structure (it stays loaded)")
    close_up = place.style in ("ligand", "sticks", "spheres") and place.residues <= 40
    if place.style not in ("whole", "chain"):
        base._clean(obj)          # one highlight at a time
    # a see-through ribbon keeps the context without hiding the sticks
    cmd.set("cartoon_transparency", 0.65 if close_up else 0.0, obj)
    if close_up:
        teach.note("set cartoon_transparency, 0.65, %s" % obj, "see-through ribbon")
        cmd.hide("nonbonded", "%s and solvent" % obj)      # water crosses are clutter here
        teach.note("hide nonbonded, %s and solvent" % obj, "hide the water molecules")
    if place.style not in ("whole", "chain"):
        # a surface or putty tube from an earlier "look" would hide the close-up
        if cmd.count_atoms("%s and rep surface" % obj):
            cmd.hide("surface", obj)
            teach.note("hide surface, %s" % obj, "the surface would cover the close-up")
            uncovered = True
        cmd.show("cartoon", "%s and polymer" % obj)
        cmd.cartoon("automatic", obj)
    if place.style == "whole":
        cmd.orient(obj)
        teach.note("orient %s" % obj, "turn and zoom so the whole structure fills the view")
        seen.append("the whole structure, turned to fill the window")
    elif place.style == "chain":
        cmd.show("cartoon", "%s and polymer" % sel)
        cmd.orient(sel)
        teach.note("orient %s" % sel, "fill the view with this chain")
        seen.append("this chain in the middle of the window; the other chains are still "
                    "there around it")
    elif place.style == "cartoon":           # a domain
        cmd.show("cartoon", "%s and polymer" % obj)
        cmd.color("grey80", "%s and polymer and elem C" % obj)
        cmd.color(place.colour, "(%s) and elem C" % sel)
        cmd.orient(sel)
        cmd.zoom(sel, 4)
        teach.note("color grey80, %s and polymer and elem C" % obj, "everything else grey")
        teach.note("color %s, <domain residues> and elem C" % place.colour)
        teach.note("zoom <domain residues>, 4", "4 A of space around it")
        seen.append("the domain as a %s ribbon; the rest of the protein is grey"
                    % place.colour)
    elif place.style == "ligand":
        pocket = "byres (%s and polymer within 4 of (%s))" % (obj, sel)
        cmd.show("sticks", "(%s) or (%s)" % (sel, pocket))
        cmd.color("yellow", "(%s) and elem C" % sel)
        cmd.util.cnc("(%s) or (%s)" % (sel, pocket))
        cmd.orient(sel)
        cmd.zoom(sel, 8)
        teach.note("show sticks, (%s) or (byres polymer within 4 of it)" % sel)
        teach.note("color yellow, (%s) and elem C" % sel, "ligand carbons yellow")
        teach.note("zoom %s, 8" % sel, "8 A of space around the ligand")
        count = cmd.count_atoms("(%s) and name CA" % pocket)
        seen.append("the ligand as yellow sticks in the centre")
        seen.append("the %d protein residues touching it as thinner sticks "
                    "(blue = nitrogen, red = oxygen)" % count)
        seen.append("the rest of the protein as a faint, see-through ribbon")
    else:                                    # residues: sites, variants, one residue
        count = place.residues
        if place.style == "spheres":
            cmd.show("spheres", "(%s) and name CA" % sel)
            cmd.set("sphere_scale", 0.7, "(%s) and name CA" % sel)
            cmd.color(place.colour, "(%s) and name CA" % sel)
            teach.note("show spheres, <these residues> and name CA", "one ball per residue")
            seen.append("one %s ball at each of the %d positions" % (place.colour, count))
        else:
            cmd.show("sticks", "(%s) and (sidechain or name CA)" % sel)
            cmd.color(place.colour, "(%s) and elem C" % sel)
            cmd.util.cnc("(%s)" % sel)
            teach.note("show sticks, <these residues> and (sidechain or name CA)")
            teach.note("color %s, <these residues> and elem C" % place.colour)
            seen.append("%d residue%s as %s sticks (blue = nitrogen, red = oxygen, "
                        "yellow = sulfur)" % (count, "" if count == 1 else "s", place.colour))
        if count <= 12:
            cmd.label("(%s) and name CA" % sel, '"%s%s" % (resn, resi)')
            teach.note('label <these residues> and name CA, "%s%s" % (resn, resi)',
                       "write each residue's name next to it")
            seen.append("a name tag next to each: amino acid + its number in the chain")
        cmd.orient(sel)
        cmd.zoom(sel, 8 if count <= 12 else 4)
        teach.note("zoom <these residues>, 8", "8 A of space around them")
        seen.append("the surrounding protein as a faint, see-through ribbon"
                    if close_up else "the ribbon of the whole protein")
    cmd.deselect()
    return seen, hidden, uncovered


def _describe(obj, place, seen, hidden, uncovered):
    _say("")
    _say("NOW SHOWING  %s of %s" % (place.title[:1].lower() + place.title[1:], obj))
    if place.detail and place.style not in ("whole",):
        _say("             %s" % place.detail)
    _say("ON SCREEN    " + seen[0])
    for line in seen[1:]:
        _say("             " + line)
    if hidden:
        _say("             (%s %s hidden so nothing overlaps;  enable all  shows %s again)"
             % (", ".join(hidden[:3]), "is" if len(hidden) == 1 else "are",
                "it" if len(hidden) == 1 else "them"))
    if uncovered:
        _say("             (the surface was hidden so you can see inside;  look ligand  "
             "or  look hydrophobicity  brings it back)")
    if place.kind in MEANING:
        _say("WHAT IT IS   " + MEANING[place.kind])
    if place.style != "whole":
        _say("GO BACK      focus all  (whole structure)   |   focus  (list of places)")
    else:
        _say("GO CLOSER    focus  (list of places)")
    if not _tipped[0]:
        _tipped[0] = True
        _say("MOUSE        left-drag = rotate  |  scroll = zoom  |  middle-drag = move")


def focus(*args, **kwargs):
    """
DESCRIPTION

    Zoom to an important place of a structure and say what is on screen.
    "about" lists the places; "focus" alone lists them again.

USAGE

    focus                          list the key places
    focus <number | name> [object]
    focus <chain/residue>          e.g.  focus H/57
    focus <place> of <protein>     loads the protein first if needed

EXAMPLES

    focus 1
    focus active site
    focus ligand
    focus HEM
    focus chain B
    focus A/87
    focus active site of thrombin
    focus all                      back to the whole structure
    """
    words = base._raw(args, kwargs).replace(",", " ").split()
    label = "focus %s" % " ".join(words)
    state = {}

    def before():
        typed, place_words = "", list(words)
        # "focus active site of haemoglobin": the words after of/in name the protein
        lowered = [w.lower() for w in words]
        for marker in ("of", "in", "on", "for"):
            if marker in lowered[1:-1]:
                cut = len(lowered) - 1 - lowered[::-1].index(marker)
                named = " ".join(words[cut + 1:])
                if base._known_object(named):
                    typed, place_words = named, words[:cut]
                else:
                    state.update(obj=None, query=named, need=True, identity=None,
                                 text=" ".join(words[:cut]))
                    return True
                break
        if typed:
            obj = base._known_object(typed)
            state.update(obj=obj, text=" ".join(place_words),
                         identity=base._identify(obj), need=obj not in _places)
            return True
        for size in range(min(3, len(words)), 0, -1):      # "... p53" at the end?
            if base._known_object(" ".join(words[-size:])) and len(words) > size:
                typed, place_words = " ".join(words[-size:]), words[:-size]
                break
        if not typed:
            for size in range(min(3, len(words) - 1), 0, -1):   # "p53 ..." at the start?
                if base._known_object(" ".join(words[:size])):
                    typed, place_words = " ".join(words[:size]), words[size:]
                    break
        if not typed and words and base._known_object(" ".join(words)):
            typed, place_words = " ".join(words), []       # "focus p53" = list its places
        everyday = set(word for words_ in LIGAND_WORDS.values() for word in words_)
        name_words = [w for w in place_words if w.lower() not in PLACE_WORDS
                      and not w.isdigit() and not _RESIDUE.match(w)
                      and w.upper() not in LIGAND_WORDS and w.lower() not in everyday]
        if not typed and name_words and not base._molecules():
            # nothing is loaded and the request names a protein: fetch it first
            state.update(obj=None, query=" ".join(name_words), need=True, identity=None,
                         text=" ".join(w for w in place_words if w not in name_words))
            return True
        if not typed and _current[0] in base._molecules():
            obj = _current[0]
        else:
            obj, message = base._pick_object(typed, "focus 1")
            if obj is None:
                base._header(label)
                _say(message)
                _say(report.LINE)
                return False
        state.update(obj=obj, text=" ".join(place_words), identity=base._identify(obj),
                     need=obj not in _places)
        return True

    def job():   # annotations are fetched once per object
        if state.get("query"):
            state["fetched"] = base._fetch_query(state["query"])
            if state["fetched"]["error"]:
                return None, {}
            state["identity"] = state["fetched"]["identity"]
        if not state["need"]:
            return None, {}
        return fetch_annotation(state["identity"])

    def done(answer, error):
        if error:
            return base._unexpected(error)
        base._header(label)
        if state.get("fetched"):
            if state["fetched"]["error"]:
                _say(state["fetched"]["error"])
                _say(report.LINE)
                return
            state["obj"] = base._place_fetched(state["fetched"])
        obj = state["obj"]
        if state["need"]:
            build(obj, answer[0], answer[1])
        _current[0] = obj
        places = _places[obj]
        if not state["text"]:
            _say(listing(obj))
            _say(report.LINE)
            return
        place = _match(places, state["text"]) or _residue_place(obj, state["text"])
        wanted = " ".join(state["text"].lower().split())
        asked_kind = None
        for kind, kind_words in KIND_WORDS.items():
            if wanted in kind_words:
                asked_kind = kind
        if asked_kind and (place is None or place.kind != asked_kind):
            # the structure has none of what was asked for: say so plainly, and
            # for "active site" / "binding site" show the nearest working part
            names = {"active": "active site", "binding": "binding site",
                     "variant": "disease variants", "domain": "domains",
                     "disulfide": "disulfide bonds", "ptm": "modified residues",
                     "ligand": "bound molecule"}
            order = {"active": ("binding", "ligand"), "binding": ("active", "ligand"),
                     "ligand": ("binding", "active")}.get(asked_kind, ())
            stand_in = [p for kind in order for p in places if p.kind == kind]
            _say("")
            if asked_kind == "active":
                why = ("an active site is where an ENZYME does its chemistry, and UniProt "
                       "does not annotate one for this protein")
            elif asked_kind == "ligand":
                why = "nothing but protein (and water) is in this structure"
            else:
                why = ("UniProt annotates none in the part of the protein that this "
                       "structure contains")
            _say("NOTE         %s has no %s: %s." % (obj, names[asked_kind], why))
            if stand_in:
                place = stand_in[0]
                _say("             Its working part is shown instead: %s."
                     % place.title.lower())
            else:
                _say("")
                _say(listing(obj))
                _say(report.LINE)
                return
        if place is None:
            _say("I don't know a place called '%s' in %s. These are the ones I know:"
                 % (state["text"], obj))
            _say("")
            _say(listing(obj))
            _say(report.LINE)
            return
        _apply(obj, place)
        _say(report.LINE)

    worker.submit(job, done, before)


def _onto_copies(annotation, copies):
    """A whole assembly: chains C, D ... are copies of A, B ... - give them the same
    UniProt mapping and features (on a copy, so the cached annotation is untouched)."""
    import copy
    annotation = copy.deepcopy(annotation)
    by_original = {}
    for new, old in copies.items():
        by_original.setdefault(old, []).append(new)
    for protein in annotation.proteins:
        for chain_map in list(protein.chains):
            for new in by_original.get(chain_map.chain, []):
                twin = copy.copy(chain_map)
                twin.chain = new
                protein.chains.append(twin)
        for feature in protein.features:
            for old, labels in list(feature.placed.items()):
                for new in by_original.get(old, []):
                    feature.placed[new] = list(labels)
    return annotation


def fetch_annotation(identity):
    """(annotation or None, ligand names) for an object identity; never raises."""
    if identity is None or identity[0] in ("mol", "built"):
        return None, {}
    annotation, names = None, {}
    try:
        if identity[0] == "pdb":
            names = rcsb.ligand_names(rcsb.entry_detail(identity[1]))
            annotation = feat.annotate(identity[1])
            copies = base._assembled.get(identity[1].upper()) or {}
            if annotation is not None and copies:
                annotation = _onto_copies(annotation, copies)
        else:
            annotation = feat.annotate_model(identity[1])
    except (feat.FeatureError, http.NetworkError):
        pass
    return annotation, names


def register():
    cmd.extend("focus", focus)
    cmd.extend("p_focus", focus)
