"""Tier 4 - the learning layer: tour (+ next / prev), quiz, surprise."""
from __future__ import annotations

import json
import os
import random

from pymol import cmd

from ..core import card, http, rank, rcsb, report
from ..core import features as feat
from . import commands as base
from . import teach, worker

_say = base._say
_tour = {"object": None, "scenes": [], "index": -1}     # scenes: [(name, title, narration)]
_quiz = {"entry": None, "object": None, "hints": 0, "last": None}
_surprised = []                                           # PDB ids already shown this session


def _data(name):
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "core", "data", name)
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _do(command, comment=""):
    """Run one PyMOL command line and record it for teach mode."""
    cmd.do(command, echo=0)
    teach.note(command, comment)


# --------------------------------------------------------------------- tour
def _secondary_structure(obj):
    counts = {"H": 0, "S": 0, "": 0}
    cmd.iterate("%s and polymer and name CA" % obj,
                "counts[ss if ss in ('H', 'S') else ''] += 1", space={"counts": counts})
    if counts["H"] + counts["S"] == 0 and counts[""] > 0:
        cmd.dss(obj)          # the file had no helix/sheet records: let PyMOL assign them
        counts = {"H": 0, "S": 0, "": 0}
        cmd.iterate("%s and polymer and name CA" % obj,
                    "counts[ss if ss in ('H', 'S') else ''] += 1", space={"counts": counts})
    total = float(sum(counts.values()) or 1)
    return int(sum(counts.values())), 100 * counts["H"] / total, 100 * counts["S"] / total


def _fold_words(helix, strand):
    if helix >= 40 and strand < 10:
        return "mostly alpha-helices"
    if strand >= 30 and helix < 15:
        return "mostly beta-strands"
    if helix < 10 and strand < 10:
        return "little regular secondary structure"
    return "a mix of alpha-helices and beta-strands"


_COLOUR_WORDS = {"lightblue": "light blue", "palegreen": "pale green", "wheat": "wheat",
                 "lightpink": "light pink", "palecyan": "pale cyan",
                 "lightorange": "light orange", "bluewhite": "blue-white", "limon": "lime"}


def _base_view(obj):
    cmd.disable("tour_polar")
    cmd.set("cartoon_transparency", 0, obj)     # a `focus` close-up may have set it
    cmd.label(obj, "")
    _do("hide everything, %s" % obj, "start each scene from a blank slate")
    _do("show cartoon, %s and polymer" % obj)
    _do("color grey80, %s and polymer" % obj, "neutral background colour")


def _store(number, title, message, lines):
    """Save the current view as a tour stop.  ``message`` is the one-liner PyMOL
    shows in the viewer; ``lines`` is the explanation: [(label, text), ...]."""
    name = "tour_%d" % number
    cmd.scene(name, "store", message=message)
    teach.note("scene %s, store" % name, "remember view + colours as scene %d" % number)
    _tour["scenes"].append((name, title, lines))


def _say_stop(index):
    """Print one stop: what is on screen, how it is drawn, what it means."""
    scenes = _tour["scenes"]
    name, title, lines = scenes[index]
    _say("[%d/%d] %s" % (index + 1, len(scenes), title.upper()))
    for label, text in lines:
        _say("   %-12s %s" % (label, text))


def _names_in(selection, limit=6):
    """'HIS57, ASP102, SER195' for the residues of a selection."""
    found = []
    cmd.iterate("(%s) and name CA" % selection, "found.append(resn + resi)",
                space={"found": found})
    seen = []
    for item in found:
        if item not in seen:
            seen.append(item)
    text = ", ".join(seen[:limit])
    return text + (" and %d more" % (len(seen) - limit) if len(seen) > limit else "")


def tour(*args, **kwargs):
    """
DESCRIPTION

    A guided tour of a loaded structure, built from PyMOL scenes:
    overall fold -> domains -> functional sites -> ligand pocket -> disease
    variants (steps without data are skipped).  One line of narration each.
    Step through with  next  and  prev.

USAGE

    tour [object | protein name]
    tour <number>                  jump to that stop of the running tour

EXAMPLES

    tour 1PPB           (a loaded object)
    tour haemoglobin    (finds and loads the best match first)
    next
    """
    typed = base._raw(args, kwargs).strip()
    label = "tour %s" % typed
    state = {}
    if typed.isdigit() and not base._known_object(typed):
        # "tour 3" = jump to stop 3 of the tour that is running
        return _step(int(typed) - 1, label, absolute=True)

    def before():
        obj, message = base._pick_object(typed, "tour")
        if obj is None and typed:
            state["query"] = typed          # a protein name: find and load it first
            return True
        if obj is None:
            base._header(label)
            _say(message)
            _say(report.LINE)
            return False
        state["obj"] = obj
        state["identity"] = base._identify(obj)
        return True

    def job():   # annotations make the tour richer but are never essential
        info = {"annotation": None, "names": {}, "title": "", "note": ""}
        if state.get("query"):
            state["fetched"] = base._fetch_query(state["query"])
            if state["fetched"]["error"]:
                return info
            state["identity"] = state["fetched"]["identity"]
        identity = state["identity"]
        if identity is None:
            info["note"] = ("(The object is not named after a PDB ID, so the tour has no "
                            "UniProt annotations to show.)")
            return info
        try:
            if identity[0] == "pdb":
                detail = rcsb.entry_detail(identity[1])
                info["names"] = rcsb.ligand_names(detail)
                info["title"] = ((detail or {}).get("struct") or {}).get("title") or ""
                info["annotation"] = feat.annotate(identity[1])
            else:
                info["annotation"] = feat.annotate_model(identity[1])
        except feat.FeatureError:
            pass
        except http.NetworkError:
            info["note"] = "(Offline: touring without UniProt annotations.)"
        return info

    def done(info, error):
        if error:
            return base._unexpected(error)
        base._header(label)
        if state.get("fetched"):
            if state["fetched"]["error"]:
                _say(state["fetched"]["error"])
                _say(report.LINE)
                return
            state["obj"] = base._place_fetched(state["fetched"])
        if base._needs_protein(state["obj"], "tour"):
            return
        if info["note"]:
            _say(info["note"])
        base._snapshot(label)
        _build_tour(state["obj"], info)
        _say(report.LINE)

    worker.submit(job, done, before)


def _build_tour(obj, info):
    for name in cmd.get_scene_list() or []:
        if name.startswith("tour_"):
            cmd.scene(name, "clear")
    _tour.update(object=obj, scenes=[], index=-1)
    annotation = info["annotation"]
    plan = annotation.selections() if annotation else []
    protein_name = ""
    if annotation and annotation.proteins:
        protein_name = annotation.proteins[0].protein.name

    # the tour shows one structure: everything else is switched off (not deleted)
    others = [n for n in cmd.get_names("objects") if n != obj]
    if others:
        _do("disable all", "hide the other objects (they stay loaded)")
        _do("enable %s" % obj)
    cmd.delete("tour_polar")

    # 1. overall fold
    residues, helix, strand = _secondary_structure(obj)
    chains = [c for c in cmd.get_chains(obj)
              if cmd.count_atoms("%s and polymer and chain %s" % (obj, c) if c
                                 else "%s and polymer" % obj)]
    cmd.set("cartoon_transparency", 0, obj)
    cmd.label(obj, "")
    _do("hide everything, %s" % obj, "start from a blank slate")
    _do("show cartoon, %s and polymer" % obj, "ribbons show the fold")
    if len(chains) > 1:
        _do("util.cbc %s" % obj, "one colour per chain")
    else:
        _do("util.chainbow %s" % obj, "rainbow: blue N-terminus to red C-terminus")
    _do("orient %s" % obj)
    what = protein_name or " ".join(info["title"].split())[:50] or obj
    if annotation and len(annotation.proteins) > 1:
        distinct = []
        for item in annotation.proteins:
            if item.protein.name not in distinct:
                distinct.append(item.protein.name)
        what = " + ".join(distinct[:2]) + (" + more" if len(distinct) > 2 else "")
    if len(chains) > 1:
        colours = ("each of the %d chains has its own colour, so you can tell them apart"
                   % len(chains))
    else:
        colours = ("a rainbow along the chain: blue = its start (N-terminus), then green "
                   "and yellow, red = its end (C-terminus)")
    _store(1, "Overall fold", "Overall fold of %s" % what, [
        ("You see", "the whole of %s as a ribbon that follows the protein backbone" % what),
        ("Shapes", "spirals = alpha-helices, flat arrows = beta-strands, thin tubes = "
                   "loops"),
        ("Colours", colours),
        ("Facts", "%d chain%s, %d residues - %s (%.0f%% helix, %.0f%% strand)"
                  % (len(chains), "" if len(chains) == 1 else "s", residues,
                     _fold_words(helix, strand), helix, strand)),
        ("Try", "drag with the left mouse button to turn it; scroll to zoom"),
    ])

    # 2. domains
    domains = [item for item in plan if item[1] == "domain"]
    if domains:
        _base_view(obj)
        described = []
        for suffix, key, title, placed, colour, style in domains:
            selection = feat.selection_text(obj, placed)
            cmd.color(colour, selection)
            short = selection if len(selection) < 60 else selection[:56] + "...)"
            teach.note("color %s, %s" % (colour, short), title)
            chains_text = "chain%s %s" % ("" if len(placed) == 1 else "s",
                                          ", ".join(sorted(placed)))
            described.append("%s, %s (%s)" % (title.split(": ", 1)[-1] or title,
                                              chains_text, _COLOUR_WORDS.get(colour, colour)))
        _do("orient %s" % obj)
        _store(len(_tour["scenes"]) + 1, "Domains",
               "Domains: %s" % ", ".join(described[:4]), [
                   ("You see", "the same ribbon, now grey, with each domain painted in its "
                               "own colour"),
                   ("Colours", "; ".join(d.replace(" (", " = ").rstrip(")")
                                         for d in described[:6])),
                   ("It means", "a domain is a stretch of the chain that folds, and often "
                                "works, as a unit"),
                   ("Look for", "where one colour ends and the next begins - that is a "
                                "domain boundary; grey parts belong to no annotated domain"),
               ])

    # 3. functional sites
    sites = [item for item in plan if item[1] in ("active", "binding")]
    site_parts, site_text = [], []
    for suffix, key, title, placed, colour, style in sites:
        selection = feat.selection_text(obj, placed)
        if cmd.count_atoms(selection) == 0:
            continue
        site_parts.append((selection, colour, title))
        first_chain = sorted(placed)[0]
        one_copy = "%s and chain %s and (%s)" % (obj, first_chain, selection)
        site_text.append("%s sticks = %s (%s)" % (colour, title.lower(),
                                                  _names_in(one_copy)))
    if site_parts:
        _base_view(obj)
        for selection, colour, title in site_parts:
            cmd.show("sticks", "(%s) and (sidechain or name CA)" % selection)
            cmd.color(colour, "(%s) and elem C" % selection)
            short = selection if len(selection) < 60 else selection[:56] + "...)"
            teach.note("show sticks, (%s) and (sidechain or name CA)" % short, title)
            teach.note("color %s, (%s) and elem C" % (colour, short))
            cmd.util.cnc("(%s)" % selection)
        teach.note("util.cnc <site residues>", "N blue, O red; carbons keep their colour")
        first = site_parts[0][0]
        cmd.orient(first)
        cmd.zoom(first, 10)
        teach.note("orient %s" % (first if len(first) < 60 else first[:56] + "...)"),
                   "turn so the site faces you")
        teach.note("zoom %s, 10" % (first if len(first) < 60 else first[:56] + "...)"),
                   "10 A of margin around the site")
        _store(len(_tour["scenes"]) + 1, "Functional sites",
               "Functional sites: %s" % "; ".join(site_text)[:90], [
                   ("You see", "a close-up: a few residues drawn as sticks, atom by atom, "
                               "inside the grey ribbon"),
                   ("Colours", "; ".join(site_text)),
                   ("Atoms", "on the sticks, blue = nitrogen, red = oxygen, yellow = "
                             "sulfur"),
                   ("It means", "active-site residues do the chemistry; binding-site "
                                "residues hold a partner molecule or metal ion"),
                   ("Look for", "how residues far apart in the chain end up side by side "
                                "in space"),
               ])

    # 4. ligand pocket
    ligands = base._ligands_in(obj)
    real = [item for item in ligands if item[1] not in rcsb.ADDITIVES]
    if real:
        chain, resn, resi, atoms = real[0]
        ligand = base._residue_selection(obj, chain, resn, resi)
        pocket = "byres (%s and polymer within 4 of (%s))" % (obj, ligand)
        _base_view(obj)
        _do("show sticks, (%s) or (%s)" % (ligand, pocket), "ligand and the residues around it")
        _do("color yellow, (%s) and elem C" % ligand, "ligand carbons yellow")
        cmd.util.cnc("(%s) or (%s)" % (ligand, pocket))
        teach.note("util.cnc <ligand and pocket>", "N blue, O red; carbons keep their colour")
        if cmd.count_atoms(pocket):
            _do("distance tour_polar, (%s), (%s), mode=2" % (ligand, pocket),
                "dashed lines = polar contacts")
            cmd.hide("labels", "tour_polar")
        _do("orient (%s)" % ligand, "turn so the ligand lies flat")
        _do("zoom (%s), 8" % ligand)
        count = cmd.count_atoms("(%s) and name CA" % pocket)
        full = info["names"].get(resn, ("", ""))[0]
        others = sorted(set(item[1] for item in real) - set([resn]))
        lines = [
            ("You see", "the bound molecule %s in the centre, as yellow sticks" % resn),
            ("Around it", "the %d protein residues that touch it, as thinner grey sticks "
                          "(blue = nitrogen, red = oxygen)" % count),
            ("Dashes", "yellow dashed lines = likely hydrogen bonds between ligand and "
                       "protein"),
        ]
        if full:
            lines.insert(1, ("It is", full[:60]))
        lines.append(("Look for", "how snugly the pocket wraps the ligand, and which "
                                  "residues the dashes reach"))
        if others:
            lines.append(("Also here", "other bound molecules: %s  (see them with  "
                                       "pocket %s %s)" % (", ".join(others[:4]), obj,
                                                           others[0])))
        _store(len(_tour["scenes"]) + 1, "Ligand",
               "Ligand %s in its pocket of %d residues" % (resn, count), lines)

    # 5. disease variants
    variants = [item for item in plan if item[1] == "variant"]
    if variants:
        selection = feat.selection_text(obj, variants[0][3])
        count = cmd.count_atoms("(%s) and name CA" % selection)
        if count:
            _base_view(obj)
            cmd.show("spheres", "(%s) and name CA" % selection)
            cmd.set("sphere_scale", 0.7, "(%s) and name CA" % selection)
            cmd.color("purple", "(%s) and name CA" % selection)
            teach.note("show spheres, <variant residues> and name CA",
                       "one sphere per mutated position")
            _do("orient %s" % obj)
            _store(len(_tour["scenes"]) + 1, "Disease variants",
                   "%d disease-variant positions (purple balls)" % count, [
                       ("You see", "the whole grey ribbon again, with a purple ball at "
                                   "%d position%s" % (count, "" if count == 1 else "s")),
                       ("It means", "at each ball, a known mutation of that residue "
                                    "causes a disease"),
                       ("Look for", "where the balls cluster - often near the active "
                                    "site, a binding site or the protein's core"),
                       ("Details", "features %s types=variant  lists every one" % obj),
                   ])

    cmd.deselect()
    _tour["index"] = 0
    cmd.scene(_tour["scenes"][0][0], "recall", animate=0)
    _say("")
    _say("Tour of %s: %d stops." % (obj, len(_tour["scenes"])))
    for number, (name, title, lines) in enumerate(_tour["scenes"], 1):
        _say("  %d. %s" % (number, title))
    _say("At each stop the picture changes and the text says what you are looking at.")
    _say("")
    _say_stop(0)
    _say("")
    _say("Type  next  to continue,  prev  to go back, or  tour 3  to jump to stop 3.")
    if others:
        _say("(Other objects are hidden during the tour;  enable all  shows them again.)")
    if len(_tour["scenes"]) == 1:
        _say("(Only one stop: this structure has no annotated domains, sites, ligands or "
             "variants to show.)")


def _step(delta, label, absolute=False):
    def action():
        base._header(label)
        scenes = _tour["scenes"]
        if not scenes or _tour["object"] not in cmd.get_names("objects"):
            _say("No tour is running. Start one with  tour <object>  (after loading a "
                 "structure, e.g.  grab 1PPB).")
            _say(report.LINE)
            return
        index = delta if absolute else _tour["index"] + delta
        if absolute and not 0 <= index < len(scenes):
            _say("This tour has stops 1 to %d. Type e.g.  tour %d" % (len(scenes),
                                                                     len(scenes)))
            _say(report.LINE)
            return
        if index >= len(scenes):
            _say("That was the last stop of the tour of %s.  prev  goes back;  "
                 "tour %s  starts again." % (_tour["object"], _tour["object"]))
            _say(report.LINE)
            return
        if index < 0:
            _say("This is the first stop. Type  next  to continue.")
            _say(report.LINE)
            return
        _tour["index"] = index
        name = scenes[index][0]
        cmd.scene(name, "recall")
        teach.note("scene %s, recall" % name, "jump to a stored scene (PyMOL animates it)")
        _say("")
        _say_stop(index)
        _say("")
        if index + 1 < len(scenes):
            _say("next = %s   |   prev = back" % scenes[index + 1][1].lower())
        else:
            _say("This is the last stop.  prev = back  |  focus = zoom to any place "
                 "yourself")
        _say(report.LINE)

    worker.submit(lambda: None, lambda result, error: action())


def tour_next(*args, **kwargs):
    """
DESCRIPTION

    Go to the next stop of the tour started with "tour".

USAGE

    next
    """
    _step(+1, "next")


def tour_prev(*args, **kwargs):
    """
DESCRIPTION

    Go back one stop in the tour started with "tour".

USAGE

    prev
    """
    _step(-1, "prev")


# --------------------------------------------------------------------- quiz
QUIZ_OBJECT = "mystery"
FOLD_CLASSES = ["mainly alpha", "mainly beta", "alpha-beta"]


def _questions():
    return ["  1. What protein is this?",
            "  2. Is a ligand bound - and what is it?",
            "  3. Which fold class: mainly alpha, mainly beta, or alpha-beta (both)?"]


def quiz(*args, **kwargs):
    """
DESCRIPTION

    Loads a famous structure as the object "mystery" without naming it and
    asks three questions.  Look around, then check yourself.

USAGE

    quiz                 a new mystery structure
    quiz hint            a hint (up to three, each more revealing)
    quiz guess <name>    check your answer to question 1
    quiz reveal          show the answers, with explanations

EXAMPLE

    quiz
    quiz hint
    quiz guess myoglobin
    quiz reveal
    """
    words = base._raw(args, kwargs).split()
    verb = words[0].lower() if words else "new"
    if verb in ("new", "start", "next"):
        return _quiz_new()

    def action():
        base._header("quiz %s" % " ".join(words))
        entry = _quiz["entry"]
        if entry is None:
            _say("No quiz is running. Type  quiz  to get a mystery structure.")
        elif verb == "hint":
            hints = entry["hints"]
            if _quiz["hints"] >= len(hints):
                _say("No more hints. Make a guess with  quiz guess <name>  or see the "
                     "answers with  quiz reveal.")
            else:
                _quiz["hints"] += 1
                _say("Hint %d of %d: %s" % (_quiz["hints"], len(hints),
                                            hints[_quiz["hints"] - 1]))
        elif verb == "guess":
            guess = " ".join(words[1:]).lower().strip()
            if not guess:
                _say("Guess a name, e.g.  quiz guess myoglobin")
            elif any(key in guess for key in entry["accept"]):
                _say("Correct - it is %s. Now check the other two answers with  "
                     "quiz reveal." % entry["name"])
            elif guess in ("alpha", "beta", "mainly alpha", "mainly beta", "alpha-beta"):
                _say("That sounds like an answer to question 3. quiz guess checks question "
                     "1 (the protein's name);  quiz reveal  shows all three answers.")
            else:
                _say("Not that one. Try  quiz hint,  guess again, or  quiz reveal.")
        elif verb in ("reveal", "answer", "answers", "solution"):
            _quiz_reveal(entry)
        else:
            _say("Use  quiz,  quiz hint,  quiz guess <name>  or  quiz reveal.")
        _say(report.LINE)

    worker.submit(lambda: None, lambda result, error: action())


def _quiz_new():
    entries = _data("quiz.json")
    choices = [e for e in entries if e["pdb"] != _quiz["last"]] or entries
    entry = random.choice(choices)

    def job():
        try:
            return base._download_pdb(entry["pdb"]), None
        except http.NetworkError as err:
            return None, base._network_message(err)

    def done(answer, error):
        if error:
            return base._unexpected(error)
        path, trouble = answer
        base._header("quiz")
        if trouble:
            _say(trouble)
            _say(report.LINE)
            return
        cmd.delete(QUIZ_OBJECT)
        cmd.load(path, QUIZ_OBJECT)
        base._objects.pop(QUIZ_OBJECT, None)      # about/features must not give it away
        cmd.hide("everything", QUIZ_OBJECT)
        cmd.show("cartoon", "%s and polymer" % QUIZ_OBJECT)
        cmd.util.chainbow(QUIZ_OBJECT)
        cmd.show("sticks", "%s and organic" % QUIZ_OBJECT)
        cmd.show("nb_spheres", "%s and inorganic" % QUIZ_OBJECT)
        cmd.orient(QUIZ_OBJECT)
        teach.note("load <a structure file>, %s" % QUIZ_OBJECT,
                   "the file name would give the answer away")
        teach.note("hide everything, %s" % QUIZ_OBJECT)
        teach.note("show cartoon, %s and polymer" % QUIZ_OBJECT)
        teach.note("util.chainbow %s" % QUIZ_OBJECT, "rainbow along each chain")
        teach.note("show sticks, %s and organic" % QUIZ_OBJECT, "small molecules, if any")
        teach.note("orient %s" % QUIZ_OBJECT)
        _quiz.update(entry=entry, object=QUIZ_OBJECT, hints=0, last=entry["pdb"])
        residues, helix, strand = _secondary_structure(QUIZ_OBJECT)
        chains = len(cmd.get_chains("%s and polymer" % QUIZ_OBJECT))
        _say("A mystery structure is loaded as the object '%s': %d chain%s, %d residues."
             % (QUIZ_OBJECT, chains, "" if chains == 1 else "s", residues))
        _say("ON SCREEN    ribbons coloured blue (start of each chain) to red (its end); "
             "sticks = anything that is not protein")
        _say("")
        _say("Turn it around, look for spirals (helices), arrows (strands) and sticks, "
             "then answer:")
        for line in _questions():
            _say(line)
        _say("")
        _say("Useful while you look:  look chains %s  |  look ligand %s  |  count_atoms "
             "%s and organic" % ((QUIZ_OBJECT,) * 3))
        _say("When ready:  quiz hint  |  quiz guess <name>  |  quiz reveal")
        _say(report.LINE)

    worker.submit(job, done)


def _quiz_reveal(entry):
    obj = _quiz["object"]
    pdb_id = entry["pdb"]
    _say("Answers:")
    _say("  1. %s  (PDB entry %s)" % (entry["name"], pdb_id))
    if entry["ligand"]:
        _say("  2. Yes: %s - %s." % (entry["ligand"], entry["ligand_name"]))
    else:
        _say("  2. %s." % (entry["ligand_name"][:1].upper() + entry["ligand_name"][1:]))
    _say("  3. %s: %s." % (entry["fold"][:1].upper() + entry["fold"][1:], entry["fold_note"]))
    _say("")
    for fact in entry["facts"]:
        _say("  - " + fact)
    if obj in cmd.get_names("objects"):
        residues, helix, strand = _secondary_structure(obj)
        _say("  - In this model %.0f%% of the residues are in helices and %.0f%% in "
             "strands." % (helix, strand))
        if pdb_id not in cmd.get_names("objects"):
            cmd.set_name(obj, pdb_id)
            base._objects[pdb_id] = ("pdb", pdb_id)
            teach.note("set_name %s, %s" % (obj, pdb_id), "rename the object")
            _say("")
            _say("The object is now called %s. Explore it:  about %s  |  tour %s"
                 % (pdb_id, pdb_id, pdb_id))
    _say("Another one?  quiz")
    _quiz.update(entry=None, object=None, hints=0)


# ----------------------------------------------------------------- surprise
def surprise(*args, **kwargs):
    """
DESCRIPTION

    Loads a random structure from Prospector's list of about 75 classic,
    textbook structures and prints its info card.

USAGE

    surprise
    """
    classics = rank.classics()
    remaining = [p for p in sorted(classics) if p not in _surprised]
    if not remaining:
        del _surprised[:]
        remaining = sorted(classics)
    pdb_id = random.choice(remaining)
    _surprised.append(pdb_id)

    def job():
        try:
            path = base._download_pdb(pdb_id)
        except http.NetworkError as err:
            return None, None, base._network_message(err)
        try:
            text = card.entry_card(pdb_id)
        except http.NetworkError:
            text = None
        return path, text, None

    def done(answer, error):
        if error:
            return base._unexpected(error)
        path, text, trouble = answer
        base._header("surprise")
        if trouble:
            _say(trouble)
            _say(report.LINE)
            return
        if pdb_id not in cmd.get_names("objects"):
            cmd.load(path, pdb_id)
        base._objects[pdb_id] = ("pdb", pdb_id)
        teach.note("fetch %s" % pdb_id, "download and load the entry")
        _do("hide everything, %s" % pdb_id)
        _do("show cartoon, %s and polymer" % pdb_id)
        _do("util.cbc %s" % pdb_id, "one colour per chain")
        _do("show sticks, %s and organic" % pdb_id, "ligands as sticks")
        _do("show nb_spheres, %s and inorganic" % pdb_id, "metal ions")
        _do("orient %s" % pdb_id)
        _say("")
        _say("Surprise! %s - %s" % (pdb_id, classics[pdb_id]))
        _say("(picked at random from Prospector's list of classic structures)")
        _say("ON SCREEN    the whole structure: ribbons = protein, one colour per chain; "
             "sticks = bound molecules; small balls = metal ions")
        _say("")
        _say(text or "(The info card could not be fetched; try  about %s  later.)" % pdb_id)
        _say("")
        _say("Next:  tour %s   |   features %s   |   surprise" % (pdb_id, pdb_id))
        _say(report.LINE)

    worker.submit(job, done)


def register():
    table = {"tour": tour, "next": tour_next, "prev": tour_prev, "quiz": quiz,
             "surprise": surprise}
    for name, function in table.items():
        cmd.extend(name, function)
        cmd.extend("p_" + name, function)
    try:
        cmd.auto_arg[0]["tour"] = cmd.auto_arg[0]["zoom"]
        cmd.auto_arg[0]["quiz"] = [lambda: cmd.Shortcut(["hint", "guess", "reveal"]),
                                   "action", " "]
    except Exception:
        pass
