"""PyMOL commands.  Each one is registered with cmd.extend, so it works from
the PyMOL command line, from scripts and from the Prospector panel."""
from __future__ import annotations

import os
import re
import textwrap

from pymol import cmd

from ..core import alphafold as af_api
from ..core import cache, card, chemcomp, http, rcsb, report, search, store, uniprot
from ..core import features as feat
from . import teach, worker

VERSION = "2.9.3"
_listeners = []   # panel callbacks: function(command_name, search.Result)


def _say(text):
    """Print to the PyMOL console; never fail on consoles with odd encodings.
    The closing dashed line of a command also releases its teach-mode notes."""
    if text == report.LINE:
        notes = teach.take()
        if _held["on"]:
            _release(short=True)
        for line in notes:
            _print(line)
        _print(text)
        for listener in list(_output_listeners):     # the panel shows it too
            try:
                listener(_full["command"], list(_full["raw"]))
            except Exception:
                pass
        return
    _full["raw"].extend(str(text).split("\n"))        # unwrapped, for the panel
    text = "\n".join(_wrap(line) for line in str(text).split("\n"))
    if _held["on"]:
        _held["lines"].extend(text.split("\n"))
        return
    _full["lines"].extend(text.split("\n"))
    _print(text)


# ------------------------------------------------------------- short output
# In short mode the text of an action command is held back until its closing
# line, then only the answer is printed;  more  shows the rest.
_SHORTENED = set("grab alphafold focus scout pocket features whatif mutate measure compare "
                 "viewas look tour next prev back clear snap report surprise batch "
                 "molecule choose".split())
_KEEP_ALWAYS = ("NOTE", "CAUTION", "WARNING", "CHECK")
_KEEP_NEVER = ("MOUSE", "GO BACK", "WHAT NEXT", "WHEN DONE", "GO CLOSER", "Tip", "FILE",
               "OPEN LATER", "THE COMMAND")
_BLOCK = re.compile(r"^([A-Z][A-Za-z]*(?: [A-Za-z?]+){0,3}?)(?::|\s{2,})\s*\S")
_output_listeners = []    # function(command, lines), called when a command has finished
_held = {"on": False, "lines": []}
_full = {"command": "", "lines": [], "raw": [], "hidden": 0}


def short_mode():
    return bool(store.get_setting("short", False))


def _shorten(lines, blocks_kept=2, lines_each=2, always=()):
    """The answer-first part of a command's output: its notes, then its first
    two labelled blocks, two lines each.  Returns (kept lines, number hidden).
    The panel asks for a little more (``blocks_kept``, ``lines_each``) and for
    the verdict lines (``always``)."""
    lines = [line for line in lines if line.strip()]
    blocks, current = [], None
    for line in lines:
        match = _BLOCK.match(line)
        if match:
            current = [match.group(1).strip(), [line]]
            blocks.append(current)
        elif current is None or not line.startswith(" "):
            current = ["", [line]]
            blocks.append(current)
        else:
            current[1].append(line)
    if not any(label for label, _ in blocks):          # a plain message or a table
        kept = lines[:6]
        return kept, len(lines) - len(kept)
    kept, normal = [], 0
    for label, body in blocks:
        if label.startswith(_KEEP_NEVER):
            continue
        if label.startswith(_KEEP_ALWAYS) or (always and label.startswith(always)):
            kept.extend(body[:max(3, lines_each)])
        elif normal < blocks_kept or not label and normal < blocks_kept + 1 \
                and len(kept) < 2:
            kept.extend(body[:lines_each])
            normal += 1
    return kept, len(lines) - len(kept)


def _release(short):
    """Print what was held back for the running command."""
    lines, _held["lines"], _held["on"] = _held["lines"], [], False
    _full["lines"].extend(lines)
    if not short:
        for line in lines:
            _print(line)
        return
    kept, hidden = _shorten(lines)
    if hidden <= 2:                 # not worth hiding
        kept, hidden = [line for line in lines if line.strip()], 0
    for line in kept:
        _print(line)
    _full["hidden"] = hidden
    if hidden:
        _print("   (+%d more line%s:  more )" % (hidden, "" if hidden == 1 else "s"))


def more(*args, **kwargs):
    """
DESCRIPTION

    Show the full text of the last command.  Useful in short mode
    (prospector short on), where commands print only the answer.

USAGE

    more
    """
    def action():
        if not _full["lines"]:
            _print("Prospector: there is nothing more to show yet.")
            return
        _print("")
        _print("Prospector: %s   (full text)" % _full["command"])
        for line in _full["lines"]:
            _print(line)
        _print(report.LINE)

    worker.submit(lambda: None, lambda result, error: action())


_LABEL = re.compile(r"^(\s*[A-Z][A-Za-z ]{1,11}\s{2,}|\s*\[\d+/\d+\] |\s*-\s|\s*\d+\.\s)")


def _wrap(line, width=80):
    """Fold a line that is wider than the PyMOL console, keeping it aligned
    under its own text (so 'CHECK        long sentence' stays a tidy block)."""
    if len(line) <= width:
        return line
    lead = len(line) - len(line.lstrip())
    match = _LABEL.match(line)
    hang = len(match.group(0)) if match else (lead if lead else 3)
    return textwrap.fill(line.strip(), width - 2, initial_indent=" " * lead,
                         subsequent_indent=" " * min(hang, 30),
                         break_long_words=False, break_on_hyphens=False)


def _print(text):
    try:
        print(text)
    except UnicodeEncodeError:
        print(text.encode("ascii", "replace").decode("ascii"))


def _raw(args, kwargs, joiner=" "):
    """Rebuild what the user typed.  PyMOL splits arguments at commas and turns
    a leading `word=` into a keyword, so the pieces are joined back together."""
    parts = [str(a) for a in args]
    for key, value in kwargs.items():
        if key.startswith("_"):      # PyMOL passes _self internally
            continue
        parts.append("%s=%s" % (key, value))
    return joiner.join(p.strip() for p in parts if str(p).strip())


_choices = []          # what typing 1, 2, 3 ... runs: the last numbered list
_choice_labels = []    # and what each line was, to say what is being loaded


def _set_choices(result):
    try:
        found = report.choices(result)
        if result.kind == "ambiguous":
            labels = ["%s, %s" % (g.label(1), g.organisms_text(1).lower())
                      for g in result.resolution.choices]
        else:
            labels = [hit.entry.title[:50] for hit in result.hits]
    except Exception:
        found, labels = [], []
    if found:
        _choices[:] = found
        _choice_labels[:] = labels[:len(found)] + [""] * (len(found) - len(labels))


def number(index):
    """Typing 1, 2, 3 ... runs that line of the last numbered list."""
    if not _choices:
        _say("A number picks a line from the last numbered list, e.g. after  grab myosin"
             "  or  find hemoglobin.  There is no list yet.")
        return
    if not 1 <= index <= len(_choices):
        _say("The last list has %d line%s: type a number from 1 to %d."
             % (len(_choices), "" if len(_choices) == 1 else "s", len(_choices)))
        return
    label = _choice_labels[index - 1] if index <= len(_choice_labels) else ""
    _say("Prospector: %d  ->  %s%s" % (index, _choices[index - 1],
                                       "   (%s)" % label if label else ""))
    cmd.do(_choices[index - 1])


def _notify(command_name, result):
    if command_name in ("find", "grab"):
        _set_choices(result)
    for listener in list(_listeners):
        try:
            listener(command_name, result)
        except Exception:
            pass


def _unexpected(error):
    if _held["on"]:
        _release(short=False)
    _say("Prospector hit an unexpected internal error (PyMOL is fine). Details:\n" + error)


# --------------------------------------------------------------------- find
def find(*args, **kwargs):
    """
DESCRIPTION

    Find the best structure for a protein name, gene, UniProt ID or PDB ID.
    Prints a ranked table, explains why #1 won, and gives the command to run.
    Nothing is loaded - use "grab" for that.

USAGE

    find <query> [organism=human] [method=xray|em|nmr] [res<2.0]
                 [ligand=ATP] [year>2015] [limit=10] [mode=text]

EXAMPLES

    find hemoglobin
    find p53 organism=human res<2.0
    find spike protein sars-cov-2 method=em
    find kinase ligand=STI
    """
    raw = _raw(args, kwargs)

    def done(result, error):
        if error:
            return _unexpected(error)
        text = report.render(result, "find")
        if teach.enabled():
            closing = "\n" + report.LINE
            if text.endswith(closing):
                text = text[:-len(closing)]
            _say(text)
            _say("\nTeach mode: find only asks the databases (UniProt, RCSB) - nothing in "
                 "PyMOL changed.\n  The line after 'Run this:' is the PyMOL command that "
                 "loads the structure.")
            _say(report.LINE)
        else:
            _say(text)
        _notify("find", result)

    worker.submit(lambda: search.find(raw), done)


# --------------------------------------------------------------------- grab
def _download_pdb(pdb_id):
    path = os.path.join(cache.structures_dir(), "%s.cif" % pdb_id.upper())
    return http.download(rcsb.DOWNLOAD_URL % pdb_id.upper(), path)


def _download_model(model):
    path = os.path.join(cache.structures_dir(), "%s.cif" % model.object_name())
    return http.download(model.cif_url, path)


def _download_molecule(molecule):
    path = os.path.join(cache.structures_dir(), "ligand_%s.sdf" % molecule.identifier)
    return http.download(molecule.sdf_url(), path)


_last_click = [None]      # (object, chain, number, name) of the newest click


def _clicked_residue():
    """(object, chain, number, three-letter name) of the residue the user
    clicked in the 3D window - PyMOL puts a click into the selection "sele" -
    or None when nothing, or more than one residue, is selected."""
    try:
        if "sele" not in cmd.get_names("selections") or not cmd.count_atoms("sele"):
            return None
        if cmd.count_atoms("byres sele") != cmd.count_atoms("byres first sele"):
            # several residues are selected (each click adds one): the panel,
            # when open, knows which was clicked last
            last = _last_click[0]
            if last and cmd.count_atoms("sele and model %s and resi %s%s" % (
                    last[0], last[2], " and chain %s" % last[1] if last[1] else "")):
                return last
            return None
        found = []
        cmd.iterate("first sele", "found.append((model, chain, resi, resn))",
                    space={"found": found})
        return found[0] if found else None
    except Exception:
        return None


def _is_molecule(obj):
    """True for an object that is a small molecule on its own (no chain)."""
    identity = _objects.get(obj)
    if identity and identity[0] == "mol":
        return True
    return bool(obj) and obj in cmd.get_names("objects") and \
        cmd.count_atoms("%s and polymer" % obj) == 0 and cmd.count_atoms(obj) > 0


def _needs_protein(obj, command):
    """Commands made for protein chains call this first.  Prints why and
    returns True when ``obj`` is a small molecule."""
    identity = _objects.get(obj)
    if identity and identity[0] == "built":
        kind = identity[1]
        if "peptide" in kind or "sheet" in kind:
            if command in ("features", "similar", "tour"):
                _say("%s is a %s you built: it has no database record, so  %s  has nothing "
                     "to look up. about %s  describes it." % (obj, kind, command, obj))
                _say(report.LINE)
                return True
            return False
        if command in ("scout", "measure", "explain", "choose", "viewas", "look") \
                and "NA" in kind:
            return False
        _say("%s is a %s you built, not a protein, so  %s  has nothing to work on."
             % (obj, kind, command))
        _say("FOR IT       about %s  |  viewas %s spheres  |  measure  |  choose  |  snap"
             % (obj, obj))
        _say(report.LINE)
        return True
    if not _is_molecule(obj):
        return False
    _say("%s is a small molecule, not a protein, so  %s  has nothing to work on "
         "(it needs a chain of amino acids)." % (obj, command))
    _say("FOR A MOLECULE   about %s  |  viewas %s sticks  |  measure  |  choose  |  snap"
         % (obj, obj))
    identity = _objects.get(obj)
    if identity and identity[0] == "mol" and not identity[1].startswith("CID"):
        _say("IN A PROTEIN     grab protein ligand=%s   loads a structure that has it "
             "bound; then  pocket  shows what holds it" % identity[1])
    _say(report.LINE)
    return True


def _load_molecule(molecule, path, quiet=False):
    """Load a small molecule as sticks.  Returns the object name."""
    name = molecule.object_name()
    if name not in cmd.get_names("objects"):
        _snapshot("loading %s" % name, allow_empty=True)
        cmd.load(path, name)
    _objects[name] = ("mol", molecule.identifier)
    _last[0] = name
    cmd.hide("everything", name)
    cmd.show("sticks", name)
    cmd.show("spheres", name)
    cmd.set("sphere_scale", 0.22, name)
    cmd.set("stick_radius", 0.14, name)
    cmd.set("valence", 1, name)
    cmd.util.cbag(name)
    cmd.orient(name)
    teach.note("load %s, %s" % (molecule.sdf_url(), name), "the molecule's 3D file")
    teach.note("show sticks, %s" % name, "bonds as sticks")
    teach.note("util.cbag %s" % name, "colour by element, carbon green")
    if not quiet:
        _say("ON SCREEN    the molecule as balls and sticks: green = carbon, red = oxygen, "
             "blue = nitrogen, white = hydrogen, yellow = sulfur, orange = phosphorus; "
             "a double stick = a double bond")
        _say("             This is %s, not a measurement."
             % ("a computed shape from PubChem" if molecule.source == "PubChem"
                else "an idealised shape from the PDB's dictionary of molecules"))
    return name


def _load_pdb(pdb_id, path):
    name = pdb_id.upper()
    if name in cmd.get_names("objects"):
        _say("%s is already loaded." % name)
    else:
        _snapshot("loading %s" % name, allow_empty=True)
        cmd.load(path, name)
    _objects[name] = ("pdb", name)
    _last[0] = name
    cmd.zoom(name)
    teach.note("fetch %s" % name, "download this entry from the PDB and load it")
    teach.note("zoom %s" % name, "centre the view on it")


def _load_model(model, path):
    name = model.object_name()
    if name not in cmd.get_names("objects"):
        _snapshot("loading %s" % name, allow_empty=True)
        cmd.load(path, name)
    _objects[name] = ("af", model.accession)
    # AlphaFold stores its per-residue confidence (pLDDT) in the B-factor column.
    cmd.color("0x0053D6", "%s and b > 90" % name)
    cmd.color("0x65CBF3", "%s and b > 70 and not b > 90" % name)
    cmd.color("0xFFDB13", "%s and b > 50 and not b > 70" % name)
    cmd.color("0xFF7D45", "%s and not b > 50" % name)
    cmd.zoom(name)
    teach.note("load %s, %s" % (model.cif_url, name), "AlphaFold model straight from its URL")
    teach.note("color 0x0053D6, %s and b > 90" % name, "pLDDT is in the B-factor column (b)")
    teach.note("color 0x65CBF3, %s and b > 70 and not b > 90" % name, "confident")
    teach.note("color 0xFFDB13, %s and b > 50 and not b > 70" % name, "low")
    teach.note("color 0xFF7D45, %s and not b > 50" % name, "very low")
    confidence = ""
    if model.mean_plddt is not None:
        confidence = ", mean pLDDT %.0f" % model.mean_plddt
    _say("")
    _say("LOADED       AlphaFold model %s (object name: %s)" % (model.model_id, name))
    _say("             %s%s, %s" % (model.name, " (gene %s)" % model.gene if model.gene
                                    else "", model.organism))
    _say("             residues %s-%s%s" % (model.start, model.end, confidence))
    _say("             This is a computer PREDICTION, not an experiment.")
    _say("")
    _say("ON SCREEN    the predicted chain as a ribbon, coloured by how sure the "
         "prediction is:")
    bands = dict((label.split(" (")[0], value) for label, value in model.fractions)
    for colour, label, key in (("dark blue ", "very sure (pLDDT > 90)", "very high"),
                               ("light blue", "confident (70-90)", "confident"),
                               ("yellow    ", "unsure (50-70)", "low"),
                               ("orange    ", "very unsure (< 50), often floppy", "very low")):
        share = bands.get(key)
        _say("                %s  %-34s %s" % (colour, label, "%3.0f%% of residues"
                                              % (100 * share) if share is not None else ""))
    _say("             Trust the blue parts; treat yellow and orange as 'position "
         "unknown'.")
    _say("")
    _say("WHAT NEXT    about %s  |  features %s  |  tour %s" % (name, name, name))


def grab(*args, **kwargs):
    """
DESCRIPTION

    Find the best structure (same search and filters as "find"), load it, and
    print a one-line summary.  If the name is ambiguous nothing is loaded and
    the choices are shown.  If no experimental structure exists, the
    AlphaFold model is loaded instead.

USAGE

    grab <query> [filters]

EXAMPLES

    grab hemoglobin
    grab 1TUP
    grab p53 method=xray res<1.5
    """
    raw = _raw(args, kwargs)

    def job():
        result = search.find(raw)
        path, problem = None, None
        try:
            if result.kind == "results":
                path = _download_pdb(result.top.pdb_id)
            elif result.kind == "alphafold":
                path = _download_model(result.model)
            elif result.kind == "molecule":
                path = _download_molecule(result.molecule)
        except http.NetworkError as err:
            problem = str(err)
        return result, path, problem

    def done(answer, error):
        if error:
            return _unexpected(error)
        result, path, problem = answer
        if result.kind == "results" and path:
            top = result.top
            _header("grab %s" % raw)
            _say("")
            _say("LOADED       %s - %s" % (top.pdb_id, report._cut(top.entry.title, 58)))
            for line in report.facts(top, result):
                _say("             " + line)
            _load_pdb(top.pdb_id, path)
            call = top.pdb_id
            generic = result.query.text.lower() in ("protein", "proteins", "structure",
                                                    "structures", "anything", "any")
            if generic:      # "grab protein ligand=CLR": the word is not a name
                result.notes = [n for n in result.notes if "whole family" not in n]
            if result.mode != "id" and not generic:
                call = result.query.text
                _remember(raw, top.pdb_id)
                _remember(call, top.pdb_id)
                _say("             From now on you can call it '%s' (or %s)."
                     % (call, top.pdb_id))
            _say("")
            _say("ON SCREEN    the whole structure: ribbons = protein chains, sticks = "
                 "bound molecules, small crosses = water")
            _say("")
            _say("WHAT NEXT    about %s  (what it is + key places to zoom to)" % call)
            _say("             tour %s   (guided walk through it)" % call)
            _say("             features %s  |  pocket %s  |  look publication"
                 % (call, call))
            if len(result.hits) > 1:
                _say("")
                _say("Why this one? It ranked best of %d matching structures (next: %s)."
                     % (result.total, ", ".join(h.pdb_id for h in result.hits[1:4])))
                _say("See the whole list and the reasoning with:  find %s" % raw)
            if result.mode != "id" and not generic and result.query.kind == "auto" \
                    and not result.query.has_filters() \
                    and " " not in result.query.text.strip() and (
                        result.resolution is None or result.resolution.top_tier < 3):
                _say("MOLECULE?    if you meant the small molecule of that name:  "
                     "molecule %s" % result.query.text)
            for note in result.notes:
                _say("Note: " + note)
            _say(report.LINE)
        elif result.kind == "molecule" and path:
            molecule = result.molecule
            _header("grab %s" % raw)
            _say("")
            for line in report.molecule_lines(molecule, "LOADED       "):
                _say(line)
            name = _load_molecule(molecule, path)
            for text in (raw, result.query.text, molecule.title(), molecule.name):
                _remember(text, name)
            _say("             From now on you can call it '%s' (or %s)."
                 % (result.query.text, name))
            _say("")
            _say("WHAT NEXT    about %s  (the full card)  |  viewas %s spheres  (its real "
                 "bulk)" % (result.query.text, result.query.text))
            if molecule.entries:
                _say("             grab protein ligand=%s   a protein with it bound (%d "
                     "structures have it)" % (molecule.comp_id, molecule.entries))
            _say("NOT THIS?    proteins named after it:  find %s kind=protein"
                 % result.query.text)
            for note in result.notes:
                _say("Note: " + note)
            _say(report.LINE)
        elif result.kind == "alphafold" and path:
            _say(report.render(result, "grab"))
            _load_model(result.model, path)
            _remember(raw, result.model.object_name())
            _remember(result.query.text, result.model.object_name())
        else:
            _say(report.render(result, "grab"))
            if problem:
                _say("The structure file could not be downloaded: %s" % problem)
        _notify("grab", result)

    worker.submit(job, done)


# ----------------------------------------------------------------- molecule
def molecule(*args, **kwargs):
    """
DESCRIPTION

    Load a small molecule on its own - a steroid, lipid, sugar, nucleotide,
    drug or cofactor - from the PDB's dictionary of molecules, and show it
    as balls and sticks.  `grab <name>` does the same when the name belongs
    to a molecule and not to a protein; use `molecule` to insist.

    Only molecules that were seen in some PDB structure are available.
    Membranes and bilayers are not single molecules and cannot be loaded.

USAGE

    molecule <name | PDB ligand code>

EXAMPLES

    molecule cholesterol
    molecule ATP
    molecule caffeine
    molecule HEM
    """
    raw = _raw(args, kwargs)
    if not raw.strip():
        worker.submit(lambda: None, lambda r, e: (
            _header("molecule"), _say(""),
            _say("Type a molecule's name after it, e.g.  molecule cholesterol   "
                 "molecule ATP   molecule caffeine"), _say(report.LINE)))
        return
    grab(raw + " kind=molecule")


# ---------------------------------------------------------------- alphafold
def alphafold(*args, **kwargs):
    """
DESCRIPTION

    Load the AlphaFold predicted model for a protein, coloured by confidence
    (pLDDT), and print the colour legend.  Works for any protein in AlphaFold
    DB, whether or not an experimental structure exists.

USAGE

    alphafold <UniProt ID | protein name>

EXAMPLES

    alphafold P04637
    alphafold MYMK organism=human
    """
    raw = _raw(args, kwargs)

    def job():
        text = raw.strip()
        if not text:
            return "Give a UniProt ID or a protein name, e.g.  alphafold P04637", None, None
        try:
            if uniprot.ACCESSION_RE.match(text):
                accession = text.upper()
            else:
                result = search.find(text)
                if result.kind == "ambiguous" or result.kind in ("error", "suggest"):
                    return report.render(result, "alphafold"), None, None
                if result.protein is None:
                    return ("Could not identify one protein for '%s'. Try the gene name or "
                            "a UniProt ID, e.g.  alphafold P04637" % text), None, None
                accession = result.protein.accession
            model = af_api.prediction(accession)
            if model is None:
                return "AlphaFold DB has no model for %s." % accession, None, None
            return None, model, _download_model(model)
        except http.NetworkError as err:
            return "Network problem: %s" % err, None, None

    def done(answer, error):
        if error:
            return _unexpected(error)
        message, model, path = answer
        _header("alphafold %s" % raw)
        if message:
            _say(message)
        else:
            _load_model(model, path)
        _say(report.LINE)

    worker.submit(job, done)


# ------------------------------------------------- which structure is this?
_objects = {}   # object name -> ("pdb", PDB id) or ("af", UniProt accession)
_PDB_NAME = re.compile(r"^[1-9][A-Za-z0-9]{3}$")
_AF_NAME = re.compile(r"^AF[_-]([A-Za-z0-9]+)[_-]F\d+", re.I)


def _molecules(enabled_only=False):
    names = cmd.get_names("objects", enabled_only=1 if enabled_only else 0) or []
    return [n for n in names if cmd.get_type(n) == "object:molecule"]


def _identify(name):
    """("pdb", id) / ("af", accession) for an object, or None if unknown."""
    if name in _objects:
        return _objects[name]
    match = _AF_NAME.match(name)
    if match:
        return ("af", match.group(1).upper())
    if _PDB_NAME.match(name):
        return ("pdb", name.upper())
    return None


_undo = []          # [(scene name, what was about to happen, {object: transparency})]
_undo_count = [0]


def _snapshot(label, allow_empty=False):
    """Remember how the screen looks now, so `back` can return to it."""
    names = _molecules()
    if not names and not allow_empty:
        return
    _undo_count[0] += 1
    scene = None
    if names:
        scene = "_prospector_undo_%d" % _undo_count[0]
        cmd.scene(scene, "store")
    settings = dict((name, cmd.get("cartoon_transparency", name)) for name in names)
    settings["*grid"] = (cmd.get("grid_mode"),
                         dict((name, cmd.get("grid_slot", name)) for name in names))
    # everything that exists now: whatever appears later is removed by `back`
    settings["*names"] = set(cmd.get_names("all"))
    settings["*last"] = _last[0]
    _undo.append((scene, label, settings))
    while len(_undo) > 10:
        old = _undo.pop(0)[0]
        if old:
            cmd.scene(old, "clear")


def _leave_panels():
    """End the two-panel view of `whatif`: one panel, groups dissolved, the
    changed copies hidden (they stay loaded)."""
    cmd.set("grid_mode", 0)
    names = cmd.get_names("all")
    for name in ("whatif_bonds_before", "whatif_bonds_after"):
        if name in names:
            cmd.delete(name)
    for group in ("whatif_original", "whatif_changed"):
        if group in cmd.get_names("all"):
            members = cmd.get_object_list("(%s)" % group) or []
            for member in members:
                cmd.ungroup(member)
                if group == "whatif_changed":
                    cmd.disable(member)
            cmd.delete(group)


def _clean(obj):
    """Remove the highlights earlier commands left on the protein (sticks, balls,
    lines, name tags), so the next view shows one thing only.  Bound molecules
    and the ribbon are left alone; `back` can always return to the old view."""
    for rep in ("sticks", "spheres", "lines"):
        cmd.hide(rep, "%s and polymer" % obj)
    cmd.label(obj, "")


_backs = [0]


def _restore(silent=False):
    """Go back one step.  Returns the text to print."""
    if not _undo:
        return ("Nothing to undo: no Prospector command has changed the picture yet "
                "(or you are already back at the start).")
    scene, label, settings = _undo.pop()
    before = settings.pop("*names", None)
    last = settings.pop("*last", None)
    grid_mode, slots = settings.pop("*grid", ("0", {}))
    if str(grid_mode) in ("0", "off"):
        _leave_panels()
    removed = []
    if before is not None:          # take away what was added since
        for name in cmd.get_names("all"):
            if name not in before and not name.startswith("_prospector") \
                    and name in cmd.get_names("all"):
                cmd.delete(name)
                removed.append(name)
                _objects.pop(name, None)
        for alias in [a for a, o in _aliases.items() if o in removed]:
            del _aliases[alias]
        if _last[0] in removed:
            _last[0] = last if last in cmd.get_names("objects") else None
    loaded = _molecules()
    for name in loaded:
        cmd.label(name, "")
    if scene:
        cmd.scene(scene, "recall", animate=0)
    if str(grid_mode) in ("0", "off"):
        _leave_panels()                      # e.g. leave the two panels of `whatif`
    else:
        cmd.set("grid_mode", grid_mode)
    for name, value in settings.items():
        if name in loaded:
            cmd.set("cartoon_transparency", value, name)
            cmd.set("stick_radius", 0.25, name)
    if scene:
        cmd.scene(scene, "clear")
    cmd.deselect()
    if silent:
        return ""
    teach.note("scene <name>, recall", "Prospector stored a scene before each change")
    _backs[0] += 1
    text = ("WENT BACK    to how it looked before:  %s\n             %s"
            % (label, "back  again undoes one more step (%d left)." % len(_undo)
               if _undo else "That was the first change - nothing more to undo."))
    if removed:
        text += "\n             Removed what was added since: %s%s" % (
            ", ".join(removed[:6]), " ..." if len(removed) > 6 else "")
    if _backs[0] == 1:
        text += ("\n             Colours, shapes and zoom are restored, and anything added "
                 "since (structures, measurements, changed copies) is taken away.")
    return text


_aliases = {}   # what the user typed (lower case) -> the object it loaded


_last = [None]   # the object the most recent command worked on


def _remember(text, obj):
    """Remember that ``text`` (e.g. "kinase") loaded the object ``obj`` (e.g. 3HEC),
    so later commands can say "about kinase" instead of "about 3HEC"."""
    text = " ".join((text or "").lower().split())
    if obj:
        _last[0] = obj
    if text and obj:
        _aliases[text] = obj


def _known_object(typed):
    """The loaded object called ``typed``, or the one a search for ``typed``
    loaded earlier.  None if there is no such object."""
    wanted = " ".join(typed.lower().split())
    if not wanted:
        return None
    loaded = _molecules()
    for name in loaded:
        if name.lower() == wanted:
            return name
    target = _aliases.get(wanted)
    return target if target in loaded else None


def _pick_object(typed, command="about"):
    """Find the object the user means.  Returns (object name or None, message)."""
    loaded = _molecules()
    if typed:
        known = _known_object(typed)
        if known:
            _last[0] = known
            return known, ""
        hint = " Loaded objects: %s." % ", ".join(loaded) if loaded else \
            " Nothing is loaded yet - try:  grab %s" % typed
        return None, "There is no object called '%s'.%s" % (typed, hint)
    if not loaded:
        return None, "Nothing is loaded yet. Load a structure first, e.g.  grab hemoglobin"
    enabled = _molecules(enabled_only=True)
    if len(enabled) == 1:
        _last[0] = enabled[0]
        return enabled[0], ""
    if len(loaded) == 1:
        _last[0] = loaded[0]
        return loaded[0], ""
    if _last[0] in loaded:          # the structure you were just working on
        return _last[0], ""
    return None, ("Several objects are loaded (%s). Say which one, e.g.  %s %s"
                  % (", ".join(loaded), command, loaded[0]))


def _split_options(raw):
    """'1PPB STI cutoff=5' -> (['1PPB', 'STI'], {'cutoff': '5'})"""
    words, options = [], {}
    last_key, separator = None, " "
    for piece in re.split(r"(\s*,\s*|\s+)", raw.strip()):
        if not piece:
            continue
        if piece.strip() in ("", ","):
            separator = "," if "," in piece else " "
            continue
        if "=" in piece:
            key, value = piece.split("=", 1)
            last_key = key.strip().lower()
            options[last_key] = value.strip()
        elif last_key and separator == ",":
            options[last_key] += "," + piece     # types=active,binding
        else:
            words.append(piece)
            last_key = None
    return words, options


_recent = []        # the last few Prospector commands, for the AI's context


def _header(text):
    teach.clear()
    if _held["on"]:                 # an earlier command never closed: show it all
        _release(short=False)
    _recent.append(" ".join(text.split())[:120])
    del _recent[:-12]
    first = text.split()[0] if text.split() else ""
    if first != "more":
        _full["command"], _full["lines"], _full["hidden"] = " ".join(text.split()), [], 0
        _full["raw"] = []
    _print("")
    _print("\n".join(_wrap("Prospector: " + text).split("\n")))
    _held["on"] = first in _SHORTENED and short_mode()
    return
    _say("")
    _say("Prospector: " + text)


# ---------------------------------------- a name where an object is expected
def _fetch_query(query, strict=False):
    """Worker-thread half of "the user typed a protein name instead of a loaded
    object": find the best structure and download it.  Returns a dict with
    either "error" (text to print) or "identity" + what is needed to load it.
    Words that cannot be a name (a keyword, a number, "and") are never
    searched for: loading some protein called that would be a wrong guess."""
    from . import typed
    why = typed.not_a_name(query, strict)
    if why:
        return {"error": "%s, so I did not search for it.\nIf it IS a protein or a "
                         "molecule, load it by name:  grab %s" % (why[:1].upper() + why[1:],
                                                                  query)}
    result = search.find(query + " limit=3")
    try:
        if result.kind == "results":
            top = result.top
            return {"error": None, "identity": ("pdb", top.pdb_id), "name": top.pdb_id,
                    "query": query, "text": result.query.text,
                    "path": _download_pdb(top.pdb_id), "model": None,
                    "line": "Loaded the best match for '%s':  %s - %s"
                            % (query, top.pdb_id, report._cut(top.entry.title, 36))}
        if result.kind == "alphafold":
            model = result.model
            return {"error": None, "identity": ("af", model.accession),
                    "query": query, "text": result.query.text,
                    "name": model.object_name(), "path": _download_model(model),
                    "model": model,
                    "line": "'%s' has no experimental structure, so its AlphaFold "
                            "model is used." % query}
        if result.kind == "molecule":
            molecule = result.molecule
            return {"error": None, "identity": ("mol", molecule.identifier),
                    "query": query, "text": result.query.text,
                    "name": molecule.object_name(),
                    "path": _download_molecule(molecule), "model": None,
                    "molecule": molecule,
                    "line": "Loaded the small molecule '%s':  %s (%s)"
                            % (query, molecule.title(),
                               molecule.comp_id or "from PubChem")}
    except http.NetworkError as err:
        return {"error": _network_message(err)}
    _set_choices(result)
    return {"error": "'%s' is not a loaded object, so I searched for it:\n%s"
                     % (query, report.render(result, "find"))}


def _place_fetched(fetched):
    """GUI-thread half: load what _fetch_query downloaded.  Returns the object name."""
    name = fetched["name"]
    if fetched.get("molecule") is not None:
        _load_molecule(fetched["molecule"], fetched["path"], quiet=True)
        _remember(fetched["molecule"].title(), name)
    elif name not in cmd.get_names("objects"):
        _snapshot("loading %s" % name, allow_empty=True)
        cmd.load(fetched["path"], name)
    _objects[name] = fetched["identity"]
    _remember(fetched["query"], name)
    _remember(fetched["text"], name)
    _say(fetched["line"])
    if fetched["model"] is None:
        teach.note("fetch %s" % name, "download this entry from the PDB and load it")
    else:
        _load_model(fetched["model"], fetched["path"])     # colours by confidence
    _say("")
    return name


# -------------------------------------------------------------------- about
def about(*args, **kwargs):
    """
DESCRIPTION

    Print an info card for a structure: protein, gene, organism, function
    (from UniProt), method, resolution, R-free, every chain, ligands with
    their full names, and the primary citation with its DOI.

USAGE

    about [object | PDB ID]

EXAMPLES

    about            (the structure that is loaded)
    about 1PPB
    about haemoglobin   (finds and loads the best match first)
    """
    raw = _raw(args, kwargs).strip()
    obj, identity = None, None
    query = {}

    def before():
        nonlocal obj, identity
        obj, message = _pick_object(raw)
        identity = _identify(obj) if obj else None
        if obj is None and raw and _PDB_NAME.match(raw):
            identity = ("pdb", raw.upper())     # not loaded, but a card is still possible
        elif obj is None and raw:
            query["text"] = raw                 # a protein name: look it up
            return True
        elif obj is None:
            _header("about %s" % raw)
            _say(message)
            _say(report.LINE)
            return False
        if identity and identity[0] == "built":
            from . import build
            _header("about %s" % (raw or obj))
            _say("")
            for line in build.card(obj) or ["%s is a %s you built." % (obj, identity[1])]:
                _say(line)
            _say(report.LINE)
            return False
        if identity is None:
            _header("about %s" % raw)
            _say("I don't know which database entry the object '%s' comes from. Rename "
                 "it to its PDB ID, or load it with  grab <name>." % obj)
            _say(report.LINE)
            return False
        return True

    def job():
        nonlocal identity
        try:
            if query:
                # a name that is not loaded: find the best match and load it,
                # like every other command does
                fetched = _fetch_query(query["text"])
                if fetched["error"]:
                    return fetched["error"]
                query["fetched"] = fetched
                identity = fetched["identity"]
            if (obj or query.get("fetched")) and identity[0] != "mol":   # key places
                from . import places
                query["places"] = places.fetch_annotation(identity)
            if identity[0] == "mol":
                found = chemcomp.get(identity[1])
                if found is None:
                    return "No molecule with the identifier %s was found." % identity[1]
                if found.source == "PDB":
                    found.entries = chemcomp.count_entries(identity[1])
                return card.molecule_card(found)
            if identity[0] == "pdb":
                text = card.entry_card(identity[1])
                return text or "The PDB has no entry called %s." % identity[1]
            model = af_api.prediction(identity[1])
            if model is None:
                return "AlphaFold DB has no model for %s." % identity[1]
            return card.model_card(model)
        except http.NetworkError as err:
            return _network_message(err)

    def done(text, error):
        if error:
            return _unexpected(error)
        nonlocal obj
        _header("about %s" % (obj or raw or identity[1]))
        if query.get("fetched"):
            obj = _place_fetched(query["fetched"])
        _say(text)
        if obj and _is_molecule(obj):
            _say("")
            _say("WHAT NEXT    viewas %s spheres  (its real bulk)  |  viewas %s sticks  |  "
                 "measure  |  choose  |  snap" % (obj, obj))
        elif obj and obj in cmd.get_names("objects"):
            from . import places
            annotation, names = query.get("places") or (None, {})
            places.build(obj, annotation, names)
            _say("")
            _say(places.listing(obj))
            _say("")
            _say("WHAT NEXT    focus 1  (zoom to a place)  |  tour %s  (guided walk)" % obj)
        elif not obj:
            _say("")
            _say("This structure is not loaded. Load it to see and zoom to its key "
                 "places:   grab %s" % (raw or identity[1]))
        _say(report.LINE)

    worker.submit(job, done, before)


def _network_message(err):
    if err.offline:
        return ("You seem to be offline: %s This command needs the internet the first "
                "time it is used on a structure." % err)
    return "Network problem: %s Please try again in a moment." % err


# ----------------------------------------------------------------- features
def features(*args, **kwargs):
    """
DESCRIPTION

    Map UniProt annotations onto a loaded structure as named, coloured
    selections: active sites, binding sites, domains, disulfides, PTMs and
    disease variants.  UniProt residue numbers are converted to the numbers
    used in the structure with SIFTS, so offsets and insertion codes are
    handled.  Every mapped residue is checked against the atoms really there.

USAGE

    features [object] [types=active,binding,domain,disulfide,ptm,variant]

    types=all also adds "site" (cleavage sites and similar).

EXAMPLES

    features 1PPB
    features 4HHB types=binding,variant
    features thrombin    (finds and loads the best match first)
    """
    words, options = _split_options(_raw(args, kwargs, ","))
    # PyMOL splits "types=active,disulfide" at the comma; put stray type names back
    known = set(feat.GROUP_KEYS) | set(feat._ALIASES) | set(["all"])
    stray = [w for w in words[1:] if w.lower() in known]
    if stray:
        words = [w for w in words if w not in stray]
        options["types"] = ",".join([options.get("types", "")] + stray).strip(",")
    typed = " ".join(words)
    label = "features %s" % " ".join(words + ["%s=%s" % kv for kv in options.items()])
    obj, identity, wanted = None, None, None
    pending = {}

    def before():
        nonlocal obj, identity, wanted
        obj, problem = _pick_object(typed, "features")
        if obj is None and typed:
            try:                               # a protein name: find and load it first
                wanted = feat.parse_types(options.get("types", ""))
                pending["query"] = typed
                return True
            except feat.FeatureError as err:
                problem = str(err)
        if obj is not None:
            identity = _identify(obj)
            if identity is None:
                problem = ("I don't know which database entry the object '%s' comes from, "
                           "so I cannot fetch its annotations. Load it with  grab <name>."
                           % obj)
            else:
                try:
                    wanted = feat.parse_types(options.get("types", ""))
                    problem = ""
                except feat.FeatureError as err:
                    problem = str(err)
        if problem:
            _header(label)
            _say(problem)
            _say(report.LINE)
            return False
        return True

    def job():
        nonlocal identity
        try:
            if pending:
                pending["fetched"] = _fetch_query(pending["query"])
                if pending["fetched"]["error"]:
                    return None, pending["fetched"]["error"]
                identity = pending["fetched"]["identity"]
            if identity[0] == "pdb":
                return feat.annotate(identity[1], wanted), None
            return feat.annotate_model(identity[1], wanted), None
        except feat.FeatureError as err:
            return None, str(err)
        except http.NetworkError as err:
            return None, _network_message(err)

    def done(answer, error):
        if error:
            return _unexpected(error)
        nonlocal obj
        annotation, trouble = answer
        _header(label)
        if pending.get("fetched") and not pending["fetched"]["error"]:
            obj = _place_fetched(pending["fetched"])
        if _needs_protein(obj, "features"):
            return
        if trouble:
            _say(trouble)
            _say(report.LINE)
            return
        _snapshot(label)
        _apply_features(obj, annotation)
        _say("")
        _say("DETAILS  (UniProt position -> chain/residue in this structure)")
        _say(feat.render(annotation, obj))
        _say("")
        _say("Residue numbers converted with SIFTS (source: %s)." % annotation.source)
        _say(report.LINE)

    worker.submit(job, done, before)


def _apply_features(obj, annotation):
    # 1. what is really in the object?  (chain, residue label) -> residue name
    present = {}
    cmd.iterate("%s and polymer and name CA" % obj, "present[(chain, resi)] = resn",
                space={"present": present})

    # 2. selections + colours; domains first so that sites are drawn on top
    plan = annotation.selections()
    plan.sort(key=lambda item: 0 if item[5] == "cartoon" else 1)
    created = []
    for suffix, key, title, placed, colour, style in plan:
        name = "%s_%s" % (obj, suffix)
        cmd.select(name, feat.selection_text(obj, placed), enable=0)
        residues = cmd.count_atoms("%s and name CA" % name)
        if residues == 0:
            cmd.delete(name)
            continue
        expression = feat.selection_text(obj, placed)
        if len(expression) > 70:
            expression = expression[:66] + "...)"
        teach.note("select %s, %s" % (name, expression), title)
        if style == "cartoon":
            cmd.color(colour, "%s and elem C" % name)
            teach.note("color %s, %s and elem C" % (colour, name), "colour the carbons")
        elif style == "spheres":
            cmd.show("spheres", "%s and name CA" % name)
            cmd.set("sphere_scale", 0.6, "%s and name CA" % name)
            cmd.color(colour, "%s and name CA" % name)
            teach.note("show spheres, %s and name CA" % name, "one ball per residue")
            teach.note("color %s, %s and name CA" % (colour, name))
        else:
            cmd.show("sticks", "%s and (sidechain or name CA)" % name)
            cmd.color(colour, "%s and elem C and (sidechain or name CA)" % name)
            teach.note("show sticks, %s and (sidechain or name CA)" % name,
                       "side chains as sticks")
            teach.note("color %s, %s and elem C and (sidechain or name CA)" % (colour, name))
        shape = {"cartoon": "ribbon", "spheres": "balls"}.get(style, "sticks")
        created.append((name, residues, "%s %s" % (colour, shape), title))
    cmd.deselect()

    if created:
        _say("")
        _say("HIGHLIGHTED ON %s" % obj)
        _say("   %-18s %-26s %-9s %s" % ("Drawn as", "What", "Residues", "Selection name"))
        for name, residues, colour, title in created:
            _say("   %-18s %-26s %-9d %s" % (colour, title[:26], residues, name))
        _say("ON SCREEN    the structure with these parts coloured; everything else is "
             "unchanged")
        _say("ZOOM IN      focus active site  |  focus variants  |  focus  (all places)")
        _say("             or use a selection name in any PyMOL command:  zoom %s"
             % created[-1][0])
    else:
        _say("")
        _say("NOTHING TO HIGHLIGHT: none of the annotated features fall inside the part "
             "of the\nprotein present in this object.")

    # 3. verify the mapping against the atoms
    checks = annotation.checks()
    if checks:
        match, absent, differ = 0, 0, []
        for chain, label, expected, position in checks:
            resn = present.get((chain, label))
            if resn is None:
                absent += 1
            elif not expected or feat.THREE_TO_ONE.get(resn.upper()) == expected:
                match += 1
            else:
                differ.append("%s/%s is %s (UniProt %d: %s)"
                              % (chain, label, resn, position, expected))
        _say("")
        _say("CHECK        %d of %d mapped residues are in the model with the amino acid "
             "UniProt expects." % (match, len(checks)))
        if absent:
            _say("             %d have no atoms here (disordered, or not built by the "
                 "authors)." % absent)
        if differ:
            _say("             %d differ - usually engineered mutations: %s%s"
                 % (len(differ), "; ".join(differ[:4]), " ..." if len(differ) > 4 else ""))


# ------------------------------------------------------------------- pocket
def _ligands_in(obj):
    """[(chain, resn, resi, heavy atoms)] for everything that is not polymer or water."""
    counts = {}
    cmd.iterate("%s and not polymer and not solvent and not elem H" % obj,
                "counts[(chain, resn, resi)] = counts.get((chain, resn, resi), 0) + 1",
                space={"counts": counts})
    found = [(chain, resn, resi, n) for (chain, resn, resi), n in counts.items()]
    found.sort(key=lambda item: (item[1] in rcsb.ADDITIVES, -item[3], item[0], item[2]))
    return found


def _residue_selection(obj, chain, resn, resi):
    label = ("\\" + resi) if resi.startswith("-") else resi
    chain_part = "chain %s and " % chain if chain else ""
    return "%s and %sresn %s and resi %s" % (obj, chain_part, resn, label)


def pocket(*args, **kwargs):
    """
DESCRIPTION

    Show a ligand-binding pocket: list the ligands, select the protein
    residues within 4 A of one, draw the polar contacts and zoom in.

USAGE

    pocket [object] [ligand code] [cutoff=4.0]

EXAMPLES

    pocket 1IEP            (picks the largest real ligand)
    pocket 1IEP STI
    pocket 4HHB HEM cutoff=5
    pocket haemoglobin   (finds and loads the best match first)
    """
    raw = _raw(args, kwargs)
    words, options = _split_options(raw)
    label = "pocket %s" % raw
    obj, identity, code, cutoff = None, None, "", 4.0
    pending = {}

    def before():
        nonlocal obj, identity, code, cutoff
        # what was typed?  [object or protein name] [ligand code]
        typed, code = "", ""
        if words:
            everything, but_last = " ".join(words), " ".join(words[:-1])
            if _known_object(everything):
                typed = everything
            elif len(words) > 1 and _known_object(but_last):
                typed, code = but_last, words[-1].upper()
            elif len(words) == 1 and _molecules() and \
                    re.match(r"^[A-Za-z0-9]{1,3}$", words[0]):
                code = words[0].upper()         # just a ligand code, e.g. "pocket HEM"
            else:                               # a protein name: find and load it first
                name_words = list(words)
                if len(name_words) > 1 and re.match(r"^[A-Z0-9]{1,3}$", name_words[-1]):
                    code = name_words.pop()
                pending["query"] = " ".join(name_words)
        obj, message = _pick_object(typed, "pocket")
        if pending:
            obj, message = "pending", ""
        try:
            cutoff = float(options.get("cutoff", 4.0))
            if not 1.0 <= cutoff <= 15.0:
                raise ValueError
        except ValueError:
            obj, message = None, ("cutoff must be a distance in Angstrom between 1 and "
                                  "15, e.g. cutoff=4.0")
        if obj is None:
            _header(label)
            _say(message)
            _say(report.LINE)
            return False
        identity = None if pending else _identify(obj)
        return True

    def job():   # ligand full names: nice to have, never essential
        nonlocal identity
        if pending:
            pending["fetched"] = _fetch_query(pending["query"])
            if pending["fetched"]["error"]:
                return {}
            identity = pending["fetched"]["identity"]
        if identity and identity[0] == "pdb":
            try:
                return rcsb.ligand_names(rcsb.entry_detail(identity[1]))
            except http.NetworkError:
                return {}
        return {}

    def done(names, error):
        nonlocal obj
        if error:
            return _unexpected(error)
        _header(label)
        if pending:
            if pending["fetched"]["error"]:
                _say(pending["fetched"]["error"])
                _say(report.LINE)
                return
            obj = _place_fetched(pending["fetched"])
        if _needs_protein(obj, "pocket"):
            return
        _snapshot(label)
        _show_pocket(obj, code, cutoff, names or {})
        _say(report.LINE)

    worker.submit(job, done, before)


def _show_pocket(obj, code, cutoff, names):
    ligands = _ligands_in(obj)
    if not ligands:
        _say("%s contains no ligands (nothing besides protein/nucleic acid and water)."
             % obj)
        return
    seen = {}
    for chain, resn, resi, atoms in ligands:
        seen.setdefault(resn, []).append("%s/%s" % (chain or "-", resi))
    sizes = dict((resn, atoms) for chain, resn, resi, atoms in ligands)
    listing = ["ALL LIGANDS IN %s" % obj]
    for resn in seen:
        name = names.get(resn, ("", ""))[0]
        note = "  [buffer/ion/additive]" if resn in rcsb.ADDITIVES else ""
        where = ", ".join(seen[resn][:3]) + (" +%d" % (len(seen[resn]) - 3)
                                            if len(seen[resn]) > 3 else "")
        listing.append("   %-4s %3d atoms  at %-18s %s%s" % (resn, sizes[resn], where,
                                                           name[:36], note))

    if code:
        chosen = [item for item in ligands if item[1].upper() == code]
        if not chosen:
            _say("There is no ligand '%s' in %s. Choose one of: %s"
                 % (code, obj, ", ".join(seen)))
            _say("")
            for line in listing:
                _say(line)
            return
    else:
        chosen = ligands           # already sorted: real ligands first, biggest first
    chain, resn, resi, atoms = chosen[0]

    ligand_name = "%s_ligand" % obj
    pocket_name = "%s_pocket" % obj
    polar_name = "%s_polar" % obj
    teach.note("select %s, %s" % (ligand_name, _residue_selection(obj, chain, resn, resi)),
               "the ligand")
    teach.note("select %s, byres (%s and polymer within %g of %s)"
               % (pocket_name, obj, cutoff, ligand_name),
               "whole residues with any atom within %g A" % cutoff)
    teach.note("show sticks, %s or %s" % (ligand_name, pocket_name))
    teach.note("color yellow, %s and elem C" % ligand_name, "ligand carbons yellow")
    teach.note("distance %s, %s, %s, mode=2" % (polar_name, ligand_name, pocket_name),
               "mode=2 draws polar contacts")
    teach.note("zoom %s or %s, 3" % (ligand_name, pocket_name), "3 A of margin")
    cmd.delete(polar_name)
    cmd.select(ligand_name, _residue_selection(obj, chain, resn, resi), enable=0)
    cmd.select(pocket_name, "byres (%s and polymer within %g of %s)"
               % (obj, cutoff, ligand_name), enable=0)
    cmd.show("sticks", "%s or %s" % (ligand_name, pocket_name))
    cmd.color("yellow", "%s and elem C" % ligand_name)
    if cmd.count_atoms(pocket_name):
        cmd.distance(polar_name, ligand_name, pocket_name, mode=2)
        cmd.hide("labels", polar_name)
    cmd.zoom("%s or %s" % (ligand_name, pocket_name), 3)
    cmd.deselect()

    residues = []
    cmd.iterate("%s and name CA" % pocket_name, "residues.append((chain, resn, resi))",
                space={"residues": residues})
    full = names.get(resn, ("", ""))[0]
    contacts = _polar_contacts(ligand_name, pocket_name)
    # 1. what you are looking at
    _say("")
    _say("NOW SHOWING  the pocket around %s in chain %s of %s" % (resn, chain or "-", obj))
    if full:
        _say("             %s = %s" % (resn, full[:60]))
    _say("             %d residues lie within %g A of it; %d possible hydrogen bond%s"
         % (len(residues), cutoff, len(contacts), "" if len(contacts) == 1 else "s"))
    _say("             yellow sticks = ligand, other sticks = pocket, dashes = polar "
         "contacts")
    if resn in rcsb.ADDITIVES:
        _say("             (This is only a buffer component / ion: the structure has no "
             "real ligand.)")
    # 2. the details
    _say("")
    _say("POCKET RESIDUES")
    by_chain = {}
    for res_chain, res_name, res_number in residues:
        item = "%s%s" % (res_name, res_number)
        if item not in by_chain.setdefault(res_chain, []):   # two shapes = one residue
            by_chain[res_chain].append(item)
    residues = [(c, i) for c, items in by_chain.items() for i in items]
    for res_chain in sorted(by_chain):
        _say(textwrap.fill(" ".join(by_chain[res_chain]), 78,
                           initial_indent="   chain %s: " % (res_chain or "-"),
                           subsequent_indent="            "))
    if not residues:
        _say("   none within %g A - try a larger distance, e.g. cutoff=6" % cutoff)
    _say("")
    if contacts:
        _say("POSSIBLE HYDROGEN BONDS (N or O atoms closer than 3.5 A)")
        for distance, text in contacts[:12]:
            _say("   %s  %.2f A" % (text, distance))
        if len(contacts) > 12:
            _say("   ... +%d more" % (len(contacts) - 12))
    else:
        _say("POSSIBLE HYDROGEN BONDS: none (no N or O atoms closer than 3.5 A)")
    _say("")
    for line in listing:
        _say(line)
    others = [r for r in seen if r != resn]
    if others:
        _say("   Show another one with:  pocket %s %s" % (obj, others[0]))


def _polar_contacts(ligand_sel, pocket_sel, limit=3.5):
    polar = "(elem N or elem O)"
    first = cmd.get_model("%s and %s" % (ligand_sel, polar)).atom
    second = cmd.get_model("%s and %s" % (pocket_sel, polar)).atom
    found = []
    for a in first:
        for b in second:
            d2 = sum((a.coord[i] - b.coord[i]) ** 2 for i in range(3))
            if d2 <= limit * limit:
                found.append((d2 ** 0.5, "%s %-4s - %s%s %s" % (a.resn, a.name, b.resn,
                                                               b.resi, b.name)))
    found.sort()
    return found



# --------------------------------------------------------------- prospector
_HELP = """Prospector %s - work with protein structures by name

 NEW HERE?
   start                   first steps and how to use the mouse
   ask <anything>          say it in your own words, e.g.
                           ask where does the drug bind in p53
   askai <anything>        the same, but straight to the AI (needs a key)
   status                  what is loaded and how it is drawn right now
   glossary <word>         what a term means (RMSD, resolution, chain ...)

 GET A STRUCTURE
   grab <name>             find the best structure and load it
   find <name>             only list the candidates (pick=2 loads row 2)
   alphafold <name>        the predicted model, when no experiment exists
   molecule <name>         a small molecule on its own: cholesterol, ATP, caffeine
   build <what> ...        make a model: dna, rna, helix, sheet, membrane, micelle,
                           or the whole biological assembly of a protein
   batch <file.txt>        load a list of names
   surprise                a random classic structure

 UNDERSTAND IT
   about <name>            what it is, and its key places
   tour <name>             guided walk (then: next, prev, tour 3)
   focus <place>           zoom to a key place, chain or residue; focus all
   explain <A/57>          everything about one residue
   sequence <name>         the sequence with helices and strands marked
   features <name>         colour the annotated sites and domains
   pocket <name>           the bound molecule and what holds it
   scout <name> <what>     spot things: histidine, disulfides, salt bridges
   choose <description>    select in plain words: histidines within 5 of the zinc

 EXPERIMENT
   whatif <name> <R248Q>   change a residue: original and changed side by side
   mutate <name> <R248Q>   the same change in a single view
   measure <A/57> <A/195>  a distance, drawn and explained
   compare <a> <b>         superpose two structures, colour the differences
   similar <name>          related proteins by sequence

 CHANGE THE PICTURE
   viewas <style>          lines, sticks, ballstick, spheres, cartoon, surface
   look <preset>           colour schemes: chains, rainbow, hydrophobicity ...
   back                    go back one step (back 3 = three steps)
   clear                   empty the 3D window (back restores it); clear text
   snap <file>             save a high-resolution picture

 KEEP YOUR WORK
   session save <name>     save everything (later: session open <name>)
   report <name>           one page with pictures, to print or hand in
   history                 your past searches
   fav                     your favourite structures (fav add, fav rm)

 PRACTISE
   quiz                    guess a mystery structure (quiz hint, quiz reveal)
   prospector teach on     show the PyMOL commands behind every action
   prospector short on     less text: only the answer (then: more)
   prospector ai           optional: let an AI read what the rules miss

 MORE
   info <command>          details and examples, e.g.  info scout
   info filters            narrowing a search: organism=human res<2.0 ...
   prospector panel        the same things as buttons

 <name> is a protein name (p53), a PDB ID (4HHB), or nothing at all for the
 structure you are working on."""

_FILTER_HELP = """FILTERS  narrow down  find  and  grab . Add them after the name,
separated by spaces:

   organism=human          only this organism (human, mouse, yeast, "e. coli")
   method=xray             how it was solved: xray, em (cryo-EM) or nmr
   res<2.0                 resolution better than 2.0 A (smaller = sharper)
   ligand=ATP              must contain this bound molecule (PDB 3-letter code)
   year>2015               released after 2015 (also year<2000)
   limit=5                 how many rows  find  shows (default 10)
   pick=2                  take row 2 of the list instead of the top one
   mode=text               search the text of all entries instead of by protein

EXAMPLES
   find p53 organism=human res<2.0
   grab kinase ligand=STI
   grab thrombin pick=2"""


def prospector(*args, **kwargs):
    """
DESCRIPTION

    Prospector settings and help.

USAGE

    prospector [help | panel | teach on|off | cache info | cache clear | async on|off]
    prospector ai [key | on | off | test | model <name> | forget | learned | unlearn]
    prospector short on|off      print only the answer; `more` shows the rest
    """
    words = _raw(args, kwargs).lower().split()

    def action():   # runs in turn, so "teach on" affects only the commands after it
        if not words or words[0] in ("help", "version"):
            _say(_HELP % VERSION)
        elif words[0] == "panel":
            open_panel()
        elif words[0] == "cache" and words[1:2] == ["clear"]:
            _say("Prospector: removed %d cached answers." % cache.clear())
        elif words[0] == "cache":
            info = cache.info()
            _say("Prospector cache: %d answers, %.1f MB, in %s\nDownloaded structures: %s"
                 % (info["files"], info["bytes"] / 1e6, info["path"], cache.structures_dir()))
        elif words[0] == "teach" and words[1:2] in (["on"], ["off"]):
            teach.set_enabled(words[1] == "on")
            _say("Prospector: teach mode %s.%s" % (
                words[1], " Every action now prints the PyMOL commands it ran, with comments."
                if words[1] == "on" else ""))
        elif words[0] == "teach":
            _say("Prospector: teach mode is %s. Use  prospector teach on  or  prospector "
                 "teach off." % ("on" if teach.enabled() else "off"))
        elif words[0] == "short" and words[1:2] in (["on"], ["off"]):
            store.set_setting("short", words[1] == "on")
            _say("Prospector: short mode %s.%s" % (
                words[1], " Commands now print only the answer; type  more  for the "
                "full text of the last one." if words[1] == "on" else
                " Commands print their full explanation again."))
        elif words[0] == "short":
            _say("Prospector: short mode is %s. Use  prospector short on  or  prospector "
                 "short off." % ("on" if short_mode() else "off"))
        elif words[0] == "ai":
            from . import aiask
            aiask.settings(words[1:])
        elif words[0] == "async" and words[1:2] in (["on"], ["off"]):
            worker.set_async(words[1] == "on")
            _say("Prospector: background mode %s." % words[1])
        else:
            _say("Prospector: I don't know '%s'. Type  prospector help" % " ".join(words))

    worker.submit(lambda: None, lambda result, error: action())


def info(*args, **kwargs):
    """
DESCRIPTION

    List every Prospector command with one line each, or explain one command
    in detail.

USAGE

    info                 the list
    info <command>       details and examples for one command
    info filters         how to narrow a search

EXAMPLES

    info scout
    info filters
    """
    topic = _raw(args, kwargs).strip().lower()

    def action():
        _say("")
        if not topic:
            _say(_HELP % VERSION)
        elif topic in ("filters", "filter"):
            _say(_FILTER_HELP)
        else:
            name = topic.split()[0]
            name = {"next": "tour", "prev": "tour", "teach": "prospector"}.get(name, name)
            function = cmd.keyword.get(name, [None])[0]
            text = getattr(function, "__doc__", None) if function else None
            ours = getattr(function, "__module__", "").startswith(__name__.rsplit(".", 1)[0])
            if not text or not ours:
                import difflib
                known = [k for k, v in cmd.keyword.items()
                         if getattr(v[0], "__module__", "").startswith(
                             __name__.rsplit(".", 1)[0]) and not k.startswith("p_")]
                close = difflib.get_close_matches(name, known, n=3, cutoff=0.6)
                _say("'%s' is not a Prospector command.%s" % (
                    topic, " Did you mean: %s?" % ", ".join(close) if close else ""))
                _say("Type  info  for the list, or  ask <what you want>.")
            else:
                _say("COMMAND  %s" % name)
                for line in text.strip("\n").splitlines():
                    line = line.rstrip()
                    if line.strip() in ("DESCRIPTION", "USAGE", "EXAMPLES", "EXAMPLE",
                                        "EXAMPLE FILE"):
                        _say("")
                        _say({"DESCRIPTION": "WHAT IT DOES", "USAGE": "HOW TO TYPE IT"}
                             .get(line.strip(), line.strip()))
                    elif line.strip():
                        _say("   " + line.strip())
        _say(report.LINE)

    worker.submit(lambda: None, lambda result, error: action())


# -------------------------------------------------------------------- panel
_panel = None


def open_panel():
    """Open (or raise) the Prospector panel."""
    global _panel
    try:
        from . import panel
    except Exception as err:
        _say("Prospector: the panel needs PyMOL's Qt interface (%s). "
             "All commands still work from the command line." % err)
        return None
    if _panel is None:
        _panel = panel.Panel()
        _listeners.append(_panel.on_result)
        _output_listeners.append(_panel.on_output)
        _dock[0] = panel.dock(_panel)       # inside PyMOL's window, if there is one
        _panel.holder = _dock[0]
        _panel.set_size(big=False)
    if _dock[0] is not None:
        _dock[0].show()
        _dock[0].raise_()
    else:
        _panel.show()
        _panel.raise_()
    return _panel


_dock = [None]


def register():
    """Register every command with PyMOL (plain name and p_ alias)."""
    commands = {"find": find, "grab": grab, "alphafold": alphafold, "about": about,
                "features": features, "pocket": pocket, "prospector": prospector, "more": more, "molecule": molecule,
                "info": info}
    for name, function in commands.items():
        cmd.extend(name, function)
        cmd.extend("p_" + name, function)
    for index in range(1, 21):        # 1, 2, 3 ... pick from the last numbered list
        cmd.extend(str(index), lambda *args, _index=index, **kwargs: number(_index))
    try:   # tab-completion for the sub-commands of `prospector`
        cmd.auto_arg[0]["prospector"] = [
            lambda: cmd.Shortcut(["help", "panel", "teach", "cache", "async"]),
            "subcommand", " "]
    except Exception:
        pass
    try:   # tab-completion of object names
        for name in ("about", "features", "pocket"):
            cmd.auto_arg[0][name] = cmd.auto_arg[0]["zoom"]
            cmd.auto_arg[0]["p_" + name] = cmd.auto_arg[0]["zoom"]
    except Exception:
        pass
    from . import ask, places, tier3, tier4, tools
    tier3.register()
    tier4.register()
    places.register()
    tools.register()
    from . import choose
    choose.register()
    from . import build
    build.register()
    ask.register()
    _tidy_everything()
    worker.init()


def _tidy_everything():
    """Every command first has what was typed tidied (see typed.py), so the
    same request works in the usual spellings and word orders."""
    import functools
    from . import typed

    def wrap(name, function):
        joiner = " , " if name == "compare" else " "

        @functools.wraps(function)
        def tidied(*args, **kwargs):
            text = _raw(args, kwargs, joiner)
            try:
                text = typed.tidy(name, text)
            except Exception:
                pass                       # never let tidying break a command
            return function(text) if text else function()
        tidied._prospector_tidied = True
        return tidied

    for name in typed._PLAIN + ("compare",):
        for key in (name, "p_" + name):
            entry = cmd.keyword.get(key)
            if entry and not getattr(entry[0], "_prospector_tidied", False):
                entry[0] = wrap(name, entry[0])
