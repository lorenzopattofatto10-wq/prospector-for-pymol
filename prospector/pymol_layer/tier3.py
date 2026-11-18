"""Tier 3 commands: compare, similar, batch, look, snap, history, fav."""
from __future__ import annotations

import os
import re
import textwrap

from pymol import cmd

from ..core import cache, http, rcsb, report, search, store
from ..core import features as feat
from ..core import similar as sim
from . import commands as base
from . import teach, worker

_say = base._say
_titles = {}     # PDB id -> title, remembered from search results (for `fav`)


def _on_result(command_name, result):
    """Listener: every find/grab goes into the history."""
    top = result.top.pdb_id if result.kind == "results" and result.top else ""
    if result.kind == "alphafold" and result.model:
        top = result.model.model_id
    for hit in result.hits:
        _titles[hit.pdb_id] = hit.entry.title
    if result.kind != "error" or result.raw:
        store.add_history(command_name, result.raw, result.kind, top)


def _in_order(action):
    """Run a quick GUI-thread action after all earlier commands have finished."""
    worker.submit(lambda: None, lambda result, error: action() if not error
                  else base._unexpected(error))


def _protein_ca(obj, chain=""):
    selection = "%s and polymer and name CA and alt +A" % obj
    if chain:
        selection += " and chain %s" % chain
    return selection


# ------------------------------------------------------------------ compare
def _split_compare(raw):
    options = dict((k.lower(), v) for k, v in re.findall(r"(\w+)\s*=\s*(\S+)", raw))
    rest = re.sub(r"\w+\s*=\s*\S+", " ", raw).strip().strip(",").strip()
    if "," in rest:
        parts = [p.strip() for p in rest.split(",") if p.strip()]
    elif re.search(r"\s+vs\.?\s+", rest, re.I):
        parts = [p.strip() for p in re.split(r"\s+vs\.?\s+", rest, flags=re.I) if p.strip()]
    else:
        parts = rest.split()
    return parts, options


def _split_chain(text):
    match = re.match(r"^(.+?):([A-Za-z0-9]{1,4})$", text.strip())
    if match:
        return match.group(1).strip(), match.group(2)
    return text.strip(), ""


def compare(*args, **kwargs):
    """
DESCRIPTION

    Load two structures, superpose the second onto the first, report the
    RMSD and colour both by how far apart matching residues are.

    The method is chosen from the sequence identity of the two proteins:
    "align" (>= 30% identical), "super" (15-30%), "cealign" (below 15%).

USAGE

    compare <query1> <query2> [method=auto|align|super|cealign]

    A query is a protein name, a PDB ID or a loaded object.  Add :CHAIN to
    use one chain.  Separate multi-word names with a comma or "vs".

EXAMPLES

    compare 1HHO 4HHB
    compare 4HHB:A 1MBN
    compare human myoglobin, horse myoglobin
    """
    raw = base._raw(args, kwargs, " , ")
    parts, options = _split_compare(raw)
    label = "compare %s" % " ".join(base._raw(args, kwargs).split())
    method = options.get("method", "auto").lower()
    problem = ""
    if len(parts) != 2:
        problem = ("compare needs exactly two structures, e.g.  compare 1HHO 4HHB   or   "
                   "compare human myoglobin, horse myoglobin")
    elif method not in ("auto", "align", "super", "cealign"):
        problem = "method must be auto, align, super or cealign (you wrote '%s')." % method
    if problem:
        _in_order(lambda: (base._header(label), _say(problem), _say(report.LINE)))
        return
    sides = [_split_chain(p) for p in parts]
    loaded = {}

    def before():
        for name in base._molecules():
            loaded[name.lower()] = name
        for text, chain in sides:          # names used in an earlier grab count too
            known = base._known_object(text)
            if known:
                loaded[text.lower()] = known
        return True

    def job():
        out = []
        for text, chain in sides:
            if text.lower() in loaded:
                out.append({"object": loaded[text.lower()], "chain": chain})
                continue
            result = search.find(text + " limit=3")
            if result.kind != "results":
                return {"error": report.render(result, "find")}
            pdb_id = result.top.pdb_id
            try:
                path = base._download_pdb(pdb_id)
            except http.NetworkError as err:
                return {"error": base._network_message(err)}
            note = ""
            if result.mode != "id":
                note = "'%s' -> %s (%s)" % (text, pdb_id, result.top.entry.title[:50])
            out.append({"object": pdb_id, "chain": chain, "path": path, "note": note})
        return {"sides": out}

    def done(answer, error):
        if error:
            return base._unexpected(error)
        if "error" in answer:
            base._header(label)
            _say("Could not pick one structure for each side:")
            _say(answer["error"])
            return
        base._header(label)
        for side in answer["sides"]:
            if side.get("note"):
                _say("Using " + side["note"])
            if side.get("path") and side["object"] not in cmd.get_names("objects"):
                cmd.load(side["path"], side["object"])
                base._objects[side["object"]] = ("pdb", side["object"])
        first, second = answer["sides"]
        if first["object"] == second["object"] and first["chain"] == second["chain"]:
            _say("Both sides are the same structure (%s). To compare two chains of it, "
                 "write e.g.  compare %s:A %s:B" % ((first["object"],) * 3))
            _say(report.LINE)
            return
        base._snapshot(label)
        _superpose(first, second, method)
        _say(report.LINE)

    worker.submit(job, done, before)


def _residue_names(selection):
    names = {}
    cmd.iterate(selection, "names[(model, index)] = (resn, chain, resi)", space={"names": names})
    return names


def _pairs(alignment):
    """[(target atom, mobile atom)] as (object, index) tuples from an alignment object."""
    pairs = []
    for column in cmd.get_raw_alignment(alignment):
        if len(column) == 2:
            pairs.append((column[0], column[1]))
    return pairs


def _superpose(first, second, method):
    target_obj, mobile_obj = first["object"], second["object"]
    same_object = target_obj == mobile_obj
    target = _protein_ca(target_obj, first["chain"])
    mobile = _protein_ca(mobile_obj, second["chain"])
    n_target, n_mobile = cmd.count_atoms(target), cmd.count_atoms(mobile)
    if n_target < 8 or n_mobile < 8:
        which = target_obj if n_target < 8 else mobile_obj
        _say("%s has too few protein residues to superpose (is the chain ID right, or is "
             "it DNA/RNA?). Chains present: %s" % (which, ", ".join(cmd.get_chains(which))))
        return
    if same_object:
        _say("Both chains are in the same object, and superposing would move the whole "
             "object. Make a copy of one chain first:")
        _say("  create chain_copy, %s and chain %s" % (mobile_obj, second["chain"] or "B"))
        _say("  compare %s:%s chain_copy" % (target_obj, first["chain"] or "A"))
        return

    # 1. how similar are the sequences?  (alignment only - nothing moves yet)
    scratch = "_prospector_aln"
    cmd.delete(scratch)
    cmd.align(mobile, target, cycles=0, transform=0, object=scratch)
    names = _residue_names("(%s) or (%s)" % (target, mobile))
    pairs = _pairs(scratch)
    identical = sum(1 for a, b in pairs if names.get(a, ("?",))[0] == names.get(b, ("!",))[0])
    identity = identical / float(min(n_target, n_mobile))
    cmd.delete(scratch)

    chosen = method
    if method == "auto":
        if identity >= 0.30:
            chosen, why = "align", "sequences are similar, so a sequence alignment is reliable"
        elif identity >= 0.15:
            chosen, why = "super", ("sequences are only distantly related, so the fit uses "
                                    "the 3D structure instead of the sequence")
        else:
            chosen, why = "cealign", ("sequences are unrelated, so only the overall fold "
                                      "can be matched (CE algorithm)")
    else:
        why = "you asked for it"
    how = ["Sequence identity: %.0f%% (%d identical of %d aligned residues; chains of %d "
           "and %d residues)." % (100 * identity, identical, len(pairs), n_target, n_mobile),
           "Method: %s - %s." % (chosen, why)]

    # 2. superpose: the second structure moves onto the first
    rmsd, used = None, 0
    ce_pairs = None
    try:
        if chosen == "align":
            out = cmd.align(mobile, target)
            rmsd, used = out[0], out[1]
        elif chosen == "super":
            out = cmd.super(mobile, target)
            rmsd, used = out[0], out[1]
        else:
            cmd.delete(scratch)
            try:
                out = cmd.cealign(target, mobile, object=scratch)
                ce_pairs = _pairs(scratch)       # the residues CE actually matched
            except TypeError:                    # older PyMOL: no object= argument
                out = cmd.cealign(target, mobile)
            rmsd, used = out["RMSD"], out["alignment_length"]
    except Exception as err:
        if chosen == "cealign":
            _say("cealign could not match these structures (%s); using super instead."
                 % str(err).strip()[:80])
            chosen = "super"
            out = cmd.super(mobile, target)
            rmsd, used = out[0], out[1]
        else:
            _say("PyMOL could not superpose these structures: %s" % str(err).strip()[:120])
            return

    # 3. residue-by-residue distances after the fit
    cmd.delete(scratch)
    if chosen == "align":
        cmd.align(mobile, target, cycles=0, transform=0, object=scratch)
    else:
        cmd.super(mobile, target, cycles=0, transform=0, object=scratch)
    pairs = ce_pairs if (chosen == "cealign" and ce_pairs) else _pairs(scratch)
    cmd.delete(scratch)
    coords = {}
    cmd.iterate_state(1, "(%s) or (%s)" % (target, mobile),
                      "coords[(model, index)] = (x, y, z)", space={"coords": coords})
    distances = []
    for a, b in pairs:
        if a in coords and b in coords:
            d = sum((coords[a][i] - coords[b][i]) ** 2 for i in range(3)) ** 0.5
            distances.append((d, a, b))
    all_rmsd = (sum(d * d for d, _, _ in distances) / len(distances)) ** 0.5 if distances \
        else 0.0

    # 4. colour both structures by distance
    bins = [("< 1 A   (same)", 0.0, 1.0, "skyblue"), ("1 - 2 A", 1.0, 2.0, "paleyellow"),
            ("2 - 4 A", 2.0, 4.0, "orange"), ("> 4 A   (different)", 4.0, 1e9, "red")]
    # cartoon only (no waters / ligands), so the comparison is easy to read
    cmd.hide("everything", "%s or %s" % (target_obj, mobile_obj))
    for obj, chain in ((target_obj, first["chain"]), (mobile_obj, second["chain"])):
        # when one chain was asked for, the other chains would only be in the way
        cmd.show("cartoon", "%s and polymer%s" % (obj, " and chain %s" % chain if chain
                                                  else ""))
    cmd.color("grey80", "(%s or %s) and polymer" % (target_obj, mobile_obj))
    counts = []
    moved = {target_obj: [], mobile_obj: []}
    for title, low, high, colour in bins:
        members = {target_obj: [], mobile_obj: []}
        for d, a, b in distances:
            if low <= d < high:
                members[a[0]].append(a[1])
                members[b[0]].append(b[1])
                if d >= 2.0:
                    moved[a[0]].append(a[1])
                    moved[b[0]].append(b[1])
        counts.append((title, colour, len(members[target_obj])))
        for obj, indices in members.items():
            for start in range(0, len(indices), 400):
                chunk = "+".join(str(i) for i in indices[start:start + 400])
                cmd.color(colour, "byres (%s and index %s)" % (obj, chunk))
    diff_name = "compare_diff"
    pieces = []
    for obj, indices in moved.items():
        for start in range(0, len(indices), 400):
            pieces.append("(%s and index %s)" % (obj, "+".join(
                str(i) for i in indices[start:start + 400])))
    cmd.delete(diff_name)
    if pieces:
        cmd.select(diff_name, "byres (%s)" % " or ".join(pieces), enable=0)
    cmd.deselect()
    cmd.zoom("byres ((%s) or (%s))" % (target, mobile), 3)

    if rmsd < 1.0:
        verdict = "essentially the same conformation"
    elif rmsd < 2.5:
        verdict = "same fold with local differences"
    elif rmsd < 5.0:
        verdict = "related fold, substantial differences"
    else:
        verdict = "very different (or only partly comparable)"
    # 1. the answer
    _say("")
    _say("RESULT       %s and %s: %s." % (target_obj, mobile_obj, verdict))
    _say("             RMSD: %.2f A  (average distance between matching atoms;" % rmsd)
    _say("             smaller = more alike)")
    if used < 0.5 * min(n_target, n_mobile):
        _say("             Caution: fewer than half of the residues could be matched, so "
             "these two")
        _say("             proteins probably do not share a fold and the RMSD means little.")
    # 2. what is on screen
    _say("")
    _say("ON SCREEN    %s was moved onto %s." % (mobile_obj, target_obj))
    _say("             Colour = how far each residue is from its partner:")
    for title, colour, count in counts:
        _say("                %-11s %-20s %4d residues" % (colour, title, count))
    _say("                %-11s %s" % ("grey80", "no partner found"))
    _say("             Only the compared chains are shown, as cartoons;  show cartoon  "
         "brings the")
    _say("             rest back.")
    # 3. details
    worst = sorted(distances, reverse=True)[:5]
    if worst and worst[0][0] >= 2.0:
        _say("")
        _say("BIGGEST DIFFERENCES")
        for d, a, b in worst:
            if d < 2.0:
                break
            ra, rb = names.get(a, ("?", "?", "?")), names.get(b, ("?", "?", "?"))
            _say("   %s %s/%s%s  <->  %s %s/%s%s   %.1f A"
                 % (a[0], ra[1], ra[0], ra[2], b[0], rb[1], rb[0], rb[2], d))
        if pieces:
            _say("   Every residue that differs by 2 A or more is in the selection  %s"
                 % diff_name)
    _say("")
    _say("HOW IT WAS DONE")
    for line in how:
        _say(textwrap.fill(line, 78, initial_indent="   ", subsequent_indent="   "))
    if chosen == "cealign":
        _say("   The RMSD covers the %d residues the CE algorithm could match." % used)
    else:
        _say("   The RMSD covers %d CA atoms, after PyMOL set aside poorly fitting pairs;"
             % used)
        _say("   over all %d matched pairs it is %.2f A." % (len(distances), all_rmsd))
    if chosen == "cealign":
        teach.note("cealign %s, %s" % (target, mobile), "target first, then the one that moves")
    else:
        teach.note("%s %s, %s" % (chosen, mobile, target),
                   "the first selection moves onto the second")
    teach.note("hide everything, %s or %s" % (target_obj, mobile_obj))
    teach.note("show cartoon, (%s or %s) and polymer" % (target_obj, mobile_obj))
    teach.note("color red, %s" % diff_name if pieces else "color grey80, polymer",
               "Prospector colours by the CA-CA distance it measured")


# ------------------------------------------------------------------ similar
def similar(*args, **kwargs):
    """
DESCRIPTION

    Find structures with a similar sequence (search run by RCSB).  Results
    are grouped so each row is a different protein.

USAGE

    similar [object] [chain=A] [identity=30] [limit=10]

EXAMPLES

    similar 1MBN
    similar 4HHB chain=B identity=20 limit=15
    similar myoglobin           (finds and loads the best match first)
    """
    words, options = base._split_options(base._raw(args, kwargs))
    label = "similar %s" % base._raw(args, kwargs)
    state = {}

    def before():
        typed = " ".join(words)
        obj, message = base._pick_object(typed, "similar")
        if obj is None and typed:
            state["query"] = typed          # a protein name: find and load it first
            obj = "pending"
        try:
            identity = float(options.get("identity", 30))
            limit = int(options.get("limit", 10))
            if not (10 <= identity <= 100 and 1 <= limit <= 50):
                raise ValueError
        except ValueError:
            obj, message = None, ("identity must be 10-100 (percent) and limit 1-50, e.g. "
                                  " similar 1MBN identity=30 limit=10")
        if obj is None:
            base._header(label)
            _say(message)
            _say(report.LINE)
            return False
        if state.get("query"):
            state.update(obj=None, chain=options.get("chain", ""), identity=identity,
                         limit=limit, identifier=None, sequence="")
            return True
        chains = [c for c in cmd.get_chains(obj)
                  if cmd.count_atoms(_protein_ca(obj, c)) >= 1]
        chain = options.get("chain", chains[0] if chains else "")
        residues = []
        cmd.iterate(_protein_ca(obj, chain), "residues.append(resn)",
                    space={"residues": residues})
        if not residues:
            base._header(label)
            _say("%s has no protein chain '%s'. Chains with protein: %s"
                 % (obj, chain, ", ".join(chains) or "none"))
            _say(report.LINE)
            return False
        state.update(obj=obj, chain=chain, identity=identity, limit=limit,
                     identifier=base._identify(obj),
                     sequence="".join(feat.THREE_TO_ONE.get(r.upper(), "X") for r in residues))
        return True

    def job():
        try:
            if state.get("query"):
                fetched = state["fetched"] = base._fetch_query(state["query"])
                if fetched["error"]:
                    return None, fetched["error"]
                state["identifier"] = fetched["identity"]
                if fetched["identity"][0] == "pdb":
                    sequences = rcsb.chain_sequences(fetched["identity"][1])
                    if not state["chain"] and sequences:
                        state["chain"] = sorted(sequences)[0]
                    state["sequence"] = sequences.get(state["chain"], "")
                else:
                    from ..core import uniprot
                    record = uniprot.details(fetched["identity"][1])[1] or {}
                    state["sequence"] = (record.get("sequence") or {}).get("value", "")
                    state["chain"] = "A"
            sequence = state["sequence"]
            identifier = state["identifier"]
            exclude = None
            if identifier and identifier[0] == "pdb":
                exclude = identifier[1]
                # the full deposited sequence is better than the residues that
                # happen to be visible in the model
                full = rcsb.chain_sequences(identifier[1]).get(state["chain"])
                if full:
                    sequence = full
            return sim.find_similar(sequence, state["identity"], state["limit"], exclude), None
        except http.NetworkError as err:
            return None, base._network_message(err)

    def done(answer, error):
        if error:
            return base._unexpected(error)
        result, trouble = answer
        base._header(label)
        if state.get("fetched") and not state["fetched"]["error"]:
            state["obj"] = base._place_fetched(state["fetched"])
        if state.get("obj") and base._needs_protein(state["obj"], "similar"):
            return
        if trouble:
            _say(trouble)
        else:
            _say(sim.render(result, "%s chain %s" % (state["obj"], state["chain"])))
            others = [f.best[2]["pdb_id"] for f in result.families
                      if f.best[2]["pdb_id"] != (state["identifier"] or ("", ""))[1]]
            if others:
                _say("")
                _say("Try:  compare %s %s" % (state["obj"], others[min(1, len(others) - 1)]))
        _say(report.LINE)

    worker.submit(job, done, before)


# -------------------------------------------------------------------- batch
def batch(*args, **kwargs):
    """
DESCRIPTION

    Load many structures at once from a text file: one protein name or PDB
    ID per line (filters allowed; lines starting with # are ignored).
    Ambiguous names are skipped and reported, never guessed.

USAGE

    batch <file.txt>

EXAMPLE FILE

    # my proteins
    4HHB
    myoglobin organism=horse
    p53 res<1.5
    """
    raw = base._raw(args, kwargs).strip().strip("\"'")
    label = "batch %s" % raw
    path = os.path.expanduser(raw)
    problem, lines = "", []
    if not raw:
        problem = "Give a text file with one name or PDB ID per line, e.g.  batch proteins.txt"
    elif not os.path.isfile(path):
        problem = ("I cannot find the file '%s'. PyMOL's current folder is %s - use  cd  "
                   "or give the full path." % (raw, os.getcwd()))
    else:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        lines.append(line)
        except OSError as err:
            problem = "Could not read '%s': %s" % (raw, err)
        if not problem and not lines:
            problem = "'%s' contains no names (empty lines and # comments are ignored)." % raw
        if len(lines) > 100:
            problem = "'%s' has %d entries; batch handles at most 100 at a time." % (
                raw, len(lines))
    if problem:
        _in_order(lambda: (base._header(label), _say(problem), _say(report.LINE)))
        return

    outcome = []      # (line, status text)
    _in_order(lambda: (base._header(label),
                       _say("Loading %d structures - PyMOL stays usable meanwhile."
                            % len(lines))))

    def make(line):
        def job():
            result = search.find(line)
            item, problem_text = None, ""
            try:
                if result.kind == "results":
                    item = base._download_pdb(result.top.pdb_id)
                elif result.kind == "alphafold":
                    item = base._download_model(result.model)
            except http.NetworkError as err:
                problem_text = str(err)
            return result, item, problem_text

        def done(answer, error):
            if error:
                outcome.append((line, "internal error"))
                return base._unexpected(error)
            result, item, problem_text = answer
            base._notify("batch", result)
            if result.kind == "results" and item:
                name = result.top.pdb_id
                if name not in cmd.get_names("objects"):
                    cmd.load(item, name)
                base._objects[name] = ("pdb", name)
                base._remember(result.query.text, name)
                outcome.append((line, "loaded " + name))
                _say("  %-28s -> %s" % (line[:28], report.one_line(result.top)[:110]))
            elif result.kind == "alphafold" and item:
                name = result.model.object_name()
                if name not in cmd.get_names("objects"):
                    cmd.load(item, name)
                base._objects[name] = ("af", result.model.accession)
                outcome.append((line, "loaded AlphaFold model " + name))
                _say("  %-28s -> no experimental structure; AlphaFold model %s"
                     % (line[:28], name))
            else:
                reasons = {"ambiguous": "ambiguous name - run  find %s  and choose" % line,
                           "suggest": "not found - did you mean '%s'?" % result.suggestion,
                           "none": "nothing found", "nofilter": "no structure passes the "
                           "filters", "error": result.message}
                reason = problem_text or reasons.get(result.kind, result.kind)
                outcome.append((line, "SKIPPED: " + reason))
                _say("  %-28s -> SKIPPED: %s" % (line[:28], reason[:100]))
        return job, done

    for line in lines:
        job, done = make(line)
        worker.submit(job, done)

    def summary():
        ok = sum(1 for _, status in outcome if status.startswith("loaded"))
        _say("Batch finished: %d of %d loaded%s." % (
            ok, len(lines), "" if ok == len(lines) else
            "; %d skipped (see the lines marked SKIPPED)" % (len(lines) - ok)))
        if ok:
            cmd.zoom()
        _say(report.LINE)
    _in_order(summary)


# --------------------------------------------------------------------- look
_HYDRO = [   # Kyte-Doolittle hydropathy, grouped into five classes
    ("prospector_hp5", [0.90, 0.45, 0.10], "ILE+VAL+LEU+PHE", "very hydrophobic"),
    ("prospector_hp4", [0.96, 0.72, 0.45], "CYS+MET+ALA", "hydrophobic"),
    ("prospector_hp3", [0.95, 0.95, 0.95], "GLY+THR+SER+TRP+TYR+PRO", "neutral"),
    ("prospector_hp2", [0.65, 0.80, 0.95], "HIS+GLU+GLN+ASP+ASN", "polar"),
    ("prospector_hp1", [0.30, 0.55, 0.90], "LYS+ARG", "very polar (charged)"),
]

PRESETS = {
    "publication": ("white background; protein as ribbons, bound molecules as sticks, "
                    "metal ions as small balls - ready for a figure", [
        "bg_color white",
        "hide everything, {o}",
        "show cartoon, {o} and polymer",
        "show sticks, {o} and organic",
        "show nb_spheres, {o} and inorganic",
        "set cartoon_fancy_helices, 1",
        "set cartoon_transparency, 0",
        "set ray_trace_mode, 1",
        "set ray_opaque_background, 0",
        "set antialias, 2",
        "set depth_cue, 0",
        "set specular, 0.2",
        "orient {o}",
    ]),
    "chains": ("protein as ribbons, every chain in its own colour; bound molecules as "
               "sticks", [
        "hide everything, {o}",
        "show cartoon, {o} and polymer",
        "util.cbc {o}",
        "show sticks, {o} and organic",
    ]),
    "rainbow": ("ribbons coloured along each chain: blue = start (N-terminus), through "
                "green and yellow, to red = end (C-terminus)", [
        "hide everything, {o}",
        "show cartoon, {o} and polymer",
        "util.chainbow {o}",
    ]),
    "hydrophobicity": ("the protein's outer surface: orange = water-avoiding "
                       "(hydrophobic) patches, blue = charged, white = in between", None),
    "putty": ("a tube instead of a ribbon: thick and red = parts that move a lot, thin "
              "and blue = rigid parts (from the B-factors)", [
        "hide everything, {o}",
        "show cartoon, {o} and polymer",
        "cartoon putty, {o}",
        "set cartoon_putty_transform, 0",
        "spectrum b, blue_white_red, {o} and polymer",
    ]),
    "secondary": ("ribbons coloured by shape: red = helices, yellow = strands "
                  "(sheets), green = loops", [
        "hide everything, {o}",
        "show cartoon, {o} and polymer",
        "color green, {o} and polymer",
        "color red, {o} and ss h",
        "color yellow, {o} and ss s",
        "show sticks, {o} and organic",
    ]),
    "elements": ("every atom as sticks, coloured by element: green = carbon, blue = "
                 "nitrogen, red = oxygen, yellow = sulfur", [
        "hide everything, {o}",
        "show sticks, {o} and not solvent",
        "show nb_spheres, {o} and inorganic",
        "util.cbag {o}",
    ]),
    "types": ("ribbons coloured by the kind of amino acid: orange = water-avoiding "
              "(hydrophobic), green = polar, blue = positive, red = negative", [
        "hide everything, {o}",
        "show cartoon, {o} and polymer",
        "color green, {o} and polymer",
        "color orange, {o} and resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO",
        "color blue, {o} and resn LYS+ARG",
        "color red, {o} and resn ASP+GLU",
        "show sticks, {o} and organic",
    ]),
    "ligand": ("a see-through grey surface over the ribbon, with bound molecules as "
               "sticks inside it - shows how buried they are", [
        "hide everything, {o}",
        "show cartoon, {o} and polymer",
        "show surface, {o} and polymer",
        "color grey80, {o} and polymer",
        "set transparency, 0.3",
        "show sticks, {o} and organic",
        "util.cbag {o} and organic",
        "show nb_spheres, {o} and inorganic",
        "zoom {o} and organic, 10",
    ]),
}
# colour schemes that need not change how the atoms are drawn: with these,
# a drawing picked before (viewas spheres) is kept and only recoloured
COLOUR_ONLY = {
    "chains": "every chain in its own colour",
    "rainbow": "coloured along each chain: blue = start (N-terminus), through green and "
               "yellow, to red = end (C-terminus)",
    "secondary": "coloured by shape: red = helices, yellow = strands (sheets), "
                 "green = loops",
    "types": "coloured by the kind of amino acid: orange = water-avoiding (hydrophobic), "
             "green = polar, blue = positive, red = negative",
    "elements": "coloured by element: green = carbon, blue = nitrogen, red = oxygen, "
                "yellow = sulfur",
}


def _drawn_as(obj):
    """The drawing chosen for the protein when it is not the usual ribbon
    ('spheres', 'surface' ...), else None."""
    try:
        total = cmd.count_atoms("(%s) and polymer" % obj)
        if not total:
            return None
        for rep in ("surface", "spheres", "sticks", "lines", "mesh", "dots", "ribbon"):
            if cmd.count_atoms("(%s) and polymer and rep %s" % (obj, rep)) > total * 0.5:
                return rep
    except Exception:
        pass
    return None


_LOOK_COMMENTS = {     # "command first-argument" -> short explanation
    "bg_color white": "white background",
    "hide everything": "start from a blank slate",
    "show cartoon": "ribbons for the protein backbone",
    "show sticks": "ligands as sticks (organic = small molecules)",
    "show nb_spheres": "metal ions as small spheres",
    "show surface": "molecular surface",
    "set cartoon_fancy_helices": "helices with thin edges",
    "set ray_trace_mode": "1 = colour plus black outlines",
    "set ray_opaque_background": "0 = transparent background in the PNG",
    "set antialias": "smoother edges",
    "set depth_cue": "0 = no fog",
    "set transparency": "see-through surface",
    "set cartoon_putty_transform": "0 = thickness follows the raw B-factor",
    "orient all": "turn the molecule to fill the window",
    "util.cbc": "colour by chain",
    "util.chainbow": "rainbow along each chain, N-terminus blue",
    "util.cbag": "colour by atom, carbons green",
    "cartoon putty": "tube whose thickness follows the B-factor",
    "spectrum b": "colour by B-factor: blue low, red high",
    "color grey80": "light grey protein",
}

def _look_comment(line):
    words = line.replace(",", " ").split()
    two = " ".join(words[:2])
    return _LOOK_COMMENTS.get(two, _LOOK_COMMENTS.get(words[0], "") if words else "")


_PRESET_ALIASES = {"pub": "publication", "chain": "chains", "cartoon": "chains",
                   "hydrophobic": "hydrophobicity", "hydro": "hydrophobicity",
                   "bfactor": "putty", "b-factor": "putty", "surface": "ligand",
                   "ligands": "ligand", "ss": "secondary", "element": "elements",
                   "atoms": "elements", "type": "types", "properties": "types"}


def look(*args, **kwargs):
    """
DESCRIPTION

    Apply a styling preset and print the PyMOL commands it used.

USAGE

    look <publication|chains|rainbow|secondary|types|elements|hydrophobicity|
          putty|ligand> [object]

EXAMPLES

    look publication
    look putty 4HHB
    look ligand haemoglobin     (finds and loads the best match first)
    """
    words, _ = base._split_options(base._raw(args, kwargs))
    label = "look %s" % " ".join(words)
    # the preset may come first or last:  look rainbow p53  =  look p53 rainbow
    known = [i for i, w in enumerate(words)
             if _PRESET_ALIASES.get(w.lower(), w.lower()) in PRESETS]
    if known and known[0] != 0:
        words = [words[known[0]]] + words[:known[0]] + words[known[0] + 1:]
    elif not known and len(words) >= 1:
        import difflib                   # look ranbow -> did you mean rainbow?
        for index, word in enumerate(words):
            close = difflib.get_close_matches(word.lower(), sorted(PRESETS), n=1, cutoff=0.75)
            if close and not base._known_object(word):
                words = ["?%s?%s" % (word, close[0])] + words[:index] + words[index + 1:]
                break

    target = " ".join(words[1:])
    state = {}

    def before():
        state["need"] = bool(target) and base._known_object(target) is None \
            and target.lower() != "all" and not words[0].startswith("?")
        return True

    def job():
        if state["need"]:        # a protein name: find and load it first
            state["fetched"] = base._fetch_query(target, strict=True)

    def action():
        base._header(label)
        fetched = state.get("fetched")
        if fetched and fetched["error"]:
            _say(fetched["error"])
            _say(report.LINE)
            return
        if fetched:
            base._place_fetched(fetched)
        preset = _PRESET_ALIASES.get(words[0].lower(), words[0].lower()) if words else ""
        if preset.startswith("?"):                     # a near miss: say which one
            _, typed_word, meant = preset.split("?")
            _say("There is no colour scheme called '%s'. Did you mean '%s'?   look %s"
                 % (typed_word, meant, meant))
            _say(report.LINE)
            return
        if preset not in PRESETS:
            if preset:
                _say("There is no preset called '%s'." % words[0])
            _say("Presets:")
            for name in sorted(PRESETS):
                _say("  look %s" % name)
                _say("       %s" % PRESETS[name][0])
            _say(report.LINE)
            return
        obj, message = base._pick_object(target, "look " + preset)
        if target.lower() == "all":
            obj = "all"
        if obj is None and len(words) < 2 and base._molecules():
            obj = "all"          # no object named: style everything
        if obj is None:
            _say(message)
            _say(report.LINE)
            return
        description, lines = PRESETS[preset]
        if preset == "hydrophobicity":
            lines = ["hide everything, {o}", "show surface, {o} and polymer"]
            for name, rgb, residues, _ in _HYDRO:
                lines.append("set_color %s, [%.2f, %.2f, %.2f]" % (name, rgb[0], rgb[1], rgb[2]))
                lines.append("color %s, {o} and polymer and resn %s" % (name, residues))
            lines.append("show sticks, {o} and organic")
        if preset == "putty":
            if obj != "all" and (base._identify(obj) or ("",))[0] == "af":
                _say("Note: in an AlphaFold model the B-factor column holds the confidence "
                     "(pLDDT), so here thick/red means CONFIDENT, not flexible.")
        kept = None
        if preset in COLOUR_ONLY:
            kept = _drawn_as(obj)
            if kept:            # spheres/sticks/surface chosen before: only recolour
                lines = [line for line in lines if not line.startswith(("hide ", "show "))]
                description = COLOUR_ONLY[preset] + (
                    " - still drawn as %s (only the colours changed)" % kept)
        base._snapshot(label)
        for name in (base._molecules() if obj == "all" else [obj]):
            cmd.set("cartoon_transparency", 0, name)   # undo a `focus` close-up
            cmd.label(name, "")
        ran = [line.replace("{o}", obj) for line in lines]
        for line in ran:
            cmd.do(line, echo=0)
        _say("")
        _say("Applied '%s' to %s." % (preset, obj if obj != "all" else "everything"))
        _say("ON SCREEN    %s" % description)
        if preset == "hydrophobicity":
            _say("COLOUR KEY   (Kyte-Doolittle scale, in five steps from orange to blue)")
            for name, _, residues, meaning in _HYDRO:
                _say("                %-22s %s" % (meaning, residues.replace("+", " ")))
        if preset == "publication":
            _say("For the final image run:  snap my_figure")
        for line in ran:
            teach.note(line, _look_comment(line))
        _say(report.LINE)

    worker.submit(job, lambda result, error: action() if not error
                  else base._unexpected(error), before)


# ------------------------------------------------------------------- viewas
# name -> (one-line description for the menu, what you then see, PyMOL commands)
STYLES = [
    ("cartoon", "ribbon that follows the backbone (the usual overview)",
     "the protein as a ribbon: spirals = alpha-helices, flat arrows = beta-strands, "
     "thin tubes = loops; bound molecules stay as sticks",
     ["show_as cartoon, {p}", "cartoon automatic, {o}", "show sticks, {l}"]),
    ("lines", "thin lines along every bond - all atoms, very light",
     "every bond as a thin line, so you see each amino acid's full shape; "
     "line ends and corners are atoms",
     ["show_as lines, {a}"]),
    ("sticks", "thicker bonds - all atoms, easy to follow",
     "every bond as a stick; corners and ends are atoms (blue = nitrogen, red = oxygen, "
     "yellow = sulfur)",
     ["show_as sticks, {a}"]),
    ("ballstick", "balls for atoms joined by thin sticks",
     "each atom as a small ball and each bond as a thin stick - the classic chemistry "
     "model",
     ["show_as sticks, {a}", "show spheres, {a}", "set stick_radius, 0.14, {o}",
      "set sphere_scale, 0.25, {o}"]),
    ("spheres", "every atom as a full-size ball (space-filling)",
     "each atom at its real size, so you see how tightly the protein is packed and "
     "what its outside looks like",
     ["show_as spheres, {a}"]),
    ("ribbon", "one thin line through the backbone",
     "a single thin line tracing the chain from residue to residue - the simplest view "
     "of its path; bound molecules stay as sticks",
     ["show_as ribbon, {p}", "show sticks, {l}"]),
    ("tube", "smooth tube along the backbone",
     "the backbone as an even, smooth tube, without helix and strand shapes; bound "
     "molecules stay as sticks",
     ["show_as cartoon, {p}", "cartoon tube, {o}", "show sticks, {l}"]),
    ("surface", "solid outer surface",
     "the solid outer skin of the molecule - what another molecule would touch; grooves "
     "and pockets show up as dents",
     ["show_as surface, {a}"]),
    ("mesh", "surface drawn as a wire net",
     "the outer surface as a net you can see through",
     ["show_as mesh, {a}"]),
    ("dots", "surface drawn as dots",
     "the outer surface as a cloud of dots you can see through",
     ["show_as dots, {a}"]),
]
STYLE_NAMES = [item[0] for item in STYLES]
_STYLE_ALIASES = {
    "line": "lines", "wire": "lines", "wireframe": "lines", "stick": "sticks",
    "licorice": "sticks", "ballandstick": "ballstick", "ball+stick": "ballstick",
    "ballsticks": "ballstick", "ballsandsticks": "ballstick", "ball": "spheres",
    "balls": "spheres", "sphere": "spheres", "cpk": "spheres", "spacefill": "spheres",
    "spacefilling": "spheres", "ribbons": "cartoon", "cartoons": "cartoon",
    "trace": "ribbon", "backbone": "ribbon", "tubes": "tube", "surfaces": "surface",
    "dot": "dots", "meshes": "mesh",
}


def _style_menu():
    _say("")
    _say("WAYS TO VIEW A STRUCTURE")
    for name, short, _, _ in STYLES:
        _say("   %-10s %s" % (name, short))
    _say("")
    _say("HOW TO USE   viewas <style>             the loaded structure")
    _say("             viewas <protein> <style>   e.g.  viewas p53 sticks")
    _say("             viewas <style>+<style>     combine, e.g.  viewas cartoon+sticks")
    _say("Colours do not change. For ready-made colour schemes see:  look")


def viewas(*args, **kwargs):
    """
DESCRIPTION

    Change how a structure is drawn: lines, sticks, balls, ribbon, surface...
    Type "viewas" alone to see every option.

USAGE

    viewas                          list the options
    viewas <style> [protein]
    viewas <protein> <style>
    viewas <style>+<style>          combine two or more

EXAMPLES

    viewas sticks
    viewas p53 ballstick
    viewas cartoon+sticks
    viewas haemoglobin surface     (finds and loads the best match first)
    """
    raw = " ".join(base._raw(args, kwargs).replace(",", " ").split())
    lowered = re.sub(r"ball(s)?\s*(and|&|\+)\s*stick(s)?", "ballstick", raw.lower())
    chosen, rest = [], []
    for word in lowered.split():
        parts = [_STYLE_ALIASES.get(p, p) for p in word.split("+") if p]
        if parts and all(p in STYLE_NAMES for p in parts):
            chosen.extend(p for p in parts if p not in chosen)
        else:
            rest.append(word)
    target = " ".join(rest)
    label = "viewas %s" % raw
    state = {}

    def before():
        state["need"] = bool(target) and bool(chosen) and \
            base._known_object(target) is None and target != "all"
        return True

    def job():
        if state["need"]:        # a protein name: find and load it first
            state["fetched"] = base._fetch_query(target, strict=True)

    def action():
        base._header(label)
        if not chosen:
            if rest and base._known_object(target) is None:
                import difflib
                close = difflib.get_close_matches(rest[-1], STYLE_NAMES
                                                  + list(_STYLE_ALIASES), n=1, cutoff=0.7)
                _say("I don't know the style '%s'.%s" % (
                    rest[-1], " Did you mean '%s'?" % _STYLE_ALIASES.get(close[0], close[0])
                    if close else ""))
            _style_menu()
            _say(report.LINE)
            return
        fetched = state.get("fetched")
        if fetched and fetched["error"]:
            _say(fetched["error"])
            _say(report.LINE)
            return
        if fetched:
            base._place_fetched(fetched)
        obj, message = base._pick_object(target, "viewas " + chosen[0])
        if obj is None:
            _say(message)
            _say(report.LINE)
            return
        base._snapshot(label)
        # undo things earlier commands may have left on this object
        cmd.set("cartoon_transparency", 0, obj)
        cmd.set("stick_radius", 0.25, obj)
        cmd.set("sphere_scale", 1.0, obj)
        cmd.cartoon("automatic", obj)
        cmd.label(obj, "")
        cmd.enable(obj)
        cmd.hide("everything", "%s and solvent" % obj)
        teach.note("hide everything, %s and solvent" % obj, "water would only clutter")
        fill = {"o": obj, "a": "%s and not solvent" % obj,
                "p": "%s and polymer" % obj,
                "l": "%s and not polymer and not solvent" % obj}
        table = dict((item[0], item) for item in STYLES)
        for index, name in enumerate(chosen):
            for line in table[name][3]:
                line = line.replace("{o}", fill["o"]).replace("{a}", fill["a"]) \
                           .replace("{p}", fill["p"]).replace("{l}", fill["l"])
                if index > 0:          # later styles are added on top of the first
                    line = line.replace("show_as ", "show ")
                cmd.do(line, echo=0)
                teach.note(line, "show_as = show only this" if line.startswith("show_as")
                           else "")
        cmd.deselect()
        _say("")
        _say("NOW VIEWING  %s as %s" % (obj, " + ".join(
            "ball and stick" if name == "ballstick" else name for name in chosen)))
        for index, name in enumerate(chosen):
            _say("%s%s" % ("ON SCREEN    " if index == 0 else "PLUS         ",
                           table[name][2]))
        _say("             (water molecules are hidden)")
        others = [n for n in STYLE_NAMES if n not in chosen]
        _say("TRY ANOTHER  viewas %s  |  viewas %s  |  viewas  (all options)"
             % (others[0], others[1]))
        _say(report.LINE)

    worker.submit(job, lambda result, error: action() if not error
                  else base._unexpected(error), before)



# --------------------------------------------------------------------- snap
def snap(*args, **kwargs):
    """
DESCRIPTION

    Save a ray-traced, high-resolution PNG of the current view.  An existing
    file is never overwritten: a number is added to the name instead.

USAGE

    snap <name> [width=2400] [dpi=300]

EXAMPLES

    snap hemoglobin_figure
    snap ~/Desktop/p53 width=3600
    """
    words, options = base._split_options(base._raw(args, kwargs))
    label = "snap %s" % base._raw(args, kwargs)

    def action():
        base._header(label)
        try:
            width = int(options.get("width", 2400))
            dpi = int(options.get("dpi", 300))
            if not (200 <= width <= 8000 and 50 <= dpi <= 1200):
                raise ValueError
        except ValueError:
            _say("width must be 200-8000 pixels and dpi 50-1200, e.g.  snap figure "
                 "width=2400 dpi=300")
            _say(report.LINE)
            return
        if not words:
            _say("Give the picture a name, e.g.  snap hemoglobin_figure")
            _say(report.LINE)
            return
        if not base._molecules():
            _say("Nothing is loaded, so there is nothing to draw. Try  grab hemoglobin  first.")
            _say(report.LINE)
            return
        name = os.path.expanduser(" ".join(words))
        if not name.lower().endswith(".png"):
            name += ".png"
        folder = os.path.dirname(name)
        if not folder:
            folder = os.getcwd()
            if not os.access(folder, os.W_OK):
                folder = os.path.join(cache.home(), "snaps")
            name = os.path.join(folder, name)
        try:
            os.makedirs(os.path.dirname(name), exist_ok=True)
        except OSError as err:
            _say("Cannot create the folder for '%s': %s" % (name, err))
            _say(report.LINE)
            return
        stem, number, path = name[:-4], 2, name
        while os.path.exists(path):
            path = "%s_%d.png" % (stem, number)
            number += 1
        view_w, view_h = cmd.get_viewport()
        height = max(1, int(round(width * float(view_h) / max(1, view_w))))
        _say("Ray tracing %d x %d pixels - PyMOL is busy for a few seconds..." % (width, height))
        try:
            cmd.png(path, width=width, height=height, dpi=dpi, ray=1, quiet=1)
        except Exception as err:
            _say("PyMOL could not render the image: %s" % str(err).strip()[:150])
            _say(report.LINE)
            return
        if os.path.exists(path):
            teach.note("png %s, width=%d, height=%d, dpi=%d, ray=1"
                       % (path, width, height, dpi), "ray=1 = ray-traced")
            _say("Saved %s (%d x %d px, %d dpi = %.1f x %.1f cm in print, %.1f MB)."
                 % (path, width, height, dpi, width / float(dpi) * 2.54,
                    height / float(dpi) * 2.54, os.path.getsize(path) / 1e6))
            _say("PyMOL equivalent:  png %s, width=%d, height=%d, dpi=%d, ray=1"
                 % (path, width, height, dpi))
        else:
            _say("PyMOL did not write the file %s (is the folder writable?)." % path)
        _say(report.LINE)

    _in_order(action)


# ----------------------------------------------------------- history / fav
def history(*args, **kwargs):
    """
DESCRIPTION

    Show your past Prospector searches (kept between PyMOL sessions).

USAGE

    history [number | clear]

EXAMPLES

    history          (last 15)
    history 50
    history clear
    """
    words = base._raw(args, kwargs).split()

    def action():
        base._header("history %s" % " ".join(words))
        if words and words[0].lower() == "clear":
            _say("Removed %d history entries." % store.clear_history())
            _say(report.LINE)
            return
        count = 15
        if words:
            try:
                count = max(1, min(200, int(words[0])))
            except ValueError:
                _say("Use  history,  history 30  or  history clear.")
                _say(report.LINE)
                return
        items = store.history()
        if not items:
            _say("No searches yet. Try  find hemoglobin")
        shown = items[-count:]
        for item in shown:
            top = (" -> " + item["top"]) if item.get("top") else " -> (%s)" % item.get("outcome")
            _say("  %s   %s %s%s" % (store.when(item.get("t")), item.get("command", "find"),
                                     item.get("query", ""), top))
        if len(items) > len(shown):
            _say("  (%d older entries; use  history %d)" % (len(items) - len(shown),
                                                           len(items)))
        _say(report.LINE)

    _in_order(action)


def fav(*args, **kwargs):
    """
DESCRIPTION

    Keep a list of favourite structures (kept between PyMOL sessions).

USAGE

    fav                      list favourites
    fav add [ID] [note]      add a PDB ID (or the loaded structure)
    fav rm <ID>              remove one

EXAMPLES

    fav add 4HHB classic deoxy haemoglobin
    fav add
    fav rm 4HHB
    """
    words = base._raw(args, kwargs).split()

    def action():
        base._header("fav %s" % " ".join(words))
        verb = words[0].lower() if words else "list"
        if verb in ("add", "+"):
            rest = words[1:]
            identifier = ""
            if rest and (rcsb.looks_like_pdb_id(rest[0]) or rest[0].upper().startswith("AF")):
                identifier = rest[0].upper()
                rest = rest[1:]
            elif rest and any(base._known_object(" ".join(rest[:n]))
                              for n in range(len(rest), 0, -1)):
                # a name used in an earlier grab, e.g.  fav add haemoglobin my note
                size = [n for n in range(len(rest), 0, -1)
                        if base._known_object(" ".join(rest[:n]))][0]
                obj = base._known_object(" ".join(rest[:size]))
                identity = base._identify(obj)
                identifier = identity[1] if identity and identity[0] == "pdb" else obj
                rest = rest[size:]
            elif found.get("id"):
                # a protein name that is not loaded: its best match was looked up
                identifier = found["id"]
                _titles[identifier] = found["title"]
                _say("'%s' -> best match %s (%s)" % (rest[0], identifier,
                                                     found["title"][:44]))
                rest = rest[1:]
            elif found.get("error"):
                _say(found["error"])
                return
            else:
                obj, message = base._pick_object("", "fav add")
                identity = base._identify(obj) if obj else None
                if identity is None:
                    _say("Say which structure, e.g.  fav add 4HHB my note   (%s)"
                         % (message or "the loaded object has no PDB ID"))
                    _say(report.LINE)
                    return
                identifier = identity[1] if identity[0] == "pdb" else obj
            note = " ".join(rest)
            added = store.add_favorite(identifier, _titles.get(identifier, ""), note)
            _say("%s %s." % ("Added" if added else "Already a favourite - updated",
                             identifier))
        elif verb in ("rm", "remove", "del", "delete", "-"):
            if len(words) < 2:
                _say("Say which one, e.g.  fav rm 4HHB")
            elif store.remove_favorite(words[1]):
                _say("Removed %s from your favourites." % words[1].upper())
            else:
                _say("%s is not in your favourites." % words[1].upper())
        elif verb in ("list", "ls"):
            items = store.favorites()
            if not items:
                _say("No favourites yet. Add one with  fav add 4HHB  (or just  fav add  "
                     "for the loaded structure).")
            for item in items:
                text = item.get("note") or item.get("title") or ""
                _say("  %-14s %s  %-46s grab %s" % (item.get("id"),
                                                    store.when(item.get("t"))[:10],
                                                    text[:46], item.get("id")))
        else:
            _say("Use  fav,  fav add [ID] [note]  or  fav rm <ID>.")
        _say(report.LINE)

    found = {}

    def before():
        rest = words[1:]
        adding = bool(words) and words[0].lower() in ("add", "+")
        if adding and rest and not rcsb.looks_like_pdb_id(rest[0]) \
                and not rest[0].upper().startswith("AF") \
                and not any(base._known_object(" ".join(rest[:n]))
                            for n in range(len(rest), 0, -1)):
            found["query"] = rest[0]
        return True

    def job():
        if found.get("query"):
            result = search.find(found["query"] + " limit=1")
            if result.kind == "results":
                found["id"] = result.top.pdb_id
                found["title"] = result.top.entry.title
            else:
                found["error"] = ("'%s' is not loaded, so I searched for it:\n%s"
                                  % (found["query"], report.render(result, "find")))

    worker.submit(job, lambda result, error: action() if not error
                  else base._unexpected(error), before)


def register():
    table = {"compare": compare, "similar": similar, "batch": batch, "look": look,
             "viewas": viewas,
             "snap": snap, "history": history, "fav": fav}
    for name, function in table.items():
        cmd.extend(name, function)
        cmd.extend("p_" + name, function)
    try:
        cmd.auto_arg[0]["look"] = [lambda: cmd.Shortcut(sorted(PRESETS)), "preset", " "]
        cmd.auto_arg[0]["viewas"] = [lambda: cmd.Shortcut(STYLE_NAMES), "style", " "]
        cmd.auto_arg[0]["similar"] = cmd.auto_arg[0]["zoom"]
    except Exception:
        pass
    if _on_result not in base._listeners:
        base._listeners.insert(0, _on_result)
