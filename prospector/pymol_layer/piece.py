"""piece: copy part of a structure (a chain, or residues of a chain) into an
object of its own, so it appears in PyMOL's object list on the right and can
be coloured, hidden or shown by itself.

In plain PyMOL the same needs two steps and its selection language:

    select /4CHA//B/40-60          (only makes a selection, pink squares)
    create 4CHA_B40-60, 4CHA and chain B and resi 40-60

People write the part many ways; all of these are read:

    piece A 19-20                  piece chymotrypsinogen A/19-20
    piece residues 19 to 40 of chain A                piece chain B
    piece A/19-20 A/30-40          piece A 19-20+30-40    piece /4CHA//A/19-20
    piece A 19-40 red              piece A 19-40 color=red
"""
from __future__ import annotations

import re

from pymol import cmd

from ..core import report
from . import commands as base
from . import teach

_say = base._say

_RANGE = re.compile(r"^(-?\d+[A-Za-z]?)(?:-(-?\d+[A-Za-z]?))?$")
_LABEL = re.compile(r"^(-?\d+)([A-Za-z]?)$")
_FILLER = set("residue residues res amino acid acids aa of from in and the a an to make "
              "into separate object objects new own its part parts piece pieces only "
              "positions position number numbers with as please copy called just put "
              "turn split extract cut out isolate take for colour color colored coloured "
              "it them that this".split())
_TO = r"(?:-|\.\.|\bto\b|\bthrough\b|\bthru\b|\buntil\b|\btill\b)"
_SHORT = 6          # a ribbon needs a few residues: shorter pieces are drawn as sticks


def _colours():
    try:
        return set(name.lower() for name, _ in cmd.get_color_indices()
                   if not name.startswith("_"))
    except Exception:
        return set("red green blue yellow orange purple pink white grey gray black "
                   "cyan magenta salmon wheat violet".split())


def _tokens(text):
    """What was typed -> words, with '19 to 40' joined into '19-40' and the
    PyMOL form '/4CHA//A/19-20' split into the structure and 'A/19-20'."""
    text = text.replace("–", "-").replace("—", "-")
    text = re.sub(r"(-?\d+[A-Za-z]?)\s*%s\s*(-?\d+[A-Za-z]?)\b" % _TO, r"\1-\2", text)
    text = re.sub(r"(?<=\d)\s*\+\s*(?=-?\d)", "+", text)
    words = []
    for word in re.split(r"[\s,;]+", text.strip()):
        if not word:
            continue
        if "/" in word and not re.match(r"^[A-Za-z]+=", word):
            leading = word.startswith("/")
            parts = word.strip("/").split("/") if leading else word.split("/")
            if leading:                                  # /model/segi/chain/resi
                parts = (parts + ["", "", "", ""])[:4]
                model, chain, resi = parts[0], parts[2], parts[3]
            elif len(parts) >= 3:                        # model/chain/resi
                model, chain, resi = parts[0], parts[-2], parts[-1]
            else:                                        # chain/resi
                model, chain, resi = "", parts[0], parts[-1]
            if model:
                words.append(model)
            if chain:
                words.append("chain=%s" % chain)
            words.extend(r for r in resi.split("+") if r)
            continue
        match = re.match(r"^([A-Za-z0-9]{1,4}):(-?\d.*)$", word)      # A:19-20
        if match:
            words.append("chain=%s" % match.group(1))
            word = match.group(2)
        if "+" in word and all(_RANGE.match(p) for p in word.split("+") if p):
            words.extend(p for p in word.split("+") if p)
        else:
            words.append(word)
    return words


def is_name(word):
    """Could this word (typed after `piece`) be part of a structure's name?"""
    low = word.lower()
    return not (len(word) <= 2 or _RANGE.match(word) or "=" in word or "/" in word
                or low in _FILLER or low in ("chain", "chains") or low in _colours())


def _number(label):
    match = _LABEL.match(label)
    return (int(match.group(1)), match.group(2)) if match else None


def _chains_of(obj):
    """{chain: [residue labels in file order]} for the polymer of ``obj``, and
    the chains of everything else that is not water."""
    residues, seen = {}, set()
    rows = []
    cmd.iterate("%s and polymer" % obj, "rows.append((chain, resi))", space={"rows": rows})
    for chain, resi in rows:
        if (chain, resi) not in seen:
            seen.add((chain, resi))
            residues.setdefault(chain, []).append(resi)
    others = set()
    cmd.iterate("%s and not polymer and not solvent" % obj, "others.add(chain)",
                space={"others": others})
    return residues, others


def _spans(labels):
    """['16', '17', ..., '146', '149', ...] -> '16-146, 149-245' (gaps shown)."""
    if not labels:
        return "-"
    spans, start, last = [], labels[0], labels[0]
    for label in labels[1:]:
        a, b = _number(last), _number(label)
        if a and b and (b[0] - a[0] <= 1):
            last = label
            continue
        spans.append(start if start == last else "%s-%s" % (start, last))
        start = last = label
    spans.append(start if start == last else "%s-%s" % (start, last))
    if len(spans) > 4:
        spans = spans[:3] + ["...", spans[-1]]
    return ", ".join(spans)


def _chain_list(residues, others):
    parts = ["%s (%s)" % (chain, _spans(labels)) for chain, labels in residues.items()]
    parts += ["%s (no protein)" % chain for chain in sorted(others - set(residues))
              if chain]
    if len(parts) > 8:
        parts = parts[:8] + ["... %d more" % (len(parts) - 8)]
    return ", ".join(parts)


def _maps(annotation):
    """{chain: sifts.ChainMap} from the UniProt annotation, when there is one."""
    found = {}
    for protein in getattr(annotation, "proteins", None) or []:
        for chain_map in protein.chains:
            found.setdefault(chain_map.chain, chain_map)
    return found


def _uniprot_of(chain_map, label):
    """The UniProt number of the residue with file label ``label``, or None."""
    for unp_start, unp_end, _seq in chain_map.segments:
        for position in range(unp_start, unp_end + 1):
            if chain_map.resi(position) == label:
                return position
    return None


def _parse(obj, words, residues, others):
    """-> (segments, colour, error).  segments: [[chain or None, [(from, to)]]]."""
    chains = list(residues) + sorted(others - set(residues))
    colours = _colours()
    segments, colour = [], None
    expect, keep = False, False

    def chain_named(word):
        if word in chains:
            return word
        same = [c for c in chains if c.lower() == word.lower()]
        return same[0] if len(same) == 1 else None

    for index, word in enumerate(words):
        low = word.lower()
        following = words[index + 1] if index + 1 < len(words) else ""
        if low in ("chain", "chains"):
            expect, keep = True, low == "chains"
            continue
        if "=" in word:
            key, _, value = word.partition("=")
            key = key.lower()
            if key in ("color", "colour"):
                if value.lower() not in colours:
                    return None, None, ("'%s' is not a colour PyMOL knows. Try red, blue, "
                                        "green, yellow, orange, magenta, cyan ..." % value)
                colour = value.lower()
                continue
            if key == "chain":
                word, low = value, value.lower()
                expect = True
            else:
                return None, None, "'%s' is not something  piece  understands." % word
        match = _RANGE.match(word)
        if match and not expect:
            start, end = match.group(1), match.group(2) or match.group(1)
            if not segments or segments[-1][2]:      # "19-40 of chain A and 50-60 ..."
                segments.append([None, [], False])
            segments[-1][1].append((start, end))
            continue
        named = chain_named(word)
        lone_a = low == "a" and word == "a" and "a" not in chains
        if named and (expect or not lone_a or _RANGE.match(following)):
            if segments and segments[-1][0] is None and not segments[-1][1]:
                segments.pop()
            if segments and segments[-1][0] is None and segments[-1][1]:
                segments[-1][0] = named          # "residues 19 to 40 of chain A"
                segments[-1][2] = True
            else:
                segments.append([named, [], False])
            expect = keep
            continue
        if expect and low not in ("and", "&"):
            return None, None, ("%s has no chain '%s'. Its chains are: %s."
                                % (obj, word, _chain_list(residues, others)))
        if re.match(r"^[A-Za-z0-9]{1,4}$", word) and (_RANGE.match(following) or (
                len(word) == 1 and word.isupper())):
            return None, None, ("%s has no chain '%s'. Its chains are: %s."
                                % (obj, word, _chain_list(residues, others)))
        if low in colours and len(low) > 2:
            colour = low
            continue
        if low in _FILLER or low in ("&",):
            continue
        return None, None, ("I did not understand '%s'. Say which part, e.g.  piece %s "
                            "%s 10-40  or  piece %s chain %s" % (
                                word, obj, chains[0] if chains else "A", obj,
                                chains[-1] if chains else "B"))
    return [s[:2] for s in segments], colour, None


def _within(labels, start, end):
    """The labels of ``labels`` from ``start`` to ``end`` (numbers, inclusive)."""
    a, b = _number(start), _number(end)
    if a is None or b is None:
        return []
    if a > b:
        a, b = b, a
    found = []
    for label in labels:
        n = _number(label)
        if n is None:
            continue
        if (a[0], a[1] if a[1] else "") <= (n[0], n[1]) <= (b[0], b[1] if b[1] else "~"):
            found.append(label)
    return found


def _escape(label):
    return ("\\" + label) if label.startswith("-") else label


def _uniprot_labels(chain_map, start, end, labels):
    """Residues ``start``-``end`` read as UniProt numbers -> file labels that exist."""
    a, b = _number(start), _number(end)
    if a is None or b is None or a[1] or b[1]:
        return []
    low, high = sorted((a[0], b[0]))
    present = set(labels)
    out = []
    for position in range(low, high + 1):
        label = chain_map.resi(position)
        if label and label in present:
            out.append(label)
    return out


def make(obj, words, extra):
    """The action of `piece` (called by tools._run with the structure found)."""
    residues, others = _chains_of(obj)
    if not residues and not others:
        _say("%s has no protein chains to take a piece of." % obj)
        return
    segments, colour, error = _parse(obj, words, residues, others)
    if error:
        _say(error)
        return
    if not segments:
        first = list(residues)[0] if residues else sorted(others)[0]
        _say("Which part of %s? Say a chain and (if you like) residue numbers, e.g." % obj)
        _say("   piece %s %s 10-40      piece %s chain %s"
             % (obj, first, obj, list(residues)[-1] if residues else first))
        _say("CHAINS       %s" % _chain_list(residues, others))
        return
    maps = _maps(extra.get("annotation"))
    parts, notes, names, described = [], [], [], []
    count = 0
    for chain, ranges in segments:
        if not ranges:                                       # a whole chain
            parts.append("chain %s" % chain)
            described.append("chain %s" % chain)
            names.append(chain)
            count += len(residues.get(chain, []))
            continue
        if chain is None:
            if len(residues) == 1:
                chain = list(residues)[0]
            else:
                having = [c for c in residues
                          if all(_within(residues[c], s, e) for s, e in ranges)]
                if not having:
                    _say("No chain of %s has residues %s." % (
                        obj, ", ".join(s if s == e else "%s-%s" % (s, e) for s, e in ranges)))
                    _say("CHAINS       %s" % _chain_list(residues, others))
                    return
                chain = having[0]
                if len(having) > 1:
                    notes.append("Chains %s have these numbers; I took chain %s. For "
                                 "another:  piece %s %s %s" % (
                                     " and ".join([", ".join(having[:-1]), having[-1]]),
                                     chain, obj, having[1],
                                     " ".join(s if s == e else "%s-%s" % (s, e)
                                              for s, e in ranges)))
        labels = residues.get(chain, [])
        if not labels:
            _say("Chain %s of %s holds no protein residues (only %s). Take the whole "
                 "chain:  piece %s chain %s" % (chain, obj, "bound molecules", obj, chain))
            return
        chosen = []
        chain_map = maps.get(chain)
        for start, end in ranges:
            asked = start if start == end else "%s-%s" % (start, end)
            found = _within(labels, start, end)
            if found:
                a, b = _number(start), _number(end)
                if a and b and not a[1] and not b[1] and \
                        len(set(_number(f)[0] for f in found)) < abs(b[0] - a[0]) + 1:
                    notes.append("Chain %s has only residues %s of %s in this file (the "
                                 "rest were not seen in the experiment, or are not part "
                                 "of this structure)." % (chain, _spans(found), asked))
                first = _uniprot_of(chain_map, found[0]) if chain_map else None
                last = _uniprot_of(chain_map, found[-1]) if chain_map else None
                if first and last and (str(first), str(last)) != (found[0], found[-1]):
                    note = ("Numbers are the file's; UniProt numbers the same residues %s."
                            % (first if first == last else "%d-%d" % (first, last)))
                    other = _uniprot_labels(chain_map, start, end, labels)
                    if other and other != found:
                        note += (" If you meant UniProt's %s:  piece %s %s %s"
                                 % (asked, obj, chain, _spans(other).replace(", ", " ")))
                    notes.append(note)
            elif chain_map:
                found = _uniprot_labels(chain_map, start, end, labels)
                if found:
                    notes.append("%s is UniProt numbering: in this file those residues "
                                 "are numbered %s." % (asked, _spans(found)))
            if not found:
                _say("Chain %s of %s has no residue%s %s." % (
                    chain, obj, "" if start == end else "s", asked))
                _say("CHAINS       %s" % _chain_list(residues, others))
                if chain_map:
                    first = _uniprot_of(chain_map, labels[0])
                    last = _uniprot_of(chain_map, labels[-1])
                    if first and last:
                        _say("UNIPROT      chain %s holds UniProt residues %d-%d of %s "
                             "(file numbers can differ from UniProt's)."
                             % (chain, first, last, chain_map.accession))
                return
            chosen.extend(f for f in found if f not in chosen)
            names.append("%s%s" % (chain, found[0] if len(found) == 1
                                   else "%s-%s" % (found[0], found[-1])))
        described.append("chain %s, residue%s %s" % (chain, "" if len(chosen) == 1 else "s",
                                                     _spans(chosen)))
        count += len(chosen)
        parts.append("(chain %s and resi %s)" % (chain, "+".join(_escape(f) for f in chosen)))
    selection = "%s and (%s) and not solvent" % (obj, " or ".join(parts))
    name = "%s_%s" % (obj, "_".join(names))
    if len(name) > 40:
        name = "%s_piece" % obj
    name = cmd.get_legal_name(name)
    existing = cmd.get_names("all")
    if name in existing and name not in _made:
        number = 2
        while "%s_%d" % (name, number) in existing:
            number += 1
        name = "%s_%d" % (name, number)
    if name in existing:
        cmd.delete(name)
    cmd.create(name, selection)
    _made[name] = obj
    identity = base._identify(obj)
    if identity:
        base._objects[name] = identity
    cmd.hide("everything", selection)                # drawn once: by the new object
    sticks = count <= _SHORT
    if sticks:
        cmd.show("sticks", "%s and not hydro" % name)
    if colour:
        cmd.color(colour, name)
        if sticks:
            cmd.util.cnc(name)
    cmd.enable(name)
    teach.note("create %s, %s" % (name, selection), "copy that part into a new object")
    teach.note("hide everything, %s" % selection,
               "hide it in the original, so it is not drawn twice")
    if colour:
        teach.note("color %s, %s" % (colour, name), "colour the new object")
    what = "; ".join(described)
    _say("MADE         a new object  %s  (%s of %s, %d residue%s)"
         % (name, what, obj, count, "" if count == 1 else "s"))
    _say("             It is in the object list on the right. %s keeps every atom: those "
         "residues are only hidden there, so you do not see them twice." % obj)
    if sticks:
        _say("             Drawn as sticks: a ribbon needs more than %d residues." % _SHORT)
    if colour:
        _say("COLOURED     %s. Another colour:  color %s, %s"
             % (colour, "blue" if colour != "blue" else "red", name))
    else:
        _say("COLOUR IT    color red, %s" % name)
    _say("HIDE IT      disable %s   (back on screen:  enable %s), or click its name "
         "on the right" % (name, name))
    for note in notes:
        _say("NOTE         %s" % note)
    _say("WHAT NEXT    viewas %s sticks  |  piece %s chain %s  |  back  (removes the "
         "new object)" % (name, obj, (list(residues) or sorted(others))[-1]))


_made = {}       # piece object -> the structure it was taken from


def piece(*args, **kwargs):
    """
DESCRIPTION

    Make part of a structure - a chain, or residues of a chain - into a new
    object of its own. It appears in the object list on the right, and can
    be coloured, hidden or shown like any object. The original structure
    keeps every atom. Numbers are the residue numbers of the file; UniProt
    numbers are recognised when the file has no such residues.

USAGE

    piece [structure] <chain> [from-to ...] [colour]
    piece [structure] chain <chain>

EXAMPLES

    piece A 19-40
    piece chymotrypsinogen A/19-20
    piece residues 19 to 40 of chain A
    piece chain B
    piece A 19-20 30-40 red
    piece /4CHA//B/40-60
    """
    from . import tools
    text = base._raw(args, kwargs)
    words = _tokens(text)
    words = [w for w in words if w.lower() not in ("piece",)]
    label = "piece %s" % text if text else "piece"
    tools._run(label, words, make, is_name=is_name, need_annotation=True)


def register():
    cmd.extend("piece", piece)
    cmd.extend("p_piece", piece)
    try:
        cmd.auto_arg[0]["piece"] = cmd.auto_arg[0]["zoom"]
    except Exception:
        pass
