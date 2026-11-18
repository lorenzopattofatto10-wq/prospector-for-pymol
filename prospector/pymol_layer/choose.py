"""`choose`: describe atoms in plain words, get a PyMOL selection - and see the
selection line PyMOL itself would need, with each piece explained.

    choose histidines within 5 of the zinc
    choose charged residues in chain A
    choose residues 100 to 120 except glycines

`translate()` is pure text -> selection expression and can be tested without
a structure; `choose()` is the PyMOL command.
"""
from __future__ import annotations

import re

from pymol import cmd

from ..core import report
from . import commands as base
from . import places, teach, tools, worker

_say = base._say


class Unknown(Exception):
    """A word in the description that `choose` has no meaning for."""


# word(s) -> (kind, expression, what it means).  Kinds: "what" things are added
# together (histidines and cysteines); "where" things narrow them down (in
# chain A, in helices).
VOCABULARY = {}


def _add(words, kind, expression, meaning):
    for word in words:
        VOCABULARY[word] = (kind, expression, meaning)


for _code, (_name, _group, _use) in tools.AA.items():
    _add((_name, _name + "s", _code.lower()), "resn", "resn %s" % _code,
         "resn %s = residues named %s (%s)" % (_code, _code, _name))
for _words, _codes, _label in [
        (("hydrophobic", "oily", "nonpolar", "apolar", "water-avoiding"),
         "ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO", "the water-avoiding amino acids"),
        (("polar", "hydrophilic"), "SER+THR+ASN+GLN+TYR+CYS+HIS", "the polar amino acids"),
        (("charged",), "ASP+GLU+LYS+ARG", "the charged amino acids"),
        (("positive", "basic", "positively charged"), "LYS+ARG+HIS",
         "the positively charged amino acids"),
        (("negative", "acidic", "negatively charged"), "ASP+GLU",
         "the negatively charged amino acids"),
        (("aromatic", "aromatics", "rings"), "PHE+TYR+TRP+HIS", "the ring amino acids"),
        (("small", "tiny"), "GLY+ALA+SER", "the smallest amino acids")]:
    _add(_words, "resn", "resn %s" % _codes, "resn %s = %s" % (_codes, _label))
_add(("helix", "helices", "alpha helix", "alpha helices", "alpha-helices"), "ss", "ss h",
     "ss h = secondary structure 'helix'")
_add(("strand", "strands", "sheet", "sheets", "beta strands", "beta sheets", "beta sheet",
      "beta strand"), "ss", "ss s", "ss s = secondary structure 'strand'")
_add(("loop", "loops", "coil", "coils", "turns"), "ss", "not ss h+s",
     "not ss h+s = neither helix nor strand, i.e. loops")
_add(("ligand", "ligands", "drug", "drugs", "inhibitor", "cofactor", "cofactors",
      "substrate", "small molecule", "small molecules", "bound molecule", "bound molecules"),
     "what", "organic", "organic = bound small molecules (not protein, not water)")
_add(("metal", "metals", "metal ion", "metal ions"), "what", "metals",
     "metals = metal atoms")
_add(("ion", "ions"), "what", "inorganic", "inorganic = ions and other inorganic atoms")
_add(("water", "waters", "water molecules", "solvent"), "what", "solvent",
     "solvent = water molecules")
_add(("protein", "polypeptide", "amino acids"), "what", "polymer.protein",
     "polymer.protein = the protein chains")
_add(("dna", "rna", "nucleic acid", "nucleic acids", "nucleotides"), "what",
     "polymer.nucleic", "polymer.nucleic = DNA and RNA chains")
_add(("everything", "anything", "all atoms"), "what", "all", "all = every atom")
_add(("backbone", "main chain", "mainchain"), "part", "backbone",
     "backbone = the N, CA, C, O atoms of each residue")
_add(("side chain", "side chains", "sidechain", "sidechains"), "part", "sidechain",
     "sidechain = everything of a residue except its backbone")
_add(("alpha carbon", "alpha carbons", "ca atoms", "c-alpha", "c-alphas"), "part", "name CA",
     "name CA = the alpha-carbon atom of each residue")
_add(("hydrogen", "hydrogens"), "part", "hydro", "hydro = hydrogen atoms")
for _word, _elem in (("sulfur", "S"), ("sulphur", "S"), ("oxygen", "O"), ("oxygens", "O"),
                     ("nitrogen", "N"), ("nitrogens", "N"), ("carbon", "C"),
                     ("carbons", "C"), ("phosphorus", "P")):
    _add((_word,), "part", "elem %s" % _elem, "elem %s = atoms of the element %s"
         % (_elem, _elem))

_FILLER = set("the a an and in of on from with that which are is to at all every each any "
              "only residue residues atom atoms amino acid acids part parts molecule "
              "molecules inside located sitting me my please select choose pick find show "
              "those ones".split())
_NEAR = [(r"\bwithin\s+(\d+(?:\.\d+)?)\s*(?:a|å|angstroms?|ångströms?)?\s*(?:of|from|to)\b",
          None),
         (r"\b(?:closer than|less than|at most|up to)\s+(\d+(?:\.\d+)?)\s*(?:a|å|angstroms?)?"
          r"\s*(?:of|from|to|away from)\b", None),
         (r"\b(?:touching|in contact with|contacting|binding|bound to|holding|that hold|"
          r"that bind|that touch|coordinating)\b", 4.0),
         (r"\b(?:near|around|close to|next to|surrounding|by)\b", 5.0)]
_NOT = r"\b(?:except|but not|without|excluding|other than|apart from|not)\b"


def _terms(text, codes, names=None):
    """The vocabulary entries in ``text``, in order.  Raises Unknown.
    ``names`` maps what the user calls a bound molecule (cholesterol) to its
    code in the structure (CLR)."""
    names = names or {}
    words = re.findall(r"[a-z0-9][a-z0-9\-']*", text)
    found, index = [], 0
    while index < len(words):
        word = words[index]
        # chain A / chains A and B
        if word in ("chain", "chains") and index + 1 < len(words):
            letters, step = [], index + 1
            while step < len(words) and (len(words[step]) == 1 or words[step] in ("and", "or")):
                if len(words[step]) == 1:
                    letters.append(words[step].upper())
                step += 1
            if letters:
                found.append(("chain", "chain %s" % "+".join(letters),
                              "chain %s = the chain%s with that letter"
                              % ("+".join(letters), "s" if len(letters) > 1 else "")))
                index = step
                continue
        # numbers: 57 / 100-120 / 100 to 120 / 10 12 and 15
        if re.match(r"^\d+(-\d+)?$", word):
            numbers, step = [word], index + 1
            if step + 1 < len(words) and words[step] in ("to", "through", "until") and \
                    words[step + 1].isdigit():
                numbers = ["%s-%s" % (word, words[step + 1])]
                step += 2
            while step < len(words) and (re.match(r"^\d+(-\d+)?$", words[step])
                                         or (words[step] in ("and", "or")
                                             and step + 1 < len(words)
                                             and words[step + 1][0].isdigit())):
                if words[step][0].isdigit():
                    numbers.append(words[step])
                step += 1
            found.append(("resi", "resi %s" % "+".join(numbers),
                          "resi %s = residue number%s %s" % (
                              "+".join(numbers), "" if len(numbers) == 1
                              and "-" not in numbers[0] else "s",
                              ", ".join(n.replace("-", " to ") for n in numbers))))
            index = step
            continue
        # the longest phrase of the vocabulary that starts here
        for size in (3, 2, 1):
            phrase = " ".join(words[index:index + size])
            if len(words[index:index + size]) == size and phrase in VOCABULARY:
                found.append(VOCABULARY[phrase])
                index += size
                break
            if len(words[index:index + size]) == size and phrase in names:
                found.append(("what", "resn %s" % names[phrase],
                              "resn %s = the bound molecule with the code %s (%s)"
                              % (names[phrase], names[phrase], phrase)))
                index += size
                break
            # a bound molecule by name (zinc) or by its code in this structure (ZN)
            if size == 1:
                code = None
                if word.upper() in codes:
                    code = word.upper()
                else:
                    for candidate, names in places.LIGAND_WORDS.items():
                        if word in names and (candidate in codes or not codes):
                            code = candidate
                            break
                if code:
                    found.append(("what", "resn %s" % code,
                                  "resn %s = the bound molecule with the code %s"
                                  % (code, code)))
                    index += 1
                    break
        else:
            if word in _FILLER:
                index += 1
                continue
            raise Unknown(word)
    return found


def _combine(terms):
    """Things of the same kind are added together; different kinds narrow each
    other down.  Returns (expression, [explanations])."""
    order, groups, notes = [], {}, []
    for kind, expression, meaning in terms:
        if kind not in groups:
            groups[kind] = []
            order.append(kind)
        if expression not in groups[kind]:
            groups[kind].append(expression)
        if meaning not in notes:
            notes.append(meaning)
    pieces = []
    for kind in order:
        items = groups[kind]
        if kind in ("resn", "chain", "resi") and len(items) > 1:
            head = items[0].split()[0]
            merged = "+".join(item.split(" ", 1)[1] for item in items)
            pieces.append("%s %s" % (head, merged))
        elif len(items) > 1:
            pieces.append("(%s)" % " or ".join(items))
        else:
            pieces.append(items[0])
    # "histidines and the zinc": two kinds of thing, both wanted
    things = [p for k, p in zip(order, pieces) if k in ("resn", "what")]
    narrow = [p for k, p in zip(order, pieces) if k not in ("resn", "what")]
    expression = things[0] if len(things) == 1 else (
        "(%s)" % " or ".join(things) if things else "")
    if len(things) > 1:
        notes.append("or = either one")
    for piece in narrow:
        expression = "%s and %s" % (expression, piece) if expression else piece
    if expression.count(" and "):
        notes.append("and = both must be true")
    return expression, notes


def _set(text, codes, names=None):
    """A description without a distance in it -> (expression, notes)."""
    parts = re.split(_NOT, text, maxsplit=1)
    wanted = _terms(parts[0], codes, names)
    if not wanted and len(parts) > 1:
        wanted = [VOCABULARY["protein"]]
    if not wanted:
        raise Unknown(text.strip() or "(nothing)")
    expression, notes = _combine(wanted)
    if len(parts) > 1:
        other, more = _set(parts[1], codes, names)
        expression = "%s and not (%s)" % (expression, other)
        notes += [n for n in more if n not in notes]
        notes.append("not = leave these out")
    return expression, notes


def translate(text, codes=(), scope=None, names=None):
    """Plain words -> (PyMOL selection expression, [explanations]).
    ``codes`` are the bound-molecule codes of the structure (ZN, HEM ...);
    ``scope`` is the object the selection is limited to.
    Raises Unknown(word) for a word it has no meaning for."""
    def inside(expression):
        return "%s and (%s)" % (scope, expression) if scope else expression

    text = " " + " ".join(text.lower().replace(",", " ").replace("?", " ").split()) + " "
    codes = set(c.upper() for c in codes)
    for pattern, default in _NEAR:
        match = re.search(pattern, text)
        if not match:
            continue
        distance = float(match.group(1)) if default is None else default
        left, right = text[:match.start()], text[match.end():]
        target, notes = _set(right, codes, names)
        whole = not re.search(r"\batoms?\b", left)
        try:
            what, more = _set(left, codes, names)
        except Unknown:
            if _terms(left, codes, names) == [] and not re.search(_NOT, left):
                what, more = "polymer", ["polymer = the protein (and DNA/RNA) chains"]
            else:
                raise
        number = ("%g" % distance)
        expression = "(%s) within %s of (%s)" % (inside(what), number, inside(target))
        if what in ("polymer", "all"):
            expression = "(%s) and not (%s)" % (expression, target)
        notes = more + [n for n in notes if n not in more]
        notes.append("within %s of = closer than %s Angstrom to" % (number, number))
        if whole:
            expression = "byres (%s)" % expression
            notes.append("byres = take whole residues, not single atoms")
        return expression, notes
    expression, notes = _set(text, codes, names)
    return inside(expression), notes


# ----------------------------------------------------------------- the command
def choose(*args, **kwargs):
    """
DESCRIPTION

    Select atoms by describing them in plain words.  Makes a PyMOL selection
    (called "chosen" unless you name it), marks it on screen, and prints the
    selection line PyMOL itself needs, with every piece explained - so you
    learn PyMOL's selection language by reading it.

USAGE

    choose <description> [as <name>]
    choose                              list the words it understands

    Things:  any amino acid, hydrophobic, polar, charged, positive, negative,
             aromatic, ligand, metals, ions, water, protein, dna, a bound
             molecule by name or code (zinc, heme, ZN)
    Where:   chain A, residues 100 to 120, helices, strands, loops,
             backbone, side chains, alpha carbons
    Near:    within 5 of ..., near ..., touching ...
    Leave out:  except ..., without ...

EXAMPLES

    choose histidines within 5 of the zinc
    choose charged residues in chain A
    choose residues 100 to 120 except glycines
    choose everything touching the ligand as pocket
    choose the helices of chain B
    """
    text = base._raw(args, kwargs)
    label = "choose %s" % text
    words = text.split()
    name = "chosen"
    if len(words) >= 2 and words[-2].lower() in ("as", "called", "named"):
        name = re.sub(r"[^A-Za-z0-9_]", "_", words[-1])
        words = words[:-2]
    obj = None
    if len(words) >= 2 and words[-2].lower() in ("in", "of") and \
            base._known_object(words[-1]):
        obj = base._known_object(words[-1])
        words = words[:-2]
    elif words and base._known_object(words[0]) and len(words) > 1:
        obj = base._known_object(words[0])
        words = words[1:]
    # "... in p53" with p53 not loaded yet: find and load it first, like every
    # other command.  Only when the words after the last in/of are not part of
    # the description itself ("in chain A", "of the zinc").
    query = ""
    if obj is None:
        lowered = [w.lower() for w in words]
        for index in range(len(words) - 2, 0, -1):
            if lowered[index] in ("in", "of", "from"):
                tail = " ".join(words[index + 1:])
                try:
                    understood = _terms(tail.lower(), set())
                except Unknown:
                    understood = None
                if understood is None and not re.search(
                        r"\b(within|near|around|touching|except|without)\b", tail.lower()):
                    query, words = tail, words[:index]
                break
    description = " ".join(words)
    fetched = {}

    def action():
        base._header(label)
        _say("")
        if not description:
            _say("Describe what you want after  choose . Words I understand:")
            _say("   THINGS   an amino acid (histidines, cys), hydrophobic, polar, charged, "
                 "positive, negative, aromatic, ligand, metals, ions, water, protein, dna, "
                 "a bound molecule (zinc, heme, ZN)")
            _say("   WHERE    chain A, residues 100 to 120, helices, strands, loops, "
                 "backbone, side chains, alpha carbons")
            _say("   NEAR     within 5 of ...   near ...   touching ...")
            _say("   NOT      except ...   without ...")
            _say("   NAME IT  ... as pocket")
            _say("EXAMPLE      choose histidines within 5 of the zinc")
            _say(report.LINE)
            return
        target = obj
        if fetched.get("result"):
            if fetched["result"]["error"]:
                _say(fetched["result"]["error"])
                _say(report.LINE)
                return
            target = base._place_fetched(fetched["result"])
        if target is None:
            target, message = base._pick_object("", "choose")
            if target is None:
                _say(message)
                _say(report.LINE)
                return
        codes = set()
        cmd.iterate("%s and not polymer and not solvent" % target, "codes.add(resn)",
                    space={"codes": codes})
        # what the user may call the bound molecules: everyday names, and the
        # names of small molecules loaded earlier (cholesterol -> CLR)
        from ..core import chemcomp
        names = dict((word, code) for word, code in chemcomp.COMMON.items()
                     if code in codes)
        for alias, loaded in base._aliases.items():
            identity = base._objects.get(loaded)
            if identity and identity[0] == "mol" and identity[1] in codes:
                names[alias] = identity[1]
        tools._ensure_ss(target)
        try:
            expression, notes = translate(description, codes, target, names)
        except Unknown as error:
            _say("I don't know what '%s' means here, so I selected nothing." % error)
            _say("Type  choose  alone for the words I understand, or write PyMOL's own "
                 "line:  select mine, <expression>")
            if codes:
                _say("Bound molecules in %s: %s" % (target, ", ".join(sorted(codes))))
            _say(report.LINE)
            return
        full = expression
        base._snapshot(label)
        try:
            count = cmd.select(name, full)
        except Exception as error:      # a selection PyMOL rejects: say so, do nothing
            _say("PyMOL did not accept the selection I built (%s)." % str(error)[:80])
            _say("PYMOL        select %s, %s" % (name, full))
            _say(report.LINE)
            return
        base._last[0] = target
        residues = []
        cmd.iterate("%s and (name CA or not polymer)" % name,
                    "residues.append((chain, resn, resi))", space={"residues": residues})
        residues = list(dict.fromkeys(residues))
        if not count:
            _say("CHOSEN       nothing in %s matches that." % target)
            cmd.delete(name)
        else:
            shown = ", ".join("%s%s%s" % (c + "/" if c else "", n, i)
                              for c, n, i in residues[:10])
            _say("CHOSEN       %d atom%s in %d residue%s of %s%s"
                 % (count, "" if count == 1 else "s", len(residues),
                    "" if len(residues) == 1 else "s", target,
                    ": %s%s" % (shown, " ..." if len(residues) > 10 else "")
                    if residues else ""))
            cmd.enable(name)
        _say("PYMOL        select %s, %s" % (name, full))
        _say("READ IT      %s = only in this structure" % target)
        for note in notes:
            _say("             %s" % note)
        if count:
            _say("ON SCREEN    small pink squares mark the chosen atoms; the selection "
                 "is called '%s' in the list on the right" % name)
            _say("USE IT       show sticks, %s  |  color yellow, %s  |  zoom %s"
                 % (name, name, name))
            _say("             deselect  (hide the squares)  |  back  (remove the "
                 "selection again)")
        teach.note("select %s, %s" % (name, full), "a named selection you can reuse")
        _say(report.LINE)

    def job():
        if query and description:
            fetched["result"] = base._fetch_query(query)

    worker.submit(job, lambda result, error: action() if not error
                  else base._unexpected(error))


def register():
    for name in ("choose", "p_choose"):
        cmd.extend(name, choose)
