"""Making sense of what was typed after a command.

People type the same request many ways:

    explain 248        explain R248        explain arginine 248
    measure 176 and 238                    measure A/176 to ZN
    whatif R248Q       whatif 248 Q        whatif arg248gln
    look rainbow p53   look p53 rainbow    scout cysteines in p53

`tidy(command, text)` turns these into the one form each command reads, and
`not_a_name(text)` says why a leftover word must not be looked up as a
protein (the old behaviour loaded a protein called "rainbow" or "and").
"""
from __future__ import annotations

import difflib
import re

from pymol import cmd

from . import commands as base

# protein names that look like "letter + number" (p53 is not proline 53)
NAME_LIKE = set("p53 p21 p63 p73 p38 p16 p27 p300 p97 p51 p66 p24 p17 p14 p19 p62 p85 p110 "
                "p120 p130 p107 h1 h2a h2b h3 h4 e1 e2 e3 c1 c3 c4 c5 c9 f1 g1 s6 l1 a1 a2 "
                "b1 b2 m1 m2 n1 n2 t4 t7 k1 d1 d2 cd4 cd8 cd20 cd28 cd40 il2 il6 il10 cox1 cox2 "
                "her2 erk1 erk2 jak2 kras hras nras".split())
FILLER = set("in of on the to and with from by at for a an between into onto vs versus "
             "is are it its this that then as".split())
# commands whose arguments are a few keywords plus (maybe) a structure's name
_PLAIN = ("scout", "look", "viewas", "explain", "measure", "whatif", "mutate", "sequence",
          "features", "pocket", "tour", "similar", "report")
_RESIDUES = ("explain", "measure", "whatif", "mutate")
_vocab = [None]


def _amino():
    """word -> one-letter code, for names (arginine), 3-letter and 1-letter codes."""
    from . import tools
    table = {}
    for code, (name, _group, _use) in tools.AA.items():
        one = [k for k, v in tools.ONE_TO_THREE.items() if v == code][0]
        for word in (name, name + "s", code.lower()):
            table[word] = one
        table[one.lower()] = one
    table.update({"aspartic": "D", "glutamic": "E"})
    return table


def vocabulary():
    """Every word some command accepts as an argument - never a protein's name."""
    if _vocab[0] is None:
        from . import places, tier3, tools
        words = set(tools._scout_words()) | set(tier3.STYLE_NAMES) | set(tier3.PRESETS)
        words |= set(tier3._PRESET_ALIASES) | set(tier3._STYLE_ALIASES)
        words |= set(places.PLACE_WORDS)
        words |= set("colour color colours colors style look view draw show hide zoom "
                     "red blue green yellow orange purple pink white grey gray black cyan "
                     "magenta first last off stop on next prev all everything none "
                     "change changes changed mutate mutation mutant replace swap make "
                     "turn become becomes spot find highlight select mark measure "
                     "distance explain describe tell compare load get open please what "
                     "which where how many residue residues amino acid acids number "
                     "position chain chains atom atoms picture image save".split())
        _vocab[0] = set(w.lower() for phrase in words for w in phrase.split())
    return _vocab[0]


def _ligand_words():
    from . import places
    return set(w for names in places.LIGAND_WORDS.values() for w in names)


def not_a_name(text, strict=False):
    """None when ``text`` could be the name of a protein or molecule; otherwise
    a sentence saying why it was not searched for.  ``strict`` (for commands
    that have their own keywords): amino acids and bound-molecule words are
    arguments too, and a near-miss of a keyword is a typo, not a name."""
    words = [w for w in text.lower().split()]
    if not words:
        return "nothing was named"
    meaningful = [w for w in words if w not in FILLER]
    if not meaningful:
        return "'%s' is not a name" % text
    known = vocabulary()
    amino = _amino()
    ligands = _ligand_words()

    def argument(word):
        if re.match(r"^-?\d+([a-z])?$", word) or len(word) == 1:
            return True
        if re.match(r"^[a-z0-9]{1,2}/[a-z]{0,3}-?\d+", word):
            return True                                  # A/248
        if word in known:
            return True
        if strict and (word in amino or word in ligands):
            return True
        return False

    if all(argument(w) for w in meaningful):
        return ("'%s' is not a loaded structure and does not look like a protein's name"
                % " ".join(meaningful))
    if strict:
        for word in meaningful:
            if argument(word) or len(word) < 5:
                continue
            close = difflib.get_close_matches(word, sorted(known | set(amino)), n=1,
                                              cutoff=0.82)
            if close:
                return "'%s' looks like a typo of '%s'" % (word, close[0])
    return None


# ------------------------------------------------------------------ residues
def _residues_of(obj):
    """{number: one-letter code} and the chain letters of a loaded structure."""
    from ..core import features as feat
    found, chains = {}, set()
    try:
        rows = []
        cmd.iterate("%s and polymer and name CA" % obj, "rows.append((chain, resi, resn))",
                    space={"rows": rows})
        for chain, resi, resn in rows:
            chains.add(chain)
            found.setdefault(resi, feat.THREE_TO_ONE.get(resn.upper(), "?"))
    except Exception:
        pass
    return found, chains


def _current():
    names = base._molecules()
    if base._last[0] in names:
        return base._last[0]
    return names[-1] if names else None


def _tidy_residues(command, words):
    """Residue mentions in every usual spelling -> the form the command reads."""
    amino = _amino()
    obj = None
    for size in range(min(3, len(words)), 0, -1):       # a structure named in the text
        for start in (0, len(words) - size):
            if base._known_object(" ".join(words[start:start + size])):
                obj = base._known_object(" ".join(words[start:start + size]))
    obj = obj or _current()
    residues, chains = _residues_of(obj) if obj else ({}, set())
    ligands = {}
    if obj and command == "measure":
        from . import places
        codes = set()
        try:
            cmd.iterate("%s and not polymer and not solvent" % obj, "codes.add(resn)",
                        space={"codes": codes})
        except Exception:
            pass
        for code in codes:
            for word in places.LIGAND_WORDS.get(code, []):
                ligands[word] = code
    changing = command in ("whatif", "mutate")
    out, index = [], 0
    while index < len(words):
        word = words[index]
        low = word.lower()
        following = words[index + 1].lower() if index + 1 < len(words) else ""
        after = words[index + 2].lower() if index + 2 < len(words) else ""
        glued_probe = re.match(r"^([a-z])(-?\d+)$", low)
        others = [w for w in words if w is not word and (
            re.match(r"^-?\d+[a-z]?$", w.lower()) or "/" in w
            or re.match(r"^[a-z]{1,3}-?\d+[a-z]{0,3}$", w.lower()))]
        is_residue = bool(glued_probe) and residues.get(glued_probe.group(2)) == \
            amino.get(glued_probe.group(1))
        if low in NAME_LIKE and (not is_residue or others):
            out.append(word)                      # p53, kras: a protein, not a residue
            index += 1
            continue
        if base._known_object(word) or "/" in word:
            out.append(word)
            index += 1
            continue
        if low in ligands:                                # "the zinc" -> ZN
            out.append(ligands[low])
            index += 1
            continue
        glued = re.match(r"^([a-z]+)(-?\d+)([a-z]*)$", low)
        if changing:
            # arg248gln / r248q written out, glued or spaced
            if glued and glued.group(1) in amino and glued.group(3) in amino:
                out.append("%s%s%s" % (amino[glued.group(1)], glued.group(2),
                                       amino[glued.group(3)]))
                index += 1
                continue
            if glued and glued.group(1) in amino and not glued.group(3):
                number = glued.group(2)
                complete = any(re.match(r"^(?:[a-z0-9]{1,4}/)?[a-z]\d+[a-z]$", w.lower())
                               for w in words if w is not word)
                step = index + 1
                if following in ("to", "into", "becomes", "->", "=>", ">", "with", "as"):
                    step += 1
                target = words[step].lower() if step < len(words) else ""
                if target in amino:
                    out.append("%s%s%s" % (amino[glued.group(1)], glued.group(2),
                                           amino[target]))
                    index = step + 1
                    continue
                # R248 alone ("becomes what?") - but only when it really is a residue
                # of the structure on screen; otherwise it may be a name (p53)
                letter = amino[glued.group(1)]
                if residues.get(number) == letter and not complete:
                    out.append(number)
                    index += 1
                    continue
                out.append(word)
                index += 1
                continue
            if low in amino and re.match(r"^-?\d+$", following):
                step = index + 2                          # arginine 248 [to] glutamine
                if after in ("to", "into", "becomes", "->", "=>", ">", "with", "as"):
                    step += 1
                target = words[step].lower() if step < len(words) else ""
                if target in amino:
                    out.append("%s%s%s" % (amino[low], following, amino[target]))
                    index = step + 1
                else:
                    out.append(following)
                    index += 2
                continue
            if re.match(r"^-?\d+$", low):                 # 248 [to] Q
                step = index + 1
                if following in ("to", "into", "becomes", "->", "=>", ">", "with", "as"):
                    step += 1
                target = words[step].lower() if step < len(words) else ""
                if target in amino and (len(target) > 1 or step == len(words) - 1
                                        or True):
                    out.append("%s%s" % (low, amino[target]))
                    index = step + 1
                    continue
            out.append(word)
            index += 1
            continue
        # explain / measure: the residue's number is what matters
        if low in amino and len(low) > 1 and re.match(r"^-?\d+[a-z]?$", following):
            out.append(following.upper())                 # arginine 248 -> 248
            index += 2
            continue
        if glued and not glued.group(3) and glued.group(1) in amino:
            number, letter = glued.group(2), amino[glued.group(1)]
            if residues.get(number) == letter or (len(glued.group(1)) == 3
                                                  and number in residues):
                out.append(number)                        # R248, arg248 -> 248
                index += 1
                continue
            if len(glued.group(1)) == 1 and glued.group(1).upper() in chains \
                    and number in residues:
                out.append("%s/%s" % (glued.group(1).upper(), number))   # A248 -> A/248
                index += 1
                continue
        out.append(word)
        index += 1
    return out


def tidy(command, text):
    """What was typed after ``command``, in the form the command reads."""
    if command == "compare":
        if "," not in text:
            text = re.sub(r"\s+(?:with|and|vs\.?|versus|to|onto|against)\s+", " , ", text,
                          count=1, flags=re.I)
        return text
    if command not in _PLAIN or not text.strip():
        return text
    if base._known_object(text):
        return text
    words = text.split()
    if command in _RESIDUES:
        words = _tidy_residues(command, words)
    noise = set("show please find highlight spot change mutate make see me display".split())
    if command in _RESIDUES:
        noise |= set("residue residues amino acid number position".split())
    kept = []
    for index, word in enumerate(words):
        if word.lower() in noise and not base._known_object(word) and len(words) > 1:
            continue
        if word.lower() in FILLER and not base._known_object(word):
            # keep a filler word that is part of a known name ("protein of ...")
            continue
        kept.append(word)
    if command == "look":
        kept = [w for w in kept if w.lower() not in ("colour", "color", "coloured",
                                                     "colored", "scheme", "style")]
    return " ".join(kept)
