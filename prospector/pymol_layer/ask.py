"""Plain-language front door: ask, start, status, glossary.

`ask` turns an everyday sentence into the Prospector command that does it,
SHOWS that command, and runs it.  It uses fixed rules, not an AI service: it is
free, works the same every time, and when it is not sure it says so and offers
choices instead of guessing.
"""
from __future__ import annotations

import difflib
import re

from pymol import cmd

from ..core import report
from . import commands as base
from . import teach, tier3, tools, worker

_say = base._say

# ----------------------------------------------------------------- glossary
GLOSSARY = {
    "colours": "on atoms shown as sticks or balls the colour is the element: red = "
               "oxygen, blue = nitrogen, yellow = sulfur, white = hydrogen, and carbon "
               "takes the colour of the highlight (orange, green ...). Ribbon colours "
               "depend on the last command: its ON SCREEN lines say what each colour "
               "means, and  status  repeats how the structure is drawn now.",
    "pdb": "the Protein Data Bank: the public archive of experimentally determined 3D "
           "structures. Every entry has a four-character ID such as 4HHB.",
    "pdb id": "the four-character code of one Protein Data Bank entry, e.g. 4HHB or 1TUP.",
    "uniprot": "the reference database of protein sequences and what is known about them "
               "(function, sites, variants). Each protein has an accession such as P04637.",
    "resolution": "how fine the detail of an experimental structure is, in Angstrom. "
                  "Smaller is sharper: below 2 A you see individual atoms clearly; above "
                  "3.5 A only the overall shape of the chain.",
    "angstrom": "the unit of length for atoms: 1 A = 0.1 nanometre. A carbon-carbon bond "
                "is 1.5 A; a hydrogen bond is about 2.8-3.2 A.",
    "rmsd": "root-mean-square deviation: the average distance between matching atoms of "
            "two superposed structures. Below 1 A = essentially identical; 1-2.5 A = same "
            "fold with local differences; above 5 A = very different.",
    "r-free": "a quality check for crystal structures: how well the model predicts data "
              "that were kept aside. Lower is better; roughly 0.20-0.25 is normal.",
    "r-work": "how well a crystal-structure model fits the data used to build it. Lower is "
              "better; it should be close to (a little below) R-free.",
    "b-factor": "how smeared out an atom is in the experiment. High values mean the atom "
                "moves or is poorly defined (loops, surface); low values mean rigid (core).",
    "plddt": "AlphaFold's own confidence in each residue of a predicted model, 0-100. "
             "Above 90 = very reliable; below 50 = do not trust the shape there.",
    "alphafold": "a program that predicts a protein's 3D structure from its sequence. An "
                 "AlphaFold model is a computed prediction, not an experiment.",
    "chain": "one continuous protein (or DNA/RNA) molecule in a structure, labelled A, B, "
             "C ... Haemoglobin has four chains.",
    "residue": "one building block of a chain: one amino acid in a protein, one nucleotide "
               "in DNA/RNA. HIS57 = the histidine at position 57.",
    "ligand": "a molecule bound to the protein that is not part of its chain: a drug, a "
              "cofactor such as heme, a substrate, or a metal ion.",
    "active site": "the pocket of an enzyme where the chemical reaction happens, and the "
                   "few residues that carry it out.",
    "binding site": "the residues that hold a partner: a ligand, a metal ion, DNA or "
                    "another protein.",
    "domain": "a part of a chain that folds into its own compact unit and often has its "
              "own job. Large proteins are built from several domains.",
    "alpha-helix": "the backbone coiled like a spring (3.6 residues per turn), held by "
                   "hydrogen bonds between residue i and i+4. Drawn as a spiral.",
    "beta-strand": "a stretched-out piece of backbone. Several strands side by side, "
                   "joined by hydrogen bonds, make a beta-sheet. Drawn as a flat arrow.",
    "loop": "an irregular stretch linking helices and strands; often flexible and often "
            "where binding and active sites sit.",
    "secondary structure": "the local shape of the backbone: alpha-helix, beta-strand or "
                           "loop.",
    "fold": "the overall 3D arrangement of a chain's helices and strands.",
    "disulfide": "a covalent bond between the sulfur atoms of two cysteines; it staples "
                 "parts of the fold together.",
    "salt bridge": "a negatively charged side chain (Asp, Glu) next to a positively "
                   "charged one (Lys, Arg); opposite charges attract.",
    "hydrogen bond": "a weak attraction between a hydrogen on N or O and another N or O, "
                     "about 2.6-3.5 A between the heavy atoms. Many together hold a fold.",
    "hydrophobic": "water-avoiding: oily side chains (Leu, Ile, Val, Phe ...) that pack in "
                   "the protein core, away from water.",
    "ptm": "post-translational modification: a chemical group added after the protein is "
           "made, e.g. phosphate, sugar or acetyl.",
    "variant": "a position where the amino acid differs between people; a disease variant "
               "is one known to cause illness.",
    "mutation": "a change of one amino acid into another, written like R248Q: arginine "
                "(R) at position 248 becomes glutamine (Q).",
    "sifts": "the service that links each residue of a PDB structure to its position in "
             "the UniProt sequence, because the two are often numbered differently.",
    "cartoon": "PyMOL's ribbon drawing of the backbone: spirals for helices, arrows for "
               "strands, tubes for loops.",
    "selection": "a named set of atoms in PyMOL (e.g. 1PPB_active) that you can use in "
                 "any command: zoom 1PPB_active, color red, 1PPB_active.",
    "object": "one loaded structure in PyMOL, listed in the panel on the right (e.g. "
              "4HHB). Prospector also lets you use the name you loaded it with.",
    "x-ray": "X-ray crystallography: the structure comes from how a protein crystal "
             "scatters X-rays. The most common method; gives a resolution value.",
    "cryo-em": "cryo-electron microscopy: the structure is reconstructed from images of "
               "frozen molecules. Good for large complexes.",
    "nmr": "nuclear magnetic resonance: the structure of a protein in solution, given as "
           "a bundle of models; it has no resolution value.",
}
_GLOSSARY_ALIASES = {
    "colors": "colours", "color": "colours", "colour": "colours",
    "colour code": "colours", "color code": "colours",
    "helix": "alpha-helix", "alpha helix": "alpha-helix", "helices": "alpha-helix",
    "strand": "beta-strand", "beta strand": "beta-strand", "sheet": "beta-strand",
    "beta sheet": "beta-strand", "beta-sheet": "beta-strand", "bfactor": "b-factor",
    "b factor": "b-factor", "temperature factor": "b-factor", "rfree": "r-free",
    "r free": "r-free", "rwork": "r-work", "r work": "r-work", "h-bond": "hydrogen bond",
    "hbond": "hydrogen bond", "hydrogen bonds": "hydrogen bond", "disulphide": "disulfide",
    "disulfide bond": "disulfide", "salt bridges": "salt bridge", "em": "cryo-em",
    "cryoem": "cryo-em", "xray": "x-ray", "crystallography": "x-ray", "a": "angstrom",
    "angstroms": "angstrom", "ligands": "ligand", "chains": "chain", "residues": "residue",
    "domains": "domain", "variants": "variant", "mutations": "mutation",
    "hydrophobicity": "hydrophobic", "confidence": "plddt", "id": "pdb id",
    "post-translational modification": "ptm", "selections": "selection",
    "objects": "object", "loops": "loop", "amino acid": "residue", "ribbon": "cartoon",
}


def _glossary_key(text):
    key = " ".join(text.lower().replace("?", " ").split())
    key = _GLOSSARY_ALIASES.get(key, key)
    return key if key in GLOSSARY else None


def glossary(*args, **kwargs):
    """
DESCRIPTION

    Explain a term Prospector or PyMOL uses, in one or two sentences.

USAGE

    glossary             list the terms
    glossary <term>

EXAMPLES

    glossary RMSD
    glossary resolution
    """
    text = base._raw(args, kwargs).strip()

    def action():
        base._header("glossary %s" % text)
        _say("")
        if not text:
            _say("TERMS I CAN EXPLAIN   (type:  glossary <term>)")
            names = sorted(GLOSSARY)
            for start in range(0, len(names), 3):
                _say("   " + "".join("%-24s" % n for n in names[start:start + 3]).rstrip())
            _say("   ... and every amino acid:  glossary histidine")
            _say(report.LINE)
            return
        key = _glossary_key(text)
        if key is None and text.upper() in tools.AA:
            name, group, use = tools.AA[text.upper()]
            _say("%s = %s (%s): %s." % (text.upper(), name, group, use))
            _say(report.LINE)
            return
        if key is None:
            for code, (name, group, use) in tools.AA.items():
                if text.lower() in (name, name + "s"):
                    _say("%s (%s, %s): %s." % (name.upper(), code, group, use))
                    _say(report.LINE)
                    return
            close = difflib.get_close_matches(text.lower(), list(GLOSSARY)
                                              + list(_GLOSSARY_ALIASES), n=3, cutoff=0.6)
            _say("I have no entry for '%s'.%s" % (text, " Did you mean: %s?" % ", ".join(
                dict.fromkeys(_GLOSSARY_ALIASES.get(c, c) for c in close)) if close else ""))
            _say("Type  glossary  to see every term.")
            _say(report.LINE)
            return
        _say("%s" % key.upper())
        _say("   " + GLOSSARY[key])
        _say(report.LINE)

    worker.submit(lambda: None, lambda result, error: action())


# -------------------------------------------------------------------- start
def start(*args, **kwargs):
    """
DESCRIPTION

    A two-minute introduction for first-time users: the five things to type
    first, and how to move the picture with the mouse.

USAGE

    start
    """
    def action():
        base._header("start")
        _say("")
        _say("WELCOME      Prospector lets you work with protein structures by name.")
        _say("             Type a line below in the box marked  PyMOL>  and press Enter.")
        _say("")
        _say("TRY THESE, IN ORDER")
        _say("   1  grab haemoglobin        loads a structure by its name")
        _say("   2  about haemoglobin       what it is, and its key places")
        _say("   3  focus 1                 zooms to the first key place")
        _say("   4  tour haemoglobin        a guided walk; type  next  to continue")
        _say("   5  viewas sticks           draws it another way;  back  undoes it")
        _say("")
        _say("NOT SURE WHAT TO TYPE?   Say it in your own words after  ask :")
        _say("   ask show me where the drug binds")
        _say("   ask how do I colour by chain")
        _say("")
        _say("THE MOUSE    left-drag = rotate   |   scroll wheel = zoom")
        _say("             middle-drag (or Ctrl + left-drag) = move sideways")
        _say("             click an atom = PyMOL names it in this window")
        _say("")
        _say("GOOD TO KNOW")
        _say("   status      what is loaded and what you are looking at")
        _say("   back        undo the last change to the picture")
        _say("   glossary    what a word means (RMSD, resolution, chain ...)")
        _say("   info        the full list of commands")
        _say(report.LINE)

    worker.submit(lambda: None, lambda result, error: action())


# ------------------------------------------------------------------- status
def _shown(obj):
    words = []
    for rep, word in (("cartoon", "ribbon"), ("sticks", "sticks"), ("spheres", "balls"),
                      ("surface", "surface"), ("lines", "lines"), ("mesh", "mesh"),
                      ("dots", "dots"), ("ribbon", "backbone trace")):
        if cmd.count_atoms("%s and rep %s" % (obj, rep)):
            words.append(word)
    return " + ".join(words) or "nothing (hidden)"


def status(*args, **kwargs):
    """
DESCRIPTION

    Where am I?  Lists what is loaded, what each structure is, how it is
    drawn, the names you can use, and what you could do next.

USAGE

    status
    """
    def action():
        base._header("status")
        _say("")
        names = base._molecules()
        if not names:
            _say("NOTHING IS LOADED YET.")
            _say("Load a structure by name:   grab haemoglobin")
            _say("New here?                   start")
            _say(report.LINE)
            return
        enabled = base._molecules(enabled_only=True)
        _say("LOADED       %d structure%s" % (len(names), "" if len(names) == 1 else "s"))
        for name in names:
            identity = base._identify(name)
            what = ("PDB entry %s" % identity[1] if identity and identity[0] == "pdb"
                    else "AlphaFold model of %s" % identity[1] if identity
                    else "not linked to a database entry")
            aliases = sorted(a for a, o in base._aliases.items() if o == name
                             and a != name.lower() and "=" not in a and "<" not in a)
            chains = len([c for c in cmd.get_chains("%s and polymer" % name)]) or 1
            residues = cmd.count_atoms("%s and polymer and name CA" % name)
            _say("   %-14s %s%s" % (name, what, "" if name in enabled else "   [hidden]"))
            _say("   %-14s %d chain%s, %d residues; drawn as %s"
                 % ("", chains, "" if chains == 1 else "s", residues, _shown(name)))
            if aliases:
                _say("   %-14s you can also call it: %s" % ("", ", ".join(aliases[:4])))
        current = base._last[0] if base._last[0] in names else None
        if current and len(names) > 1:
            _say("")
            _say("CURRENT      %s  (commands without a name act on this one)" % current)
        extras = []
        if str(cmd.get("grid_mode")) not in ("0", "off"):
            extras.append("two-panel comparison is on (whatif off  leaves it)")
        from . import tier4
        if tier4._tour["scenes"] and tier4._tour["object"] in names:
            extras.append("a tour of %s is at stop %d of %d (next / prev)"
                          % (tier4._tour["object"], tier4._tour["index"] + 1,
                             len(tier4._tour["scenes"])))
        if tier4._quiz["entry"] is not None:
            extras.append("a quiz is running (quiz hint / quiz reveal)")
        if base._undo:
            extras.append("%d step%s can be undone with  back"
                          % (len(base._undo), "" if len(base._undo) == 1 else "s"))
        extras.append("teach mode is %s" % ("on" if teach.enabled() else
                                           "off (prospector teach on  shows the PyMOL "
                                           "commands)"))
        _say("")
        _say("RIGHT NOW    " + extras[0])
        for line in extras[1:]:
            _say("             " + line)
        target = current or names[-1]
        _say("")
        _say("YOU COULD    about %s  |  tour %s  |  focus  (key places)" % (target, target))
        _say("             viewas sticks  |  scout %s histidine  |  ask <anything>" % target)
        _say(report.LINE)

    worker.submit(lambda: None, lambda result, error: action())


# ---------------------------------------------------------------------- ask
_FILLER = re.compile(
    r"\b(please|kindly|can you|could you|would you|will you|i want to|i wanna|"
    r"i would like to|i'd like to|i need to|let me|let's|lets|for me|just|now|hey|hi|"
    r"hello|prospector|pymol|thanks|thank you)\b")
_STOP = set("the a an of in on for to at this that it its my me is are was be with by from "
            "and or enzyme protein proteins structure structures molecule molecules model object "
            "all every each some any one ones there here what which where how do does i "
            "you we see look up out into onto about called named have has had having there "
            "whats wheres please most important main major key biggest "
            "top crucial relevant really very actually exactly also too so if when can "
            "could would should will did am been being them they their these those".split())
# ordinary words that are never (part of) a protein's name
_JUNK = set("atom atoms bond bonds view picture thing things stuff way good better best "
            "nice cool great again reset clear delete mistake something anything work "
            "works working action happens happen moves move part parts commands command "
            "use used using confident sure made make makes made look looks like only "
            "just version kind type sort available exist exists same different more "
            "less much many big small large long short high low quality solved "
            "determined used method new old first last next previous screen whole "
            "entire full simple please help want need know think mean means right "
            "wrong here now then well ok okay yes no not any every between near around "
            "holds hold touches touch surrounds surround attached bound stuck linked "
            "important interesting famous example examples around center centre "
            "place differently label labels sit sits goes together".split())
_VERBS = ("show", "color", "colour", "zoom", "focus", "view", "make", "load", "fetch", "get",
          "open", "find", "highlight", "measure", "save", "take", "explain", "tell",
          "compare", "start", "give", "draw", "display", "go", "hide", "spot", "mutate",
          "change", "then", "list", "search", "grab", "label", "undo", "print")

_AA_WORDS = {}
for _code, (_name, _group, _use) in tools.AA.items():
    _one = [k for k, v in tools.ONE_TO_THREE.items() if v == _code][0]
    for _word in (_name, _name + "s", _code.lower()):
        _AA_WORDS[_word] = _one
_AA_PATTERN = "|".join(sorted(_AA_WORDS, key=len, reverse=True))


def _letter(word):
    word = word.lower()
    if word in _AA_WORDS:
        return _AA_WORDS[word]
    if len(word) == 1 and word.upper() in tools.ONE_TO_THREE:
        return word.upper()
    return None


# How sure the rules are.  `_subject` notes how many leftover words it had to
# call "the protein's name"; a long leftover means the rule matched one word
# of a sentence it did not really understand.
_doubt = [0]
_named = []      # what the rules took for names of things that are not loaded
# words that describe a structure instead of naming one
_VAGUE = set("mutant mutants mutated normal healthy cancer cancers tumour tumor tumours "
             "tumors wildtype wild-type version versions variant typical usual common "
             "another other sick diseased found seen known famous".split())
_MAX_NAME_WORDS = 4
# Function words of other languages: two of them and the sentence is not
# English, so the English rules must not guess at it.
_FOREIGN = set(
    "il lo gli della dello delle degli dei del sul sulla sulle nel nella nelle che "
    "con una uno questo questa sono mostrami mostra evidenzia colora dove quali "
    "fammi vedere anche molecola proteina "                      # Italian
    "el los las del una uno muestra donde esta este proteina tambien "   # Spanish
    "les des une dans avec montre moi quel quelle "               # French
    "der das und ein eine zeige mir wo ist nicht mit dem den "    # German
    "uma os nas nos onde mostre".split())                         # Portuguese


_FOREIGN_VERBS = set("mostrami mostra evidenzia colora fammi zooma ingrandisci carica "
                     "muestra muestrame resalta zeige montre mostre dove donde "
                     "proteina molecola".split()) - set(["carica"])


def _foreign(text):
    words = re.findall(r"[a-zà-ÿ]+", text.lower())
    return len([w for w in words if w in _FOREIGN]) >= 2 or any(
        w in _FOREIGN_VERBS for w in words)


def _subject(text, *remove):
    """What is left of a sentence once the words that carried the request are
    taken out: normally the protein's name."""
    for pattern in remove:
        text = re.sub(pattern, " ", text)
    words = [w for w in re.findall(r"[a-z0-9][a-z0-9\-\.']*", text) if w not in _STOP]
    phrase = " ".join(words)
    if not phrase:
        return ""
    known = base._known_object(phrase)
    if known:
        return phrase
    everything = words
    words = [w for w in words if w not in _JUNK]
    phrase = " ".join(words)
    if not phrase:
        return ""
    for size in range(len(words) - 1, 0, -1):          # a known name inside the phrase
        for start in range(0, len(words) - size + 1):
            piece = " ".join(words[start:start + size])
            if base._known_object(piece):
                # The rest of the sentence was ignored.  One stray word is
                # fine ("p53 again"); several mean the sentence says more
                # than this rule understood ("p53 normal and its most common
                # mutated version in breast cancer").
                ignored = [w for w in everything if w not in piece.split()]
                if len(ignored) >= 2:
                    _doubt[0] = 99
                return piece
    # nothing loaded goes by this name: note how much was left over
    _named.append(phrase)
    _doubt[0] = max(_doubt[0], len(words))
    if len(words) >= 2 and any(w in _PART_WORDS or w in _AA_WORDS for w in words):
        _doubt[0] = 99      # a part of a structure was left over inside a "name"
    if any(w in _VAGUE for w in everything):
        _doubt[0] = 99      # "the mutant found in cancer" describes, it does not name
    return phrase


def _specs(text):
    """Residue addresses in a sentence: 'H/57', 'residue 57', 'his 87 in chain a'."""
    found = []
    for match in re.finditer(r"\b([a-z0-9]{1,2})/(?:[a-z]{3})?(-?\d+[a-z]?)(?:/([a-z0-9']{1,4}))?\b",
                             text):
        spec = "%s/%s" % (match.group(1).upper(), match.group(2).upper())
        if match.group(3):
            spec += "/" + match.group(3).upper()
        found.append((match.start(), spec, match.group(0)))
    taken = [m[2] for m in found]
    pattern = (r"\b(?:residue|res|position|number|aa|amino acid|%s)\s*#?\s*(\d+[a-z]?)\b"
               r"(?:\s*(?:in|of|on|from)?\s*chain\s*([a-z0-9])\b)?" % _AA_PATTERN)
    for match in re.finditer(pattern, text):
        if any(match.group(0) in t or t in match.group(0) for t in taken):
            continue
        spec = match.group(1).upper()
        if match.group(2):
            spec = "%s/%s" % (match.group(2).upper(), spec)
        found.append((match.start(), spec, match.group(0)))
    for match in re.finditer(r"\b([acdefghiklmnpqrstvwy])(\d{1,4})\b", text):    # "r248"
        word = match.group(0)
        if word in _NAME_LIKE or base._known_object(word) or any(
                word in t or t in word for t in [m[2] for m in found]):
            continue
        found.append((match.start(), match.group(2), word))
    found.sort()
    return [(spec, raw) for _, spec, raw in found]


# protein names that look like "letter + number" and are not residues
from .typed import NAME_LIKE as _NAME_LIKE      # noqa: E402  (p53 is not proline 53)


_STYLE_WORDS = dict((name, name) for name in tier3.STYLE_NAMES)
_STYLE_WORDS.update(tier3._STYLE_ALIASES)
_STYLE_WORDS.update({"ball and stick": "ballstick", "balls and sticks": "ballstick",
                     "space filling": "spheres", "space-filling": "spheres"})
_SCOUT_PHRASES = sorted(tools._scout_words() - set(["motif", "sequence", "pattern", "end",
                                                    "start", "small", "tiny", "ends",
                                                    "inside", "outside", "core", "ion",
                                                    "contacts", "stable", "ordered"]),
                        key=lambda word: (-len(word), word))


# ---------------------------------------------------------- understanding
# Sentences are normalised before the rules look at them, so the rules do not
# depend on one exact wording:
#   1. plurals         "active sites" -> "active site", "helices" -> "helix"
#   2. typos           "histidne" -> "histidine"  (only words close to a known one)
#   3. other wordings  "catalytic residues" -> "active site", "paint" -> "color"
_SYNONYMS = [
    (r"\bam+in+o[- ]?ac+ids?\b|\baminoacids?\b", "amino acid"),
    (r"\b(n[- ]?term(inus|inal)?|beginning|start|initial) (residue|amino acid)\b",
     "first amino acid"),
    (r"\b(c[- ]?term(inus|inal)?|end|final) (residue|amino acid)\b", "last amino acid"),
    (r"\b(1st)\b", "first"),
    (r"\b(have a look at|take a look at|look at|pull up|check out|put up)\b", "show"),
    (r"\b(everything|all|more|something|anything) about\b", "about"),
    (r"\b(interface|contacts?) between the (chains|subunits|monomers)\b", "interface"),
    (r"\b(the |do the |does the )?(chains|subunits|monomers) (touch|meet|stick together|"
     r"interact|contact each other|bind each other|hold together|fit together)\b",
     "interface"),
    (r"\bh[- ]bonds?\b", "hydrogen bonds"),
    (r"\b(i made a mistake|oops|that was wrong|that is wrong|not what i wanted|"
     r"undo that|wrong one)\b", "undo"),
    (r"\b(i am|i'm|im) bored\b|\b(something|anything) (cool|fun|interesting|nice|new|"
     r"random)\b", "surprise"),
    (r"\b(what now|what next|what should i do( next| now)?|what can i do (now|next)|"
     r"i (dont|don't|do not) know what to do|i am (stuck|lost)|i'm (stuck|lost))\b",
     "status"),
    (r"\b(does|do|doing) (all )?the (real |actual )?(work|job)\b|\bwhere the action "
     r"(happens|is)\b|\bbusiness end\b|\bworking part\b", "active site"),
    (r"\b(moves? (the most|a lot|around|most)|wobbl\w+|wiggl\w+|most (mobile|flexible)|"
     r"disordered)\b", "flexible"),
    (r"\b((chain|protein|sequence) (starts?|begins?)|(start|beginning) of the (chain|"
     r"protein|sequence)|n[- ]?term(inus|inal)?)\b", "n terminus"),
    (r"\b((chain|protein|sequence) (ends?|stops?|finish(es)?)|end of the (chain|protein|"
     r"sequence)|c[- ]?term(inus|inal)?)\b", "c terminus"),
    (r"\b(prettier|pretty|beautiful|nicer|nice|textbook|simplif\w+|clean(er)?( it)? up|"
     r"cleaner|tid(y|ier)( it)?( up)?)\b", "publication"),
    (r"\b(what is|what molecules? (is|are)|what else is|which molecules? (is|are)|"
     r"is anything|anything) (bound|attached|stuck)( to)?\b", "show ligands"),
    (r"\bwhat ?s\b|\bwhat's\b", "what is"),
    (r"\bwhere ?s\b|\bwhere's\b", "where is"),
    (r"\bactive sites\b", "active site"),
    (r"\b(catalytic|enzymatic|reactive|reaction) (sites?|cent(er|re)s?|residues?|triads?|"
     r"machinery|pockets?|parts?|amino acids?)\b", "active site"),
    (r"\bwhere (the )?(reaction|catalysis|chemistry) (happens|occurs|takes place|"
     r"is done)\b", "active site"),
    (r"\b(functional|important|key|critical|essential|main) (sites?|residues?|parts?|"
     r"regions?|amino acids?)\b", "active site"),
    (r"\bbinding (pockets?|clefts?|grooves?|cavity|cavities)\b", "pocket"),
    (r"\b(bound|small) (molecules?|compounds?)\b", "ligand"),
    (r"\b(compounds?|medicines?|medications?|pharmaceuticals?)\b", "ligand"),
    (r"\b(disease|pathogenic|harmful|cancer|clinical|known|common|famous) "
     r"(mutations?|variants?)\b", "disease variant"),
    (r"\b(mutations? hot ?spots?|hot ?spots?)\b", "disease variant"),
    (r"\bmutations\b", "disease variant"),
    (r"\b(sulfur|sulphur|cystine|s-s|ss|cysteine) (bridges?|bonds?|links?|pairs?)\b",
     "disulfides"),
    (r"\bdisulphides?\b", "disulfides"),
    (r"\b(ionic|ion|charge|charged|electrostatic) (bonds?|pairs?|interactions?|bridges?)\b",
     "salt bridges"),
    (r"\b(spirals?|alpha helix|alpha helices)\b", "helices"),
    (r"\b(on the (outside|surface|exterior)|outside of|exterior|solvent[- ]exposed|"
     r"exposed|outer (surface|part))\b", "exposed"),
    (r"\b(on the inside|in the (core|interior|middle)|inside of|interior|hidden|"
     r"inner (part|core))\b", "buried"),
    (r"\b(oily|greasy|water[- ]avoiding|water[- ]hating|apolar|non-?polar|lipophilic)\b",
     "hydrophobic"),
    (r"\bwater[- ]loving\b", "polar"),
    (r"\b(sugar chains?|oligosaccharides?|carbohydrates?)\b", "sugars"),
    (r"\b(magnify|enlarge|closer look at|close[- ]up (of|on)|cent(er|re) on|jump to|"
     r"bring me to|take me to|navigate to|move to|go to|get closer to|zoom in on|"
     r"zoom into|zoom in to|zoom on|focus on)\b", "zoom to"),
    (r"\b(reveal|present|visuali[sz]e|point out|point to|indicate|bring up|let me see|"
     r"i want to see|i wanna see|i would like to see|i'd like to see|can i see|"
     r"could i see|display)\b", "show"),
    (r"\b(paint|shade|tint)\b", "color"),
    (r"\b(emphasi[sz]e|mark|flag|pick out|light up)\b", "highlight"),
    (r"\b(pic|pics|photo|image|screen ?shot|figure file)\b", "picture"),
    (r"\b(amino acids?|aa|residues)\b(?=\s*#?\s*\d)", "residue"),
    (r"\b(resembles?|resembling|alike|homologous)\b", "similar"),
    (r"\b(contrast)\b", "compare"),
    (r"\b(whereabouts of|location of|position of|locate)\b", "where is"),
    (r"\bwhich (residues?|amino acids?|parts?|atoms?) (make up|form|are in|belong to|"
     r"is|are)\b", "show"),
    (r"\b(stands? for)\b", "mean"),
]
_SINGULAR = {"helices": "helix", "termini": "termini", "ions": "ions", "bases": "bases",
             "this": "this", "is": "is", "as": "as", "has": "has", "does": "does",
             "was": "was", "its": "its", "class": "class", "glycosylation sites":
             "glycosylation sites"}
_vocabulary = None


def _known_words():
    """Every word the rules can react to (used for plurals and typo repair)."""
    global _vocabulary
    if _vocabulary is None:
        words = set()
        sources = list(_SCOUT_PHRASES) + list(_STYLE_WORDS) + list(GLOSSARY) + \
            list(_GLOSSARY_ALIASES) + list(_AA_WORDS) + list(_VERBS)
        sources += ("active site binding pocket ligand drug inhibitor cofactor substrate "
                    "domain variant disease chain residue sequence similar tour quiz "
                    "surprise report session picture distance measure compare superpose "
                    "align mutate change replace explain describe highlight zoom focus "
                    "color colour rainbow hydrophobicity flexibility publication figure "
                    "structure protein feature annotation motif alphafold predicted "
                    "undo status glossary whole everything interaction contact "
                    "bridge bond site helix strand loop metal sugar glycan water "
                    "surface buried interface modified between").split()
        for source in sources:
            words.update(w for w in re.findall(r"[a-z][a-z\-]+", source.lower())
                         if len(w) >= 3)
        _vocabulary = words
    return _vocabulary


def _normalize(text):
    """Bring a sentence into the wording the rules expect (see above)."""
    known = _known_words()
    out = []
    for token in re.findall(r"\S+", text.lower().replace("?", " ").replace("!", " ")
                            .replace(",", " ")):
        word = token.strip(".;:\"()")
        if not word:
            continue
        plain = word.isalpha() and not base._known_object(word)
        if plain and word not in known:
            # 1. plural -> singular, when the singular is a word we know
            if word in _SINGULAR:
                word = _SINGULAR[word]
            elif word.endswith("ies") and word[:-3] + "y" in known:
                word = word[:-3] + "y"
            elif word.endswith("es") and word[:-2] in known:
                word = word[:-2]
            elif word.endswith("s") and word[:-1] in known:
                word = word[:-1]
            # 2. a likely typo of a word we know (long words only, close matches only)
            elif len(word) >= 6:
                close = difflib.get_close_matches(word, known, n=1, cutoff=0.86)
                if close:
                    word = close[0]
        out.append(word)
    sentence = " " + " ".join(out) + " "
    for pattern, replacement in _SYNONYMS:       # 3. other ways of saying it
        sentence = re.sub(pattern, replacement, sentence)
    return " ".join(sentence.split())


# ------------------------------------------------------------- being sure
# After a rule has picked a command, every meaningful word of the sentence
# must be accounted for: it is in the command itself, or it is one of the
# words that belong to that kind of request.  A word nobody accounted for
# means the sentence says more than the rule understood - then the rules step
# aside (the AI is asked, or the request is refused) instead of guessing.
_COVER_ANY = set(
    "show display see view look give tell make let want need please can could would "
    "protein proteins structure structures molecule molecules residue residues amino "
    "acid acids part parts thing things whole entire now then just also loaded current "
    "present screen picture chain chains".split())
_C = {
    "grab": "load fetch get open grab download bring pull draw best good version from "
            "kind small compound looks like what does bind binds binding contain "
            "contains containing bound has have with inside attached complex which "
            "any enzyme something that",
    "find": "find search list best good available exist exists version entry entries pdb "
            "options candidates from high resolution cryo-em cryo em x-ray xray crystal "
            "nmr bound containing",
    "molecule": "load fetch get open grab draw small compound chemical called named",
    "alphafold": "alphafold predicted predict prediction model load fetch get open shape "
                 "fold confident sure how",
    "about": "what who about describe information info details function role job purpose "
             "does do why important disease diseases linked associated related cause "
             "causes resolution method technique experiment organism species come comes "
             "size big large long heavy many much chain chains atoms solved determined "
             "discovered published when year old quality good reliable accurate card "
             "facts explain more kind type enzyme dimer monomer trimer tetramer oligomer "
             "kinase receptor membrane hormone complex antibody homodimer made weight "
             "mass length x-ray xray crystal subunits domains has have an for everything "
             "most main key from "
             "variant variants mutation mutations used use are there",
    "tour": "tour walk through guide guided around how work works function start take "
            "stop end quit",
    "focus": "zoom focus go take closer where located location site sites active binding "
             "ligand ligands drug inhibitor cofactor substrate domain domains variant "
             "variants disease mutation mutations chain whole everything reset view full "
             "center centre recenter out all highlight point known positions again "
             "which what does do have has in on to label sticks goes sits called only "
             "fit hotspot cancer cause find as load fetch get open grab most important "
             "main major key common",
    "explain": "explain what which why about describe role near neighbours neighbors "
               "around touch touches contacts interacts kind type surface buried exposed "
               "inside outside important does do on is find",
    "scout": "highlight find spot locate mark where list count how many which what are "
             "there have has any bonds bond bridges bridge every each called names in "
             "on is does do touch touches",
    "pocket": "pocket bind binds binding what where touch touches hold holds holding "
              "surround surrounds around near contacts interactions interaction "
              "hydrogen bonds ligand drug inhibitor substrate cofactor place neighbours "
              "neighbors coordinate coordinates which does do how doing in is the list "
              "go goes",
    "whatif": "change changes changed mutate mutated mutating mutation mutant replace "
              "replaced swap substitute what if whatif happens happen effect becomes "
              "become became to into with for side by compare model simulate try test "
              "why cause causes cancer disease first last removed remove does do was "
              "were is be would when one position its chain did original again only "
              "back panel panels single two revert undo leave close exit stop view "
              "comparison",
    "measure": "distance measure how far close between apart from to and long is are",
    "compare": "compare superpose superimpose align overlay put difference differences "
               "between and with to vs versus onto on top how different similar same "
               "rmsd is are as identical this it close related compared like",
    "similar": "similar related relatives homolog homologs homologues like resemble "
               "others other find does have what are to which",
    "viewas": "as in with draw render switch change turn representation style mode back "
              "to like only just atoms atom bonds bond backbone main chain trace "
              "skeleton every each all individual single chemical covalent by",
    "look": "color colour colours colors by style scheme each chain chains differently "
            "per rainbow hydrophobicity hydrophobic polarity secondary shape fold helix "
            "helices strands sheets element elements atom type types kind class property "
            "properties chemistry chemical publication figure paper pretty nice clean "
            "presentation poster b-factor b factor bfactor flexibility flexible mobility "
            "putty surface ligand buried how different a n c terminus start end to like "
            "textbook simple simplify tidy",
    "sequence": "sequence print letters order primary get what",
    "snap": "save take export create picture image png screenshot as file for my report "
            "high resolution called named to",
    "report": "report hand summary page write up one generate create about for my",
    "session": "save work session everything progress open load restore keep store as "
               "my all it called named saved list",
    "back": "undo back go take revert that last step steps once twice two three four "
            "five times changes things previous",
    "clear": "clear delete remove wipe close empty everything all start over again "
             "scratch new session page window text log output terminal messages clean",
    "status": "what where am which have is loaded looking status on",
    "start": "start started getting begin how new beginner tutorial introduction zoom "
             "move pan drag click select mouse controls use",
    "info": "what can do type ask say commands command list help info your",
    "quiz": "quiz test my knowledge challenge",
    "surprise": "surprise random",
    "next": "next continue stop carry on", "prev": "previous prev stop back one go",
    "more": "more full text explanation answer version details rest long say explain",
    "prospector": "teach pymol commands real raw actual behind used use what less too "
                  "much shorter answers answer output brief short mode words talking "
                  "long wordy longer full verbose again text give",
    "build": "build make create generate construct model assemble draw show see load get "
             "view display a an of with the from sequence double helix helices stranded dna "
             "rna alpha beta sheet sheets strand strands peptide peptides polypeptide lipid "
             "lipids bilayer bilayers membrane membranes layer layers micelle micelles "
             "assembly biological unit whole complete full complex cage capsid virus shell "
             "oligomer chains subunits copies all percent per cent nm nanometre "
             "nanometres nanometers form some me containing made molecules big small patch",
    "glossary": "what is are does mean means define meaning term explain stand for "
                "colour colours color colors why thing bit ball sphere stick red orange "
                "yellow blue green grey gray white",
    "color": "color colour make turn paint highlight helix helices strand strands sheet "
             "sheets loop loops coil coils ligand ligands drug inhibitor cofactor "
             "substrate water waters everything all whole chain residues",
    "hide": "hide remove turn off get rid without no water waters hydrogens hydrogen "
            "ligand ligands drug inhibitor cofactor substrate chain except but apart "
            "other than balls ribbon ribbons name tags names dashed lines bonds "
            "distances measurements measurement",
    "set": "transparent see-through translucent ghost surface ribbon thicker fatter "
           "bigger larger thinner smaller finer sticks stick bonds bond spheres sphere "
           "balls ball atoms atom",
    "rock": "spin rotate rotating rock around stop end quit rocking spinning moving "
            "turning turn",
    "turn": "turn rotate spin flip degrees degree deg around about along axis over",
    "move": "zoom closer bigger enlarge out further away smaller bit little slightly "
            "more",
    "bg_color": "background change",
    "distance": "hydrogen bonds bond draw",
    "show": "side chain chains sidechain sidechains draw as in with surface cartoon "
            "ribbon ribbons sticks stick lines line spheres sphere balls ball mesh dots "
            "dot view",
    "zoom": "zoom focus go to", "delete": "delete remove unload close measurement "
            "measurements distance distances dashed line lines last",
}
_COVER = dict((head, set(words.split())) for head, words in _C.items())
_COVER["mutate"] = _COVER["whatif"]
for _head in ("whatif", "mutate", "explain", "focus", "measure", "scout", "color"):
    _COVER[_head] = _COVER[_head] | set(_AA_WORDS)
_COVER["viewas"] = _COVER["viewas"] | set(
    w for phrase in _STYLE_WORDS for w in phrase.replace("-", " ").split())

# adjectives that change what is asked for ("the biggest pocket"): harmless in
# a protein's name, but they must be accounted for like any other word
_EMPHASIS = set("most important main major key biggest top crucial relevant".split())
# words that are never the name of a protein or a molecule
_NOT_A_NAME = set(
    "conserved conservation electrostatics electrostatic charge charges text original "
    "simulation dynamics dock docking blast background steroid steroids lipid lipids "
    "sugar sugars hormone hormones vitamin drug drugs something anything nothing "
    "everything stuff example picture image color colour bigger smaller thicker "
    "thinner faster slower better worse again more less details answer answers "
    "question help mistake wrong right screen window page view views label labels name "
    "names tag tags number numbers list table file files folder session work homework "
    "sandwich degrees angle angles rotation animation movie video interesting "
    "important common rare first last next previous new old same different similar "
    "side chains chain atoms atom bonds bond hydrogens hydrogen waters water surface "
    "region regions area areas place places spot zone section sections segment middle "
    "end ends bottom left front biggest largest smallest "
    "shape shapes size big small energy force stability stable interaction "
    "interactions properties property elements structure protein molecule".split())


def _clean(text):
    """The sentence as the rules see it: normalised, filler words removed."""
    return " " + " ".join(_FILLER.sub(" ", _normalize(text)).split()) + " "


def _unaccounted(command, text):
    """The meaningful words of ``text`` that the chosen command does not explain."""
    head = command.split()[0]
    if head == "choose":             # everything after "select" is handed over
        return []
    in_command = set(re.findall(r"[a-z0-9]+", command.lower()))
    allowed = in_command | _COVER_ANY | _COVER.get(head, set())
    blob = command.lower()
    left = []
    cleaned = _clean(text)
    # a chain the sentence names must be in the command ("... of chain A")
    for letter in re.findall(r"\bchains? ([a-z0-9])\b", cleaned):
        if not re.search(r"\bchain %s\b|\b%s/|chain=%s" % (letter, letter, letter), blob):
            left.append("chain %s" % letter.upper())
    for word in re.findall(r"[a-z0-9][a-z0-9\-'\+]*", cleaned):
        if word in allowed or (word in _STOP and word not in _EMPHASIS):
            continue
        pieces = re.findall(r"[a-z0-9]+", word)
        if pieces and word not in _EMPHASIS and all(
                p in allowed or p in _STOP or p in ("m", "s", "d", "ll", "re",
                                                              "t", "ve")
                          for p in pieces):                    # i'm, what's, don't
            continue
        if word.endswith("s") and word[:-1] in allowed or word + "s" in allowed:
            continue
        if word.isdigit() and word in blob:          # 248 inside R248Q
            continue
        if base._known_object(word):                 # what the user calls an object
            continue
        if word in _PART_WORDS and (_inside_current(word) or "").lower() in in_command:
            continue                                 # "the zinc" -> ZN
        if _molecule_code(word) and _molecule_code(word).lower() in in_command:
            continue                                 # cholesterol -> ligand=CLR
        number = re.match(r"^[a-z]{1,3}(\d+)[a-z]?$", word)       # arg248, r248
        if number and number.group(1) in blob:
            continue
        if len(word) == 1 and head in ("whatif", "mutate") and word in blob:
            continue
        if len(word) >= 3 and word in blob.replace(" ", ""):
            continue
        left.append(word)
    return left


def _molecule_code(word):
    from ..core import chemcomp
    return chemcomp.COMMON.get(word, "")


def _namey(word, text):
    """Is there positive evidence that this one word is (part of) a name?"""
    from ..core import fuzzy
    if not re.search(r"[a-z]", word):
        return False                                 # a bare number is not a name
    if re.search(r"\d", word) and re.search(r"[a-z]", word):
        return True                                  # p53, 4hhb, cdk2
    if len(word) >= 2 and re.search(r"[a-z]", word) and re.search(
            r"\b%s\b" % re.escape(word.upper()), text):
        return True                                  # typed in capitals: ATP, GFP
    system = _known_words() | _COVER_ANY | set().union(*_COVER.values())
    if word in system or word in _NOT_A_NAME or word in _STOP:
        return False
    normal = fuzzy.normalize_spelling(word)[0].lower()
    return fuzzy.known(normal) or bool(_molecule_code(word)) or len(word) <= 6


def _name_like(subject, text, explicit):
    """Could ``subject`` be the name of a protein or a molecule?  ``explicit``:
    the user said load / fetch / open, so a single odd word is taken as a name."""
    words = subject.split()
    if not words or len(words) > 5:
        return False
    if base._known_object(subject):
        return True
    if any(w in _NOT_A_NAME for w in words):
        return False
    if any(_namey(w, text) for w in words):
        return True
    return explicit and len(words) == 1          # "load heme"; not "show side chains"


def _interpret(text):
    """One request -> (command line, plain description) or (None, reason).
    The rules answer only when they are sure: a sentence in another language,
    or one where most words were left over as a "name", is not guessed at."""
    if _foreign(text):
        return None, ("'I only read English without the AI.' Switch the AI on "
                      "(prospector ai) to ask in your own language.")
    pointing = re.search(r"\b(this|that|the selected|the clicked|the chosen|selected|"
                         r"clicked) (residue|amino acid|one)\b|\bthe residue i "
                         r"(clicked|selected|picked)\b", text.lower())
    if pointing and not re.search(r"\d", text):
        clicked = base._clicked_residue()
        if clicked is None:
            return None, ("'This residue' - which one? Click a residue in the 3D window "
                          "first, or give its number (e.g. residue 248).")
        obj, chain, resi, resn = clicked
        where = "residue %s%s" % (resi, " in chain %s" % chain if chain else "")
        text = text[:pointing.start()] + where + text[pointing.end():]
    _doubt[0] = 0
    del _named[:]
    command, description = _rules(text)
    if command is not None and command.split()[0] not in ("choose", "glossary"):
        explicit = bool(re.search(r"\b(load|fetch|open|grab|download|find|search|get)\b",
                                  text.lower()))
        head = command.split()[0]
        arguments = _COVER.get(head, set()) | (tools._scout_words() if head == "scout"
                                               else set())
        for name in _named:
            # only names that really ended up in the command matter
            if not name or name not in command.lower():
                continue
            if all(w in arguments for w in name.split()):
                continue             # a second thing to scout, "stop" for the tour ...
            if not _name_like(name, text, explicit):
                return None, None
    if command is not None and _doubt[0] > _MAX_NAME_WORDS:
        return None, None
    if command is not None and _unaccounted(command, text):
        return None, None
    return command, description


def _rules(text):
    """The rule chain itself: one request -> (command, description) or (None, None)."""
    normal = _normalize(text)
    t = " " + " ".join(_FILLER.sub(" ", normal).split()) + " "

    raw = " " + " ".join(text.lower().replace("?", " ").split()) + " "

    def has(pattern):
        return re.search(pattern, t)

    if re.search(r"\b(what can (you|i|we) (do|type|ask|say)|what do you do|"
                 r"what are (the|your) commands)\b", raw):
        return "info", "list every command"
    if re.search(r"\b(teach|pymol commands|real commands|raw commands|commands behind|"
                 r"what commands|actual commands|commands (you|u) used?)\b", raw):
        return "prospector teach on", "show the PyMOL commands behind every action"

    match = re.match(r"^ (?:select|choose|pick|make a selection of|selection of) (.+) $", raw)
    if match:
        return "choose %s" % match.group(1).strip(), "select that and show PyMOL's own line"
    built = _build_request(text, raw)
    if built:
        return built

    # --- moving through things -------------------------------------------------
    if has(r"^ (next|continue|go on|next stop|carry on) $"):
        return "next", "go to the next stop of the tour"
    if has(r"^ (previous|prev|previous stop|go back one stop) $"):
        return "prev", "go back one stop of the tour"
    if has(r"\b(revert|undo|remove|leave|close|exit|stop)\b.*\b(mutations?|mutant|two "
           r"panels?|side by side|comparison view|what ?if)\b") or has(
            r"\b(original|one panel|single panel) (again|only|back)\b"):
        return "whatif off", "leave the two-panel view and show the original again"
    if has(r"\b(undo|go back|take that back|revert)\b") or t.strip() == "back":
        count = has(r"\b(twice|two|three|four|five|[2-5])\b")
        if count:
            steps = {"twice": 2, "two": 2, "three": 3, "four": 4, "five": 5}.get(
                count.group(1)) or int(count.group(1))
            return "back %d" % steps, "undo the last %d steps" % steps
        return "back", "undo the last change to the picture"
    if has(r"\b(where am i|what is loaded|what's loaded|what have i loaded|status|"
           r"what am i looking at|what do i have|which protein is this|what protein is "
           r"this|what is this|what is on (the )?screen)\b"):
        return "status", "show what is loaded and how it is drawn"
    if has(r"\b(what can (you|i) do|what can i type|list (the )?commands|all commands|"
           r"show (the )?commands)\b") or t.strip() in ("help", "commands", "info"):
        return "info", "list every command"
    if has(r"\b(get started|getting started|how do i start|where do i start|"
           r"i am new|i'm new|beginner|tutorial|introduction)\b") or t.strip() == "start":
        return "start", "show the first steps for a new user"

    # --- vocabulary -------------------------------------------------------------
    match = has(r"\b(?:what is|what's|what are|what does|define|meaning of|explain the "
                r"term|explain what)\s+(?:an?\s+|the\s+)?(.+?)(?:\s+mean|\s+means|\s+is)? $")
    if not match:
        match = has(r"^ (?:what is|what does|define)\s+(?:an?\s+|the\s+)?(.+?)"
                    r"(?:\s+mean|\s+means)? $")
    if match and _glossary_key(match.group(1)):
        key = _glossary_key(match.group(1))
        the_one_here = key in ("active site", "binding site", "ligand", "domain",
                               "variant", "chain", "resolution", "rmsd") and re.search(
            r"\b(what is|what's|what are)\s+the\b", t)
        if not the_one_here:        # "what is AN active site" = the meaning of the word
            return "glossary %s" % key, "explain the term '%s'" % match.group(1)

    # --- changing residues ------------------------------------------------------
    tokens = [w for w in text.replace(",", " ").split() if tools._MUTATION.match(w)
              and re.search(r"\d", w) and len(w) >= 3]
    change = re.search(r"\b(%s|[a-z])\s*(\d+)\s*(?:to|into|with|->|for|by|becomes?|became|as|was|"
                       r"were|is|turns into)\s*"
                       r"(?:an?\s+)?(%s|[a-z])\b" % (_AA_PATTERN, _AA_PATTERN), t)
    if re.search(r"\b(first|last)\s+(residue|amino acid|one|position)\b", t) or (
            change and base._known_object(change.group(0).split()[0])):
        change = None
    if not tokens and change and _letter(change.group(1)) and _letter(change.group(3)) \
            and has(r"\b(change|chang|mutat|replac|swap|substitut|what if|whatif|turn|"
                    r"becomes|convert)"):
        chain = re.search(r"\bchain\s*([a-z0-9])\b", t)
        token = "%s%s%s" % (_letter(change.group(1)), change.group(2),
                            _letter(change.group(3)))
        tokens = ["%s/%s" % (chain.group(1).upper(), token) if chain else token]
        t = t.replace(change.group(0), " ")
        if chain:
            t = t.replace(chain.group(0), " ")
    if tokens and (has(r"\b(change|chang|mutat|replac|swap|substitut|what if|whatif|turn|"
                       r"happen|effect|becomes|convert|side by side)") or True):
        cleaned = t
        for token in tokens:
            cleaned = cleaned.replace(token.lower(), " ")
        subject = _subject(cleaned, r"\b(change|changed|changing|mutate|mutation|mutated|"
                           r"mutant|replace|replaced|swap|substitute|what if|whatif|"
                           r"happens|happen|effect|effects|show|side by side|if|were|was|"
                           r"compare|would|when|make|model|why|does|do|causes?|caused|"
                           r"cancer|disease|bad|harmful|dangerous|simulate|became|"
                           r"become|try|test|see|what)\b")
        return ("whatif %s" % " ".join(([subject] if subject else []) + tokens),
                "change %s%s and compare original and changed version side by side"
                % (", ".join(tokens), " in %s" % subject if subject else ""))

    # "what happens if p53 changes its first amino acid", "mutate residue 248"
    verb = has(r"\b(change|changes|changed|changing|mutat\w*|replac\w*|swap\w*|"
               r"substitut\w*|alter\w*|modif\w*|remov\w*|delet\w*|los[et]\w*|knock\w*|"
               r"what if|whatif|different|broken|breaks?|damag\w*)\b")
    where = re.search(r"\b(first|last)\s+(?:residue|amino acid|one|position)\b", t) or \
        re.search(r"\b(?:residue|amino acid|position|number|%s)\s*#?\s*(-?\d+)\b"
                  % _AA_PATTERN, t)
    if verb and where and not has(r"\b(colou?r|view|style|look)\b"):
        target = re.search(r"\b(?:to|into|with|for|by|becomes?|as|was|were|is|be)\s+"
                           r"(?:an?\s+)?(%s)\b"
                           % _AA_PATTERN, t)
        chain = re.search(r"\bchain\s*([a-z0-9])\b", t)
        token = where.group(1) + (_letter(target.group(1)) if target else "")
        if chain:
            token = "%s/%s" % (chain.group(1).upper(), token)
        cleaned = t.replace(where.group(0), " ")
        for found in (target, chain):
            if found:
                cleaned = cleaned.replace(found.group(0), " ")
        subject = _subject(cleaned, verb.re.pattern, r"\b(happens?|effects?|show|side by "
                           r"side|were|was|compare|make|model|its|his|her|residue|amino "
                           r"acid|gets?|got|is|one|position)\b")
        return ("whatif %s" % " ".join(([subject] if subject else []) + [token]),
                "change the %s residue%s and compare original and changed version side "
                "by side" % (where.group(1) if not where.group(1).isdigit()
                             else "number-" + where.group(1),
                             " of %s" % subject if subject else ""))

    specs = _specs(t)

    here = base._last[0] if base._last[0] in base._molecules() else ""
    colours = (r"red|blue|green|yellow|orange|purple|pink|white|grey|gray|cyan|magenta|"
               r"black|salmon|violet|brown")
    numbers = {"once": 1, "one": 1, "twice": 2, "two": 2, "three": 3, "four": 4,
               "five": 5, "2": 2, "3": 3, "4": 4, "5": 5, "1": 1}

    # --- how much text ----------------------------------------------------------
    if has(r"\b(less text|too much text|shorter (answers?|output|text)|be brief|short "
           r"(answers?|mode|version)|less (output|words|talking)|too (long|wordy))\b"):
        return "prospector short on", "print only the answer from now on (more = the rest)"
    if has(r"\b(full (text|explanation|answer|version)|more details?|tell more|"
           r"explain more|show more|the rest|long version|say more)\b") or \
            t.strip() in ("more", "tell me more", "details"):
        return "more", "show the full text of the last command"
    if has(r"\b(long answers?|full answers? again|longer (answers?|text)|verbose)\b"):
        return "prospector short off", "print the full explanation again"

    # --- several steps back -----------------------------------------------------
    match = has(r"\b(?:undo|back|go back)\b.*?\b(once|twice|one|two|three|four|five|[1-5])"
                r"\b(?: (?:steps?|times?|things?|changes?|commands?))?")
    if match and not specs:
        count = numbers[match.group(1)]
        return ("back %d" % count if count > 1 else "back",
                "undo the last %d steps" % count if count > 1 else "undo the last step")
    if has(r"\b(revert|undo|remove|leave|close|exit|stop)\b.*\b(mutations?|mutant|two "
           r"panels?|side by side|comparison view|what ?if)\b") or has(
            r"\b(original|one panel|single panel) (again|only|back)\b"):
        return "whatif off", "leave the two-panel view and show the original again"

    # --- small molecules --------------------------------------------------------
    match = has(r"\b(?:protein|proteins|structure|structures|enzyme|enzymes|something)\b"
                r".*?\b(?:binds?|binding|contains?|containing|has|have|with|bound to|"
                r"in complex with)\s+(?:the |a |an )?([a-z0-9\-\+]+(?: [a-z]+)?)"
                r"(?: (?:bound|inside|in it|attached))? $")
    if match and has(r"^ (a |an |any |which |what |find |show |load |get |give |grab )?"
                     r"(protein|proteins|structure|structures|enzyme|something)\b"):
        from ..core import chemcomp
        thing = match.group(1).strip()
        code = chemcomp.COMMON.get(thing) or chemcomp.COMMON.get(thing.split()[0])
        if code is None and re.match(r"^[a-z0-9]{2,3}$", thing) and thing.upper() in text:
            code = thing.upper()
        if code:
            return ("grab protein ligand=%s" % code,
                    "load a protein that has %s (%s) bound" % (thing, code))
        return None, None
    match = has(r"^ (?:show |load |get |draw |see |give |display |grab |open )?(?:the |a )?"
                r"(?:small )?(?:molecule|compound|chemical)(?: of| called| named)? "
                r"([a-z0-9][a-z0-9\-\+' ]{1,40}?) $") or has(
        r"^ (?:show |load |get |draw |see |give |display |grab |open )?(?:the |a )?"
        r"([a-z0-9][a-z0-9\-\+']{1,30}(?: acid)?) (?:molecule|compound) $") or has(
        r"^ what does (?:the |a |an )?([a-z0-9][a-z0-9\-\+']{1,30}(?: acid)?) "
        r"(?:molecule )?looks? like $")
    if match and not base._known_object(match.group(1).strip()):
        name = match.group(1).strip()
        if " looks" in t or " look " in t:
            return "grab %s" % name, "find and load %s" % name
        return "molecule %s" % name, "load the small molecule %s" % name

    # --- background, thickness, turning ------------------------------------------
    match = has(r"\bbackground\b.*?\b(%s)\b" % colours) or has(
        r"\b(%s)\b.*?\bbackground\b" % colours)
    if match:
        colour = match.group(1).replace("gray", "grey")
        return "bg_color %s" % colour, "make the background %s" % colour
    match = has(r"\b(sticks?|bonds?|spheres?|balls?|atoms?)\b.*?\b(thicker|fatter|bigger|"
                r"larger|thinner|smaller|finer)\b")
    if match:
        more = match.group(2) in ("thicker", "fatter", "bigger", "larger")
        if match.group(1).startswith(("stick", "bond")):
            return ("set stick_radius, %s" % ("0.4" if more else "0.12"),
                    "make the sticks %s (normal: set stick_radius, 0.25)"
                    % ("thicker" if more else "thinner"))
        return ("set sphere_scale, %s" % ("1.3" if more else "0.5"),
                "make the balls %s (normal: set sphere_scale, 1)"
                % ("bigger" if more else "smaller"))
    if has(r"\b(stop|no more|end|quit)\b.*\b(rotat\w*|spin\w*|rock\w*|moving|turning)\b"):
        return "rock 0", "stop the rocking"
    match = has(r"\b(?:turn|rotate|spin|flip)\b.*?\b(\d{1,3})\s*(?:degrees?|deg|°)?"
                r"(?:.*?\b(?:around|about|on|along)\s+(?:the )?([xyz])\b)?")
    if match and not specs and has(r"\b(degrees?|deg|around|about|turn|rotate)\b"):
        axis = match.group(2) or "y"
        return ("turn %s, %s" % (axis, match.group(1)),
                "turn the view %s degrees around the %s axis" % (match.group(1), axis))
    if has(r"^ (flip|turn) (it|this|the (protein|structure|molecule))?( around| over)? ?$"):
        return "turn y, 180", "turn the view half-way round (180 degrees)"
    if has(r"\b(zoom in|closer|make (it|this|the \w+) bigger|bigger view|enlarge)\b") and \
            not has(r"\b(to|on|at|into)\b") and len(t.split()) <= 5:
        return "move z, 15", "move the camera closer (scroll the mouse wheel does the same)"
    if has(r"\b(zoom out|further away|make (it|this) smaller)\b") and has(
            r"\b(bit|little|slightly|more)\b") and len(t.split()) <= 6:
        return "move z, -15", "move the camera back a little"

    # --- hiding things -----------------------------------------------------------
    hide = has(r"\b(hide|remove|turn off|get rid of|without|no)\b")
    if hide:
        reps = [(r"\b(sticks?)\b", "sticks"), (r"\b(lines?)\b", "lines"),
                (r"\b(spheres?|balls?)\b", "spheres"),
                (r"\b(ribbons?|cartoon)\b", "cartoon"), (r"\b(surface)\b", "surface"),
                (r"\b(mesh)\b", "mesh"), (r"\b(dots?)\b", "dots"),
                (r"\b(labels?|name tags?|names|tags)\b", "labels"),
                (r"\b(dashes|dashed lines?|hydrogen bonds?|distances?|measurements?)\b",
                 "dashes")]
        chain = has(r"\bchain\s+([a-z0-9])\b")
        if has(r"\b(everything|all)\b.*\b(except|but|apart from|other than)\b"):
            keep = None
            if chain:
                keep = ("chain %s" % chain.group(1).upper(), "chain %s"
                        % chain.group(1).upper())
            elif has(r"\b(ligands?|drug|inhibitor|cofactor|substrate)\b"):
                keep = ("organic", "the bound molecules")
            elif has(r"\b(protein|chains?)\b"):
                keep = ("polymer", "the protein")
            if keep:
                return ("hide everything, not %s" % keep[0],
                        "hide everything except %s" % keep[1])
            return None, None
        if chain and len(t.split()) <= 5:
            return ("hide everything, chain %s" % chain.group(1).upper(),
                    "hide chain %s (show it again: show cartoon, chain %s)"
                    % (chain.group(1).upper(), chain.group(1).upper()))
        if has(r"\bhydrogens?\b") and not has(r"\bbonds?\b"):
            return "hide everything, hydro", "hide the hydrogen atoms"
        if has(r"\b(ligands?|drug|inhibitor|cofactor|substrate)\b") and len(t.split()) <= 5:
            return "hide everything, organic", "hide the bound molecules"
        for pattern, name in reps:
            if has(pattern) and len(t.split()) <= 6:
                return "hide %s" % name, "hide the %s" % name
    if has(r"^ (show|display|see|draw) (the |all )?(side ?chains?) $"):
        return ("show sticks, polymer and not backbone",
                "show every side chain as sticks (undo: hide sticks)")
    match = has(r"^ (?:show|display|draw|make|view) chain ([a-z0-9]) (?:as|in|with) "
                r"(?:a |an )?(sticks?|lines?|spheres?|balls?|surface|cartoon|ribbons?|"
                r"mesh|dots?) $")
    if match:
        rep_name = {"stick": "sticks", "line": "lines", "sphere": "spheres",
                    "ball": "spheres", "balls": "spheres", "ribbon": "cartoon",
                    "ribbons": "cartoon", "dot": "dots"}.get(match.group(2), match.group(2))
        return ("show %s, chain %s" % (rep_name, match.group(1).upper()),
                "draw chain %s as %s" % (match.group(1).upper(), rep_name))

    # --- a stretch or a list of residues ------------------------------------------
    stretch = has(r"\bresidues?\s+(\d+)(?:\s*(?:to|-|through|until|and|,)\s*|\s+)(\d+)\b"
                  r"((?:\s*(?:and|,)?\s*\d+)*)")
    if stretch and not has(r"\b(distance|measure|how far|how close|between|apart|change|"
                           r"mutat\w*|replace|swap|what if|angle)\b"):
        listed = [stretch.group(1), stretch.group(2)] + re.findall(r"\d+",
                                                                  stretch.group(3))
        is_range = bool(re.search(r"\d\s*(to|-|through|until)\s*\d", stretch.group(0)))
        words_out = "%s to %s" % (listed[0], listed[1]) if is_range and len(listed) == 2 \
            else " ".join(listed)
        selection = "resi %s" % ("%s-%s" % (listed[0], listed[1]) if is_range
                                 and len(listed) == 2 else "+".join(listed))
        colour = has(r"\b(%s)\b" % colours)
        if colour and has(r"\b(colou?r|make|turn|paint)\b"):
            name = colour.group(1).replace("gray", "grey")
            return ("color %s, %s" % (name, selection),
                    "colour residues %s %s" % (words_out, name))
        if has(r"\b(zoom|focus|go)\b"):
            return "zoom %s" % selection, "zoom to residues %s" % words_out
        if has(r"\b(show|see|display|highlight|select|mark|find|where|view)\b"):
            return ("choose residues %s" % words_out,
                    "select residues %s and show PyMOL's own line" % words_out)

    # --- deleting one thing ---------------------------------------------------------
    match = has(r"^ (?:delete|remove|unload|close) (?:the )?(.+?) $")
    if match:
        target = match.group(1).strip()
        if base._known_object(target):
            return ("delete %s" % base._known_object(target),
                    "remove %s from the session (cannot be undone)" % target)
        if re.match(r"^(last )?(measurement|distance|dashed line)s?$", target):
            return "delete measure_*", "remove the measurements"

    # --- a bound molecule by name: "what is ZN", "what is the zinc doing" -----------
    match = has(r"^ what (?:is|does) (?:the |this |that )?([a-z0-9]{2,12})"
                r"( doing| do| for| here| there)? $")
    if match and base._molecules():
        code = _inside_current(match.group(1))
        if code and code.isupper():
            if match.group(2):
                return ("pocket %s %s" % (here, code),
                        "show %s and the residues that hold it" % code)
            return "focus %s" % code, "zoom to %s and say what it is" % code

    # --- yes/no questions about the protein -------------------------------------------
    if has(r"^ is (this|it|that|the protein|this protein|[a-z0-9\-]+) (a|an) "
           r"(enzyme|dimer|monomer|trimer|tetramer|oligomer|kinase|receptor|membrane "
           r"protein|hormone|complex|antibody|homodimer)\b"):
        subject = _subject(t, r"\b(is|enzyme|dimer|monomer|trimer|tetramer|oligomer|kinase|"
                           r"receptor|membrane|hormone|complex|antibody|homodimer)\b")
        return "about %s" % subject, "show the info card (it answers this)"

    # --- measuring --------------------------------------------------------------
    if has(r"\b(distance|measure|how far|how close)\b"):
        ligand = re.findall(r"\b([A-Z0-9]{2,3})\b", text)
        ligand = [w for w in ligand if not w.isdigit() and cmd.count_atoms(
            "resn %s and not polymer" % w)] if base._molecules() else []
        for word in t.split():          # "the zinc" -> ZN
            code = _inside_current(word) if word in _PART_WORDS else None
            if code and code.isupper() and code not in ligand:
                ligand.append(code)
                t = re.sub(r"\b%s\b" % word, " ", t)
        things = [s for s, _ in specs] + ligand
        if len(things) < 2:             # "residues 248 and 273": bare numbers count too
            rest = t
            for _, spec_text in specs:
                rest = rest.replace(spec_text, " ")
            for number in re.findall(r"\b(\d+[a-z]?)\b", rest):
                if number.upper() not in things:
                    things.append(number.upper())
                    t = re.sub(r"\b%s\b" % number, " ", t)
        if len(things) >= 2:
            cleaned = t
            for _, raw in specs:
                cleaned = cleaned.replace(raw, " ")
            for word in ligand:
                cleaned = cleaned.replace(word.lower(), " ")
            subject = _subject(cleaned, r"\b(distance|measure|how far|how close|how long|"
                               r"between|apart|is|are|from|residues?|and|to)\b")
            return ("measure %s" % " ".join(([subject] if subject else []) + things[:2]),
                    "measure the distance between %s and %s" % (things[0], things[1]))
        return "measure", "show how to measure a distance"

    # --- one residue ------------------------------------------------------------
    if len(specs) == 1 and has(r"\b(explain|what is|what's|tell|describe|about|role of|"
                               r"what does|why is|important|near|neighbou?rs?|around|"
                               r"surround\w*|touch\w*|contacts?|interacts?)\b") or (
            len(specs) == 1 and has(r"\b(what|which|why|is)\b") and not has(
                r"\b(where|zoom|focus|show|go|find|see)\b")):
        subject = _subject(t.replace(specs[0][1], " "),
                           r"\b(explain|what is|what's|tell|describe|role|does|do|why|"
                           r"important|residue|chain [a-z0-9]|near|neighbou?rs?|around|"
                           r"surround\w*|touch\w*|contacts?|interacts?|with|show|see|find|"
                           r"what|which|kind|type|amino acid|is|on|surface|buried|"
                           r"exposed|inside|outside)\b")
        return ("explain %s" % " ".join(([subject] if subject else []) + [specs[0][0]]),
                "explain residue %s" % specs[0][0])

    # --- comparing --------------------------------------------------------------
    match = has(r"\b(?:compare|superpose|superimpose|align|overlay)\s+(.+?)\s+"
                r"(?:and|with|to|vs\.?|versus|onto|on|on top of|against)\s+(.+?) $") or \
        has(r"\bput\s+(.+?)\s+(?:onto|on top of)\s+(.+?) $") or \
        has(r"\brmsd\s+(?:between|of|for)\s+(.+?)\s+and\s+(.+?) $") or \
        has(r"^ is\s+(.+?)\s+(?:the same as|identical to|different from|similar to|"
            r"like)\s+(.+?) $") or \
        has(r"\bdifferences?\s+between\s+(.+?)\s+and\s+(.+?) $") or \
        has(r"\bhow (?:different|similar|close|related) (?:is|are)\s+(.+?)\s+"
            r"(?:and|from|to|compared to|compared with)\s+(.+?) $") or \
        has(r"^ (.+?)\s+(?:vs\.?|versus)\s+(.+?) $")
    if match:
        def side(words):
            if words.strip() in ("this", "it", "this one", "that", "the loaded one",
                                 "what is loaded", "current") and here:
                return here
            return _subject(words) or words.strip()
        first, second = side(match.group(1)), side(match.group(2))
        return ("compare %s, %s" % (first, second),
                "superpose %s and %s and show where they differ" % (first, second))

    if has(r"\b(colou?rs?|red|orange|yellow|blue|green|grey|gray|white)\b.*\b(mean|"
           r"means|stand for|thing|part|bit|ball|sphere|stick)s?\b") or has(
            r"\bwhy is (it|this|that|everything)\b.*\b(red|orange|yellow|blue|green|"
            r"grey|gray|white|colou?red)\b"):
        return "glossary colours", "explain what the colours mean"
    if has(r"^ (zoom|zoom in|move|move it|pan|drag|click|select|use the mouse|mouse|"
           r"mouse controls) $"):
        return "start", "show the mouse controls and the first steps"

    # --- facts about the structure: all of them are on the info card ------------
    if has(r"\b(how many (chains?|subunits?|residues?|amino acids?|atoms?|domains?)|"
           r"(which|what) chains?|"
           r"how (big|large|long|heavy|good|reliable|accurate|old)|(size|length|weight|"
           r"mass) of|organism|species|come from|comes from|resolution|quality of|"
           r"when was|what year|who (solved|made|determined|discovered|published)|"
           r"(what|which) (method|technique|experiment)|how was|x-?ray|crystal structure|"
           r"function|job of|role of|purpose of|why is|what disease|diseases?\b.*\b"
           r"(linked|associated|related|cause[sd]?)|linked to|what kind of|what type of|"
           r"what does .{1,40} do|what is .{1,30} for)\b") and not specs and not has(
            r"\b(find|search|are there|is there|load|get|fetch|grab|open|structures? "
            r"(of|for)|(high|good|best|atomic)[- ]resolution)\b"):
        subject = _subject(t, r"\b(how many|how|chains?|subunits?|residues?|amino acids?|"
                           r"atoms?|domains?|big|large|long|heavy|good|reliable|accurate|"
                           r"old|size|length|weight|mass|organism|species|comes?|"
                           r"resolution|quality|when|was|what|year|who|solved|"
                           r"determined|discovered|published|which|method|technique|"
                           r"experiment|x-?ray|crystal|function|job|role|purpose|why|"
                           r"diseases?|linked|associated|related|causes?|caused|kind|"
                           r"type|do|does|for|has|have|an|tell|show|give|info)\b")
        return "about %s" % subject, "show the info card%s (it lists these facts)" % (
            " of %s" % subject if subject else "")
    if has(r"\bhow does\b.*\b(work|function)\b") or has(r"\bhow\b.*\bworks\b"):
        subject = _subject(t, r"\b(how|does|work|works|function)\b")
        return "tour %s" % subject, "start a guided tour"

    # --- small things PyMOL does with one of its own commands -------------------
    if has(r"\b(clear|wipe|clean) (the )?(text|log|output|terminal|messages)\b"):
        return "clear text", "wipe the text in the output area"
    if has(r"\b(delete|clear|remove|wipe|close|empty) (everything|all|it all|the "
           r"(screen|page|window|view|session))\b") or t.strip() == "clear" or has(
            r"\b(start over|start again|start from scratch|new session)\b"):
        return "clear", "empty the 3D window (back restores it)"
    if has(r"\b(hide|remove|without|no|get rid of|delete)\b.*\bwaters?\b"):
        return "hide everything, solvent", "hide the water molecules"
    if has(r"\b(transparent|see[- ]through|translucent)\b") and has(r"\bsurface\b") \
            and len(t.split()) <= 5:
        return ("set transparency, 0.5",
                "make the surface see-through (back to solid: set transparency, 0)")
    if has(r"\b(transparent|see[- ]through|translucent|ghost)\b"):
        return ("set cartoon_transparency, 0.6",
                "make the ribbon see-through (back to solid: set cartoon_transparency, 0)")
    if has(r"\b(spin|rotate|rotating|turn (it )?around|rock)\b") and not specs:
        return "rock 1", "rock the structure back and forth (stop it:  rock 0)"
    if has(r"\bhydrogen bonds?\b") and has(r"\b(ligand|drug|inhibitor|substrate|"
                                             r"cofactor)\b"):
        subject = _subject(t, r"\b(hydrogen|bonds?|ligand|drug|inhibitor|substrate|"
                           r"cofactor|show|draw|see|makes?|made|forms?|does)\b")
        return "pocket %s" % subject, "show the ligand, its pocket and its hydrogen bonds"
    if has(r"\bhydrogen bonds?\b") and not has(r"\b(ligand|drug|inhibitor|substrate|"
                                                 r"cofactor)\b") and not specs:
        return ("distance hydrogen_bonds, polymer, polymer, mode=2",
                "draw the hydrogen bonds inside the protein as dashed lines "
                "(remove them: delete hydrogen_bonds)")
    colour = has(r"\b(red|blue|green|yellow|orange|purple|pink|white|grey|gray|cyan|"
                 r"magenta|black|salmon|violet|brown)\b")
    if colour and has(r"\b(colou?r|make|turn|highlight)\b"):
        targets = [(r"\bheli(x|ces)\b", "ss h", "the helices"),
                   (r"\b(strands?|sheets?)\b", "ss s", "the strands"),
                   (r"\b(loops?|coils?)\b", "ss l+''", "the loops"),
                   (r"\b(ligands?|drug|inhibitor|cofactor|substrate)\b", "organic",
                    "the bound molecules"),
                   (r"\b(waters?)\b", "solvent", "the water"),
                   (r"\b(everything|it|all|protein|structure|whole)\b", "all",
                    "everything")]
        chain = has(r"\bchain\s+([a-z0-9])\b")
        amino = has(r"\b(%s)\b" % _AA_PATTERN)
        name = {"gray": "grey", "purple": "purple", "brown": "brown"}.get(
            colour.group(1), colour.group(1))
        for word in t.split():
            code = _inside_current(word) if word in _PART_WORDS else None
            if code and code.isupper():
                return "color %s, resn %s" % (name, code), "colour %s %s" % (code, name)
        if chain:
            return ("color %s, chain %s" % (name, chain.group(1).upper()),
                    "colour chain %s %s" % (chain.group(1).upper(), name))
        if amino:
            code = tools.ONE_TO_THREE[_AA_WORDS[amino.group(1)]]
            return ("color %s, resn %s" % (name, code),
                    "colour every %s %s" % (tools.AA[code][0], name))
        for pattern, selection, description in targets:
            if has(pattern):
                return ("color %s, %s" % (name, selection),
                        "colour %s %s" % (description, name))

    # --- pictures and files -----------------------------------------------------
    if has(r"\b(save|take|make|export|create|give)\b.*\b(picture|image|screenshot|png|"
           r"photo|snapshot)\b") or has(r"\b(screenshot|snapshot)\b"):
        name = (base._last[0] or "picture") if base._molecules() else "picture"
        named = re.search(r"\b(?:as|called|named|to)\s+([A-Za-z0-9_\-]+)(?:\.png)?\s*$",
                          text.strip())
        if named and named.group(1).lower() not in ("png", "image", "picture", "file"):
            return ("snap %s" % named.group(1),
                    "save a high-resolution picture as %s.png" % named.group(1))
        return "snap %s_picture" % name, "save a high-resolution picture of the current view"
    if has(r"\b(report|hand in|hand-in|summary page|write up|write-up|one page)\b"):
        subject = _subject(t, r"\b(report|make|write|create|hand in|hand-in|summary|page|"
                           r"up|give|generate|write-up|one)\b")
        return "report %s" % subject, "write a one-page report with pictures"
    if has(r"\bsave\b.*\b(work|session|everything|progress)\b") or has(
            r"^ save (this|it|that|all)? ?$") or has(r"\b(keep|store) (my|this|the) "
                                                   r"(work|session)\b"):
        named = re.search(r"\b(?:as|called|named)\s+([A-Za-z0-9_\-]+)\s*$", text.strip())
        if named:
            return ("session save %s" % named.group(1),
                    "save the whole session under the name %s" % named.group(1))
        return "session save my_work", "save the whole session under the name my_work"
    if has(r"\b(open|load|restore)\b.*\b(work|session|progress)\b"):
        return "session", "list your saved sessions"
    if has(r"\b(teach|pymol commands|real commands|raw commands|commands behind|"
           r"what commands)\b"):
        return "prospector teach on", "show the PyMOL commands behind every action"
    if has(r"\b(sequence|amino acid order|primary structure)\b") and not has(r"\bmotif\b"):
        subject = _subject(t, r"\b(sequence|show|print|give|display|amino acid order|"
                           r"primary|get|letters)\b")
        return "sequence %s" % subject, "print the sequence with helices and strands marked"
    styled = any(re.search(r"\b%s\b" % re.escape(w), t) for w in _STYLE_WORDS)
    if has(r"\b(similar|related|relatives|homolog|homologs|homologues|looks? like)\b") \
            and not has(r"\bwhat does\b.*\blooks? like\b") \
            and not styled and not has(r"\b(publication|figure)\b"):
        subject = _subject(t, r"\b(similar|related|relatives|homologs?|homologues|find|"
                           r"show|list|look like|others|other|search|to)\b")
        return "similar %s" % subject, "find related structures by sequence"
    if has(r"\b(tour|walk me through|guide me|guided|show me around|walk through)\b"):
        subject = _subject(t, r"\b(tour|walk|through|guide|guided|show|around|give|take|"
                           r"start)\b")
        return "tour %s" % subject, "start a guided tour"
    if has(r"\b(quiz|test me|test my|challenge)\b"):
        return "quiz", "start a quiz with a mystery structure"
    if has(r"\b(surprise|random)\b"):
        return "surprise", "load a random classic structure"
    if has(r"\b(alphafold|predicted|prediction|predict)\b"):
        subject = _subject(t, r"\b(alphafold|predicted|prediction|predict|show|load|get|"
                           r"fetch|open|model|shape|fold)\b")
        return "alphafold %s" % subject, "load the AlphaFold predicted model"

    # --- how it is drawn --------------------------------------------------------
    for pattern, style, description in [
            (r"\b(every|all|each|the|individual|single) atoms?\b|\batom by atom\b",
             "ballstick", "every atom as a ball, every bond as a stick"),
            (r"\b(the|all|every|chemical|covalent) bonds\b", "sticks",
             "every bond as a stick"),
            (r"\b(backbone|main chain|trace|skeleton)\b", "ribbon",
             "only the path of the backbone, as a thin ribbon")]:
        if has(pattern) and has(r"\b(show|view|see|draw|make|only|just)\b"):
            subject = _subject(t, pattern, r"\b(show|view|see|draw|make|as|with)\b")
            return ("viewas %s" % " ".join(([subject] if subject else []) + [style]),
                    "draw %s" % description)
    for phrase in sorted(_STYLE_WORDS, key=len, reverse=True):
        if re.search(r"\b%s\b" % re.escape(phrase), t) and has(
                r"\b(show|view|display|draw|make|render|as|in|with|see|switch|change|"
                r"turn|representation|style|mode|look|looks|like|to)\b") and not has(
                r"\b(colou?r|highlight|buried|exposed|residues?|amino acids?|which|"
                r"where)\b"):
            style = _STYLE_WORDS[phrase]
            subject = _subject(t, r"\b%s\b" % re.escape(phrase),
                               r"\b(show|view|display|draw|make|render|as|see|switch|"
                               r"change|turn|representation|style|mode|atoms?|bonds?|looks?|like|"
                               r"back|to)\b")
            return ("viewas %s" % " ".join(([subject] if subject else []) + [style]),
                    "draw %s as %s" % (subject or "the structure", style))
    looks = [(r"\b(by chain|per chain|each chain|chains? (a )?different|colou?r chains)\b",
              "chains", "colour each chain differently"),
             (r"\b(rainbow|n to c|n-terminus to c|start to end)\b", "rainbow",
              "colour each chain from blue (start) to red (end)"),
             (r"\b(hydrophobicity|by hydrophob\w*|hydrophobic surface|polarity)", "hydrophobicity",
              "colour the surface by hydrophobicity"),
             (r"\b(b-factor|b factor|bfactor|flexib\w*|mobil\w*|putty)", "putty",
              "show flexibility as a tube (B-factor putty)"),
             (r"\b(secondary structure|by (shape|fold|helix|ss)|helices and "
              r"(strands|sheets))\b", "secondary",
              "colour helices red, strands yellow and loops green"),
             (r"\b(by (element|atom|atom type)|element colou?rs?|cpk)\b", "elements",
              "show every atom, coloured by element"),
             (r"\b((residue|amino acid) (type|kind|class|property|properties)|by "
              r"(property|properties|chemistry|type)|chemical (property|properties))\b",
              "types", "colour by the kind of amino acid"),
             (r"\b(publication|figure|paper|pretty|nice|clean|presentation|poster)\b",
              "publication", "apply a clean publication style"),
             (r"\b(buried ligand|ligand .*surface|surface .*ligand|how buried)\b", "ligand",
              "show the surface with the ligand inside")]
    for pattern, preset, description in looks:
        if has(pattern) and (preset in ("secondary", "elements", "types", "publication")
                             or has(r"\b(colou?rs?|look|style|make|show|view|picture|"
                                    r"see|display|give)\b")):
            subject = _subject(t, pattern, r"\b(colou?r|look|style|make|show|view|picture|"
                               r"display|by|surface|as|give|structure|colou?rs|elements?)\b")
            return ("look %s" % " ".join([preset] + ([subject] if subject else [])),
                    description)

    # --- places -----------------------------------------------------------------
    if has(r"\b(pocket|what binds|where does .* bind|binds? to|ligand interactions?|"
           r"interactions? (with|of) (the )?(ligand|drug|inhibitor|substrate|cofactor)|"
           r"drug bind|inhibitor bind|binding pocket|ligand contacts)\b"):
        subject = _subject(t, r"\b(pocket|what binds|binds?|binding|ligand|interactions?|"
                           r"drug|inhibitor|substrate|cofactor|show|where|does|contacts|"
                           r"the|see|display)\b")
        return "pocket %s" % subject, "show the binding pocket of the ligand"
    near = has(r"\b(touch\w*|holds?|holding|surround\w*|coordinat\w*|around|near|next to|"
               r"neighbou?rs? of|contacts?|interacts? with|binds?)\b\s+(?:the |a |an )?"
               r"([a-z0-9]+)\b")
    if near and (near.group(2) in _PART_WORDS or _inside_current(near.group(2))):
        code = _inside_current(near.group(2))
        obj = base._last[0] if base._molecules() else ""
        if code and len(code) <= 3 and code.isupper() and obj:
            return ("pocket %s %s" % (obj, code),
                    "show what surrounds %s in %s" % (code, obj))
        subject = _subject(t, near.re.pattern, r"\b(what|which|show|residues?|amino "
                           r"acids?|the|does|do|in|how)\b")
        return "pocket %s" % subject, "show the binding pocket of the ligand"
    places = [(r"\b(zoom out|whole (thing|structure|protein)|everything|reset (the )?view|"
               r"full view|show all of it|^ reset $|recent(er|re)|cent(er|re) (the )?"
               r"(view|it|protein|structure|molecule)|fit (it )?(to|on) (the )?screen)\b"
               r"|^ reset $", "all", "the whole structure"),
              (r"\b(active site)\b",
               "active site", "the active site"),
              (r"\b(binding sites?|metal sites?)\b", "binding sites", "the binding sites"),
              (r"\b(disease variants?|variants?)\b", "variants",
               "the disease-variant positions"),
              (r"\b(drug|inhibitor|ligand|cofactor|substrate)\b", "ligand", "the ligand"),
              (r"\b(domains?)\b", "domain", "the first domain")]
    for pattern, place, description in places:
        if has(pattern):
            subject = _subject(t, pattern, r"\b(zoom|focus|go|take|show|see|look|where|"
                               r"find|view|display|reset|bring|in|on|to|closer|what is|"
                               r"what's|what are|which|whats|highlight|point|tell|does|"
                               r"do|have|has|residues?|located|location|parts?|are|is|amino acids?|"
                               r"positions?|known|again|me|as|sticks|goes|go|sits?)\b")
            if subject:
                return ("focus %s of %s" % (place, subject),
                        "highlight and zoom to %s of %s" % (description, subject))
            return "focus %s" % place, "highlight and zoom to %s" % description
    match = has(r"\bchain\s+([a-z0-9])\b")
    if match and not specs and has(r"\b(zoom|focus|go|show|see|look|view|display)\b"):
        subject = _subject(t, r"\bchain\s+[a-z0-9]\b", r"\b(zoom|focus|go|show|see|look|"
                           r"view|display|only|just)\b")
        return ("focus chain %s%s" % (match.group(1).upper(), " " + subject if subject
                                      else ""), "zoom to chain %s" % match.group(1).upper())
    if len(specs) == 1 and has(r"\b(zoom|focus|go|show|see|look|where|find|view|display)\b"):
        subject = _subject(t.replace(specs[0][1], " "), r"\b(zoom|focus|go|show|see|look|"
                           r"where|find|view|display|residue|chain [a-z0-9])\b")
        return ("focus %s%s" % (specs[0][0], " " + subject if subject else ""),
                "zoom to residue %s" % specs[0][0])

    # --- spotting things --------------------------------------------------------
    match = re.search(r"\bmotif\s+([a-z\[\]\^]{2,30})\b", raw)
    if match:
        raw = re.search(r"motif\s+(\S+)", text, re.I).group(1)
        subject = _subject(t, r"\bmotif\s+\S+", r"\b(find|show|highlight|spot|where|"
                           r"locate|sequence|pattern|search)\b")
        return ("scout %s" % " ".join(([subject] if subject else []) + ["motif", raw]),
                "find the sequence motif %s" % raw)
    for phrase in _SCOUT_PHRASES:
        if re.search(r"\b%s\b" % re.escape(phrase), t):
            subject = _subject(t, r"\b%s\b" % re.escape(phrase),
                               r"\b(highlight|find|show|spot|locate|mark|where|display|"
                               r"see|colou?r|list|residues?|bonds?|count|how many|"
                               r"parts?|regions?|amino acids?|are|is|which|what)\b")
            return ("scout %s" % " ".join(([subject] if subject else []) + [phrase]),
                    "highlight: %s" % phrase)

    if has(r"\b(features?|annotations?|important (residues|parts|sites)|functional sites|"
           r"key residues|sites and domains)\b"):
        subject = _subject(t, r"\b(features?|annotations?|important|residues|parts|sites|"
                           r"functional|key|show|domains|highlight|colou?r|map|display)\b")
        return "features %s" % subject, "colour the annotated sites and domains"
    if has(r"\b(key places|interesting places|where can i zoom|places)\b"):
        return "focus", "list the key places you can zoom to"

    # --- "zoom to the zinc", "where is the heme in haemoglobin" -------------------
    if has(r"\b(zoom|focus|where)\b"):
        rest = re.sub(r"\b(zoom|focus|where|is|are|to|in|on|closer|located|sits?|find|"
                      r"show|go)\b", " ", t)
        words = [w for w in re.findall(r"[a-z0-9][a-z0-9\-\.']*", rest) if w not in _STOP]
        subject = " ".join(words)
        if subject and base._known_object(subject):
            return "focus all %s" % subject, "show the whole of %s" % subject
        if subject and _inside_current(subject):
            inside = _inside_current(subject)
            return "focus %s" % inside, "zoom to %s in the loaded structure" % inside
        for size in range(len(words) - 1, 0, -1):       # "<thing> <protein>" either way
            for start in (0, len(words) - size):
                name = " ".join(words[start:start + size])
                thing = " ".join(words[:start] + words[start + size:])
                if base._known_object(name) or (set(thing.split()) <= _PART_WORDS):
                    return ("focus %s of %s" % (thing, name),
                            "zoom to %s in %s" % (thing, name))

    # --- about, find, load ------------------------------------------------------
    if has(r"\b(what is|what's|tell me about|tell about|about|describe|information|info "
           r"on|details|who is|what does .* do|function of|explain)\b"):
        subject = _subject(t, r"\b(what is|what's|tell|about|describe|information|info|"
                           r"details|who is|what does|function|explain|do|does|more|"
                           r"give|show|get)\b")
        return "about %s" % subject, "show the info card%s" % (
            " of %s" % subject if subject else "")
    filters, t = _filters_in(t)
    if has(r"\b(find|search|look for|look up|which|list|are there|is there|options for|"
           r"candidates|best structure)\b") or (has(r"\bstructures? (of|for)\b") and not
                                                 has(r"\b(load|fetch|get|open|grab|show|"
                                                     r"download|see|need|want)\b")):
        subject = _subject(t, r"\b(find|search|look for|look up|which|list|are there|"
                           r"is there|options|candidates|best|structures?|of|for|"
                           r"available|exist|pdb entries|entries|version)\b")
        if subject and _name_like(subject, text, True):
            return ("find %s" % " ".join([subject] + filters),
                    "search for structures of %s" % subject)
    if has(r"\b(load|fetch|get|open|grab|download|bring up|show|display|see|view|give|"
           r"need|want|draw)\b"):
        subject = _subject(t, r"\b(load|fetch|get|open|grab|download|bring up|show|"
                           r"display|see|view|give|want|need|version|of|for|draw)\b")
        if subject:
            if base._known_object(subject):
                return "focus all %s" % subject, "show the whole of %s" % subject
            inside = _inside_current(subject)
            if inside:
                return "focus %s" % inside, "zoom to %s in the loaded structure" % inside
            wants_new = has(r"\b(load|fetch|get|open|grab|download|bring up|need)\b")
            if not _name_like(subject, text, wants_new):
                return None, None
            a_part = bool(set(subject.split()) & _PART_WORDS)   # "the heme", "the zinc"
            if wants_new or not base._molecules() or not a_part:
                return ("grab %s" % " ".join([subject] + filters),
                        "find and load the best structure of %s" % subject)
            # "show me X" with a structure on screen and X not part of it: loading
            # a different protein would be a guess, so say so instead
            return None, ("'%s' is not something I can find in the loaded structure. "
                          "To load a protein called that:  grab %s" % (subject, subject))
    # just a name: "p53"
    words = t.split()
    if 1 <= len(words) <= 3 and not any(w in _VERBS or w in _STOP for w in words):
        name = " ".join(words)
        if base._known_object(name):
            return "about %s" % name, "show the info card of %s" % name
        if not _name_like(name, text, False):
            return None, None
        return "find %s" % name, "search for structures of %s (nothing is loaded)" % name
    return None, None


from ..core.ai import PYMOL_COMMANDS as _AI_OWN
# PyMOL's own commands the RULES may run.  `delete` and `move` are here only:
# the rules use them for exactly what the user typed; the AI may not.
_PYMOL_OWN = tuple(_AI_OWN) + ("delete", "move")


def _filters_in(t):
    """Search filters said in words -> (['method=em', ...], sentence without them)."""
    found = []
    for pattern, text in [(r"\b(cryo[- ]?em|electron microscopy)\b", "method=em"),
                          (r"\b(x-?ray|crystal|crystallography)\b", "method=xray"),
                          (r"\bnmr\b", "method=nmr"),
                          (r"\b(high|good|best|atomic)[- ]resolution\b", "res<2.0")]:
        if re.search(pattern, t):
            found.append(text)
            t = re.sub(pattern, " ", t)
    match = re.search(r"\b(?:with|bound to|containing|in complex with)\s+([a-z0-9]{2,3})"
                      r"(?:\s+bound)?\b", t)
    if match and match.group(1) not in _STOP and match.group(1) not in ("dna", "rna"):
        found.append("ligand=%s" % match.group(1).upper())
        t = t.replace(match.group(0), " ")
    return found, t


_BUILD_KINDS = [
    (r"\b(double[- ]stranded dna|dna double helix|double helix|dna)\b", "dna"),
    (r"\brna\b", "rna"),
    (r"\b(beta[- ]sheets?|sheets?)\b", "sheet"),
    (r"\b(beta[- ]strands?|strands?)\b", "strand"),
    (r"\b(alpha[- ]helix|alpha[- ]helices|helix|helices)\b", "helix"),
    (r"\b(peptides?|polypeptide)\b", "peptide"),
    (r"\b(lipid bilayers?|bilayers?|membranes?|lipid layers?)\b", "membrane"),
    (r"\b(micelles?)\b", "micelle"),
]
_ASSEMBLY = (r"\b(biological (assembly|unit)|assembly|(?:whole|complete|full|entire) "
             r"(?:(?P<mid>\S+) )?(?:complex|cage|capsid|virus|shell|oligomer)|all (the )?"
             r"(chains|subunits|copies))\b")


def _build_request(text, raw):
    """'build a dna double helix of ATGCGT', 'make a membrane with 30% cholesterol',
    'show the whole ferritin cage' -> a build command, or None."""
    verb = re.search(r"\b(build|make|create|generate|construct|model|assemble|draw)\b", raw)
    assembly = re.search(_ASSEMBLY, raw)
    if assembly and (verb or re.search(r"\b(show|see|load|get|view|display)\b", raw)):
        rest = raw.replace(assembly.group(0), " %s " % (assembly.group("mid") or ""))
        subject = _subject(" %s " % rest, r"\b(build|make|create|generate|construct|model|"
                           r"assemble|draw|show|see|load|get|view|display|complete|whole|"
                           r"full|biological|cage|capsid|complex|shell|virus|unit)\b")
        here = base._last[0] if base._last[0] in base._molecules() else ""
        target = subject or here
        return ("build assembly %s" % target if target else "build assembly",
                "build the complete biological assembly%s" % (" of " + target if target
                                                              else ""))
    if not verb:
        return None
    kind = None
    for pattern, name in _BUILD_KINDS:
        if re.search(pattern, raw):
            kind = name
            break
    if kind is None:
        return None
    letters = {"dna": "ACGT", "rna": "ACGU"}.get(kind, "ACDEFGHIKLMNPQRSTVWY")
    sequences = []
    if kind not in ("membrane", "micelle"):
        for word in re.findall(r"[A-Za-z]+", text):
            if len(word) < 3 or not set(word.upper()) <= set(letters):
                continue
            # a sequence is typed in capitals; a long lower-case one only for DNA/RNA,
            # and never an ordinary word ("create" uses only amino-acid letters too)
            if word.isupper() or (kind in ("dna", "rna") and len(word) >= 6
                                  and word.lower() not in _known_words()):
                if word.upper() not in ("DNA", "RNA"):
                    sequences.append(word.upper())
    if kind == "dna" and re.search(r"\ba[- ]form\b|\ba-dna\b", raw):
        sequences.append("form=A")
    if kind in ("dna", "rna", "helix", "strand", "peptide") and sequences:
        sequences = [" ".join(s for s in sequences if not s.startswith("form="))] + \
            [s for s in sequences if s.startswith("form=")]
    extra = []
    if kind == "membrane":
        for lipid in ("popc", "dopc", "dope"):
            if re.search(r"\b%s\b" % lipid, raw):
                extra.append(lipid)
        size = re.search(r"\b(\d{1,2})\s*(?:nm|nanometres?|nanometers?|x\s*\d{1,2})", raw)
        if size:
            extra.append(size.group(1))
        chol = re.search(r"\b(\d{1,2})\s*(?:%|percent|per cent)\s*(?:of\s*)?cholesterol\b", raw)
        if chol:
            extra.append("cholesterol=%s" % chol.group(1))
        elif re.search(r"\bcholesterol\b", raw):
            extra.append("cholesterol=30")
    if kind == "micelle":
        for detergent in ("sds", "ldao"):
            if re.search(r"\b%s\b" % detergent, raw):
                extra.append(detergent)
        count = re.search(r"\b(\d{2,3})\b", raw)
        if count:
            extra.append(count.group(1))
    command = " ".join(["build", kind] + sequences + extra)
    what = {"dna": "a DNA double helix", "rna": "an RNA strand", "helix": "an alpha-helix",
            "strand": "a beta-strand", "sheet": "a beta-sheet", "peptide": "a stretched peptide",
            "membrane": "a lipid membrane", "micelle": "a detergent micelle"}[kind]
    return command, "build %s%s" % (what, " of " + " ".join(sequences) if sequences else "")


def _part_words():
    from . import places
    words = set(places.PLACE_WORDS) - set(["the", "of", "in", "to", "structure", "all",
                                            "whole", "everything", "back", "out"])
    for names in places.LIGAND_WORDS.values():
        words.update(names)
    words.update("metal metals ion ions water waters pocket cofactor inhibitor drug "
                 "substrate helix helices strand strands loop loops terminus".split())
    return words


_PART_WORDS = _part_words()


def _inside_current(subject):
    """If ``subject`` names something in the structure being worked on (a bound
    molecule such as 'heme' or 'HEM', or one of its key places), return the
    words to give to `focus`."""
    from . import places
    names = base._molecules()
    obj = base._last[0] if base._last[0] in names else (names[-1] if names else None)
    if obj is None:
        return None
    word = subject.lower()
    codes = set()
    cmd.iterate("%s and not polymer and not solvent" % obj, "codes.add(resn)",
                space={"codes": codes})
    for code in codes:
        if word == code.lower() or word in places.LIGAND_WORDS.get(code, []):
            return code
    listed = places._places.get(obj)
    if listed and places._match(listed, word) is not None:
        return subject
    return None


_EXAMPLES = [
    ("load haemoglobin", "grab haemoglobin"),
    ("what is p53", "about p53"),
    ("show the active site", "focus active site"),
    ("where does the drug bind", "pocket"),
    ("show it as sticks", "viewas sticks"),
    ("colour by chain", "look chains"),
    ("highlight the histidines", "scout histidine"),
    ("compare 1HHO and 4HHB", "compare 1HHO, 4HHB"),
    ("change arginine 248 to glutamine", "whatif R248Q"),
    ("distance between A/87 and HEM", "measure A/87 HEM"),
    ("explain residue A/248", "explain A/248"),
    ("give me a tour", "tour"),
    ("save a picture", "snap my_picture"),
    ("what does RMSD mean", "glossary rmsd"),
    ("undo", "back"),
]


def _split_steps(text):
    """'load p53 and show it as sticks' -> two requests.  'compare a and b' and
    'ball and stick' stay whole."""
    text = re.sub(r"\bball(s)?\s+and\s+stick(s)?\b", "ballstick", text, flags=re.I)
    parts = re.split(r"\s*(?:;|(?<=[a-z0-9]{3})\.\s+|\bthen\b|\band then\b|\bafter that\b|,\s*(?=(?:%s)\b))"
                     r"\s*" % "|".join(_VERBS), text)
    steps = []
    for part in parts:
        if not part.strip():
            continue
        if re.match(r"^\s*(compare|superpose|superimpose|align|overlay|differences?|"
                    r"distance|measure|how far)\b", part, re.I):
            steps.append(part.strip())
            continue
        pieces = re.split(r"\s+and\s+(?=(?:%s)\b)" % "|".join(_VERBS), part, flags=re.I)
        steps.extend(p.strip() for p in pieces if p.strip())
    return steps


def ask(*args, **kwargs):
    """
DESCRIPTION

    Say what you want in your own words.  Prospector works out which command
    does it, shows you that command, and runs it.  It uses fixed rules, not an
    AI: when it is not sure it offers choices instead of guessing.

    Start with "how do I ..." to only see the command, without running it.

USAGE

    ask <what you want>

EXAMPLES

    ask load haemoglobin and show it as sticks
    ask where does the drug bind in p53
    ask change arginine 248 to glutamine in p53
    ask what does RMSD mean
    ask how do I colour by chain
    """
    text = " ".join(base._raw(args, kwargs, ", ").split())
    label = "ask %s" % text

    def action():
        base._header(label)
        if not text:
            _say("")
            _say("Say what you want after  ask . For example:")
            for phrase, command in _EXAMPLES:
                _say("   ask %-36s (runs:  %s)" % (phrase, command))
            _say(report.LINE)
            return
        only_show = bool(re.match(r"^\s*(how (do|can|could|would|should) (i|you|we)|"
                                  r"how to|what do i type to|what is the command (to|for)|"
                                  r"which command)\b", text, re.I))
        body = re.sub(r"^\s*(how (do|can|could|would|should) (i|you|we)|how to|"
                      r"what do i type to|what is the command (to|for)|which command)\s*",
                      "", text, flags=re.I)
        steps = _split_steps(body) or [body]
        plans = []
        for step in steps:
            if plans and re.match(r"^\s*(it|that|this|them)\b", step, re.I) is None:
                pass
            command, description = _interpret(step)
            plans.append((step, command, description))
        _say("")
        if any(command is None for _, command, _ in plans):
            unclear, _, reason = [plan for plan in plans if plan[1] is None][0]
            from . import aiask
            from ..core import ai
            here = base._last[0] if base._last[0] in base._molecules() else ""
            learned = ai.learned_get(body, here)
            if learned:              # the AI answered this very sentence before
                _present_and_run([(body, command, description)
                                  for command, description in learned],
                                 only_show, via_ai="learned")
                return
            if aiask.available():
                _say("The built-in rules did not understand '%s'." % unclear)
                _say("Asking the AI (sends your sentence and the names of what is "
                     "loaded) ...")
                reference_text, context_text = aiask.reference(), aiask.context()
                with worker.next_in_line():
                    worker.submit(
                        lambda: aiask.call(body, reference_text, context_text),
                        lambda result, error: _after_ai(result, error, unclear, reason,
                                                        only_show, sentence=body))
                return
            _refuse(unclear, reason, tip=True)
            return
        _present_and_run(plans, only_show)

    worker.submit(lambda: None, lambda result, error: action() if not error
                  else base._unexpected(error))


def askai(*args, **kwargs):
    """
DESCRIPTION

    Like  ask , but the sentence goes straight to the AI: the built-in rules
    are skipped.  Use it for requests in another language, follow-ups that
    depend on what you did before, or when  ask  picked the wrong command.

    Needs a stored key (prospector ai key).  The AI only proposes command
    lines; Prospector checks them, shows them, and then runs them.
    Start with "how do I ..." to see the commands without running them.

USAGE

    askai <what you want>

EXAMPLES

    askai zooma sul sito attivo e evidenzia le cariche negative
    askai show the residues that hold the zinc and colour the rest grey
    askai the same but for chain B
    askai how do I compare normal p53 with a cancer mutant
    """
    text = " ".join(base._raw(args, kwargs, ", ").split())
    label = "askai %s" % text

    def action():
        from ..core import ai
        from . import aiask
        base._header(label)
        _say("")
        if not text:
            _say("Say what you want after  askai . It goes straight to the AI, e.g.")
            _say("   askai show the residues that hold the zinc and colour the rest grey")
            _say("   ( ask  tries the free built-in rules first.)")
            _say(report.LINE)
            return
        if not aiask.available():
            if not ai.api_key():
                _say("askai  needs an AI key, and none is stored yet.")
                _say("SET UP       prospector ai key   (then:  prospector ai test)")
            else:
                _say("The AI is switched off. Switch it on with:  prospector ai on")
            _say("MEANWHILE    ask %s   uses the built-in rules" % text)
            _say(report.LINE)
            return
        only_show = bool(re.match(r"^\s*(how (do|can|could|would|should) (i|you|we)|"
                                  r"how to|what do i type to|what is the command (to|for)|"
                                  r"which command)\b", text, re.I))
        _say("Asking the AI (sends your sentence and the names of what is loaded) ...")
        reference_text, context_text = aiask.reference(), aiask.context()
        with worker.next_in_line():
            worker.submit(lambda: aiask.call(text, reference_text, context_text),
                          lambda result, error: _after_ai(result, error, text, None,
                                                          only_show, direct=True))

    worker.submit(lambda: None, lambda result, error: action() if not error
                  else base._unexpected(error))


def _refuse(unclear, reason, tip=False, did_nothing=True):
    """The honest "I did not understand" answer, with things that do work."""
    if did_nothing:
        _say("I'm not sure what you mean by '%s', so I did nothing." % unclear)
    if reason and reason.startswith("'"):
        _say(reason)
    words = set(re.findall(r"[a-z]+", unclear.lower())) - _STOP
    scored = sorted(_EXAMPLES, key=lambda item: -len(
        words & set(re.findall(r"[a-z]+", (item[0] + " " + item[1]).lower()))))
    _say("")
    _say("THINGS I UNDERSTAND   (say it like this, or type the command on the "
         "right)")
    for phrase, command in scored[:6]:
        _say("   ask %-36s =  %s" % (phrase, command))
    _say("   More:  ask  (all examples)  |  info  (all commands)  |  start")
    if tip:
        _say("   Optional:  prospector ai  - let an AI read sentences the rules miss")
    _say(report.LINE)


def _after_ai(result, error, unclear, reason, only_show, direct=False, sentence=None):
    """GUI-thread half of the AI fallback: check what came back, then run it."""
    from ..core import ai
    if error:
        base._unexpected(error)
        return
    kind, value = result
    _say("")
    if kind == "error":
        _say("AI NOT USED  %s" % value)
        if direct:
            _say("MEANWHILE    ask %s   uses the built-in rules" % unclear[:60])
            _say(report.LINE)
            return
        _refuse(unclear, reason)
        return
    good, bad = ai.validate(value)
    if bad:
        _say("AI NOT USED  it proposed something I do not allow, so I ran nothing:")
        for command, why in bad[:3]:
            _say("             %s   (%s)" % (command[:60], why))
        if direct:
            _say(report.LINE)
            return
        _refuse(unclear, reason, did_nothing=False)
        return
    if not good:
        _say("THE AI SAYS  %s" % (value.cannot or "it cannot turn that into commands."))
        _say("             Nothing was run.   More:  ask  (examples)  |  info  (commands)")
        _say(report.LINE)
        return
    remembered = False
    if not only_show:
        here = base._last[0] if base._last[0] in base._molecules() else ""
        also = [a for a, o in base._aliases.items() if o == here]
        general = ai.generalise(sentence or unclear, good, here, also)
        loaded = list(base._molecules()) + list(base._aliases)
        if ai.learnable(sentence or unclear, general, loaded):
            ai.learned_put(sentence or unclear, general)
            remembered = True
    _present_and_run([(unclear, command, description) for command, description in good],
                     only_show, via_ai=True, remembered=remembered)


def _present_and_run(plans, only_show, via_ai=False, remembered=False):
    tag = " (learned)" if via_ai == "learned" else (" (AI)" if via_ai else "")
    for index, (step, command, description) in enumerate(plans, 1):
        command = " ".join(command.split())
        number = "%d. " % index if len(plans) > 1 else ""
        if only_show:
            _say("%sTO %s" % (number, description.upper()))
            _say("   TYPE        %s" % command)
            first = command.split()[0]
            function = cmd.keyword.get(first, [None])[0]
            doc = (getattr(function, "__doc__", "") or "")
            lines = [l.strip() for l in doc.split("DESCRIPTION")[-1].split("USAGE")[0]
                     .strip().splitlines() if l.strip()]
            if lines and first not in _PYMOL_OWN:
                _say("   WHAT IT DOES  %s" % " ".join(lines)[:230])
        else:
            _say("%sI UNDERSTOOD%s  %s" % (number, tag, description))
            _say("   RUNNING     %s" % command)
    if only_show:
        _say("")
        _say("Nothing was run. Type the line%s above yourself, or drop 'how do I' to "
             "let me run %s." % ("s" if len(plans) > 1 else "",
                                 "them" if len(plans) > 1 else "it"))
        _say(report.LINE)
        return
    _say("   (next time you can type %s directly)"
         % ("those commands" if len(plans) > 1 else "that command"))
    if via_ai == "learned":
        _say("   The AI answered this sentence before, so no request was sent. Wrong? "
             " askai  asks afresh;  prospector ai unlearn  forgets everything learned.")
    elif via_ai:
        _say("   Chosen by the AI, checked by Prospector.  back  undoes a change to "
             "the picture.")
        if remembered:
            _say("   Remembered: next time  ask  runs this sentence at once, without "
                 "the AI.")
    _say(report.LINE)
    # run them now, ahead of anything typed after this `ask`
    with worker.next_in_line():
        for _, command, _ in plans:
            parts = " ".join(command.split()).split(" ", 1)
            if parts[0] in _PYMOL_OWN:          # one of PyMOL's own commands
                cmd.do(command, echo=0)
                continue
            function = cmd.keyword[parts[0]][0]
            if len(parts) > 1:
                function(parts[1])
            else:
                function()




def register():
    table = {"ask": ask, "askai": askai, "start": start, "status": status,
             "glossary": glossary}
    for name, function in table.items():
        cmd.extend(name, function)
        cmd.extend("p_" + name, function)
