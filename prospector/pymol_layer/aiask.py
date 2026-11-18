"""The PyMOL side of the optional AI fallback for `ask` (see core/ai.py).

`reference()` is the command manual the model reads, `context()` describes
what is on screen, `settings()` is `prospector ai ...`.
"""
from __future__ import annotations

import re

from pymol import cmd

from ..core import ai, report
from . import commands as base
from . import worker

_say = base._say
_reference = [None]


def available():
    return ai.enabled()


# ----------------------------------------------------------------- reference
def _section(doc, name):
    """One section (USAGE, EXAMPLES ...) of a PyMOL-style docstring."""
    match = re.search(r"^%s\s*\n(.*?)(?=^\S|\Z)" % name, doc or "", re.M | re.S)
    if not match:
        return []
    return [line.strip() for line in match.group(1).splitlines() if line.strip()]


def reference():
    """The command manual sent to the model.  Built once, and always in the
    same order, so the service can cache it between questions."""
    if _reference[0] is not None:
        return _reference[0]
    from . import tier3, tools
    lines = []
    for name in ai.PROSPECTOR_COMMANDS:
        entry = cmd.keyword.get(name)
        doc = getattr(entry[0], "__doc__", "") if entry else ""
        what = " ".join(_section(doc, "DESCRIPTION"))[:300]
        lines.append("%s - %s" % (name, what))
        for label in ("USAGE", "EXAMPLES"):
            for line in _section(doc, label)[:14]:
                lines.append("    %s" % line)
    lines.append("")
    lines.append("WORDS THE COMMANDS ACCEPT")
    lines.append("viewas styles: %s" % ", ".join(sorted(tier3.STYLE_NAMES)))
    lines.append("look presets: %s" % ", ".join(sorted(tier3.PRESETS)))
    lines.append("scout <what>: %s, motif <pattern such as CxxC>"
                 % ", ".join(sorted(tools._scout_words())))
    lines.append("focus <place>: a number from the list `focus` prints, active site, "
                 "binding sites, variants, domain, ligand, a bound molecule's name or "
                 "code (heme, ZN), chain B, a residue (A/57), all")
    lines.append("search filters after a name (grab, find): organism=human  method=xray|"
                 "em|nmr  res<2.0  ligand=ATP  year>2015  pick=2  limit=5  mode=text")
    lines.append("")
    lines.append("PYMOL COMMANDS YOU MAY ALSO USE (display only)")
    lines.append("hide <representation>, <selection>     e.g.  hide everything, solvent")
    lines.append("show <representation>, <selection>     e.g.  show sticks, resn HIS")
    lines.append("color <colour>, <selection>            e.g.  color red, ss h")
    lines.append("set <setting>, <value>                 only: %s" % ", ".join(ai._SETTINGS))
    lines.append("distance <name>, <sel1>, <sel2>[, mode=2 for hydrogen bonds]")
    lines.append("zoom <selection>   orient <selection>   center <selection>")
    lines.append("turn <x|y|z>, <degrees>   rock   bg_color <colour>")
    lines.append("Selections use PyMOL's language: chain A, resi 50-60, resn HIS, ss h "
                 "(helix), ss s (strand), organic, inorganic, solvent, polymer, "
                 "name CA, <object name>, and/or/not, within 5 of ...")
    _reference[0] = "\n".join(lines)
    return _reference[0]


# ------------------------------------------------------------------- context
def context():
    """A few lines on what is loaded.  Call on the GUI thread."""
    from . import places
    names = base._molecules()
    if not names:
        lines = ["Nothing is loaded."]
    else:
        lines = []
        for name in names[:6]:
            identity = base._objects.get(name)
            what = "%s %s" % (identity[0].upper(), identity[1]) if identity else "unknown source"
            aliases = sorted(a for a, o in base._aliases.items() if o == name and a != name.lower())
            chains = cmd.get_chains("%s and polymer" % name)
            codes = set()
            cmd.iterate("%s and not polymer and not solvent" % name, "codes.add(resn)",
                        space={"codes": codes})
            count = cmd.count_atoms("%s and polymer and name CA" % name)
            line = "object %s (%s): %d residues, chains %s" % (
                name, what, count, ",".join(c or "-" for c in chains[:12]) or "none")
            if aliases:
                line += "; also called %s" % ", ".join(aliases[:4])
            if codes:
                line += "; bound molecules %s" % ", ".join(sorted(codes)[:12])
            if not cmd.get_object_list("(%s) and enabled" % name):
                line += "; hidden"
            lines.append(line)
            listed = places._places.get(name)
            if listed:
                lines.append("  key places of %s: %s" % (name, "; ".join(
                    "%d %s" % (i, p.title) for i, p in enumerate(listed[:14], 1))))
        if len(names) > 6:
            lines.append("... and %d more objects" % (len(names) - 6))
        if base._last[0] in names:
            lines.append("The structure worked on last: %s" % base._last[0])
    recent = [c for c in base._recent if not c.startswith(("ask ", "askai "))][-6:]
    if recent:
        lines.append("Commands run recently, oldest first: %s" % " | ".join(recent))
    return "\n".join(lines)


def call(sentence, reference_text, context_text):
    """Worker-thread half: returns ("plan", Plan) or ("error", AIError)."""
    try:
        return "plan", ai.ask(sentence, reference_text, context_text)
    except ai.AIError as error:
        return "error", error


# ------------------------------------------------------------------ settings
_ABOUT = """AI FALLBACK FOR  ask
   When the built-in rules do not understand a sentence,  ask  can send it to
   an AI model (Claude, by Anthropic) that picks Prospector commands for you.
   The rules always go first; the AI is only asked when they give up.

 WHAT IS SENT   your sentence, the list of commands, and the names of what
                is loaded (structures, chains, bound molecules, recent
                commands). Never files or structures.
 WHAT IT COSTS  roughly 1 cent per question with the default model, billed
                to your own Anthropic account. Questions the rules answer
                cost nothing.
 SAFETY         the AI only proposes command lines. Prospector runs them only
                if they are its own commands or simple display commands, and
                always shows them to you.

 SET UP         1. Create a key at  console.anthropic.com  (API keys)
                2. prospector ai key      opens a box to paste it into
 COMMANDS       askai <sentence>          skip the rules, go straight to the AI
                prospector ai on | off    use the AI, or never use it
                prospector ai test        one small question, to check the key
                prospector ai model <name>   default %s
                prospector ai forget      delete the stored key
                prospector ai learned     sentences the AI answered once; ask now
                                          runs them without it (unlearn = forget)"""


def _ask_for_key():
    try:
        from pymol.Qt import QtWidgets
        if QtWidgets.QApplication.instance() is None:
            raise RuntimeError("no window")
        text, ok = QtWidgets.QInputDialog.getText(
            None, "Prospector - AI key",
            "Paste your Anthropic API key (it starts with sk-ant-).\n"
            "It is stored on this computer only, in the Prospector folder.",
            QtWidgets.QLineEdit.Password)
        return (text if ok else ""), True
    except Exception:
        return "", False


def settings(words):
    """`prospector ai ...` - runs on the GUI thread."""
    verb = words[0] if words else ""
    base._header("prospector ai %s" % " ".join(words))
    _say("")
    if verb in ("on", "off"):
        ai.set_enabled(verb == "on")
        if verb == "off":
            _say("AI fallback is OFF.  ask  uses only the built-in rules; nothing is sent.")
        elif ai.api_key():
            _say("AI fallback is ON.  ask  consults the AI when the rules give up.")
        else:
            _say("AI fallback is ON, but no key is stored yet. Add one:  prospector ai key")
    elif verb == "key":
        if len(words) > 1:
            _say("Do not type the key on the command line: it would be saved in PyMOL's "
                 "log. Type just  prospector ai key  and paste it into the box.")
            _say("If that key was real, consider replacing it in the Anthropic Console.")
        else:
            text, had_window = _ask_for_key()
            if not had_window:
                _say("No window is available to paste the key into. Set the environment "
                     "variable ANTHROPIC_API_KEY before starting PyMOL instead.")
            elif not text.strip():
                _say("No key was entered; nothing changed.")
            else:
                ai.set_key(text)
                ai.set_enabled(True)
                _say("KEY STORED   %s  in %s (readable by you only)"
                     % (ai.masked(ai.api_key()), ai._key_path()))
                _say("CHECK IT     prospector ai test")
    elif verb == "learned":
        learned = ai.learned_all()
        if not learned:
            _say("Nothing learned yet. When the AI answers a sentence the rules did not "
                 "understand,  ask  remembers it and runs it without the AI next time.")
        else:
            _say("LEARNED FROM THE AI   %d sentence%s (kept in %s)"
                 % (len(learned), "" if len(learned) == 1 else "s", ai._learned_path()))
            for entry in list(learned.values())[-15:]:
                _say("   %s" % entry.get("sentence", "")[:70])
                _say("      ->  %s" % " ; ".join(c for c, _ in entry.get("steps", []))[:70])
            if len(learned) > 15:
                _say("   ... and %d older ones" % (len(learned) - 15))
            _say("FORGET THEM  prospector ai unlearn")
    elif verb == "unlearn":
        count = ai.learned_forget()
        _say("Forgot %d learned sentence%s." % (count, "" if count == 1 else "s"))
    elif verb == "forget":
        _say("The stored key was deleted." if ai.forget_key()
             else "There was no stored key.")
        if ai.key_source() == "environment":
            _say("(A key is still set in the environment variable ANTHROPIC_API_KEY.)")
    elif verb == "model":
        if len(words) < 2:
            _say("The model is  %s . Change it with:  prospector ai model <name>"
                 % ai.model())
            _say("   %s   best understanding (default)" % ai.DEFAULT_MODEL)
            _say("   claude-haiku-5-5  much cheaper and faster, a little less careful")
        else:
            try:
                ai.set_model(words[1])
                _say("The AI model is now  %s ." % ai.model())
            except ai.AIError as error:
                _say(str(error))
    elif verb == "test":
        if not ai.api_key():
            _say("No key is stored. Add one first:  prospector ai key")
        else:
            _say("Sending one small question to the AI ...")
            reference_text, context_text = reference(), context()

            def done(result, error):
                if error:
                    base._unexpected(error)
                    return
                kind, value = result
                if kind == "error":
                    _say("NOT WORKING  %s" % value)
                else:
                    good, bad = ai.validate(value)
                    usage = value.usage or {}
                    _say("WORKING      model %s answered:  %s" % (
                        value.model or ai.model(),
                        " ; ".join(c for c, _ in good) or value.cannot or "(nothing)"))
                    _say("             tokens: %s in (%s read from cache), %s out"
                         % (usage.get("input_tokens", "?"),
                            usage.get("cache_read_input_tokens", 0),
                            usage.get("output_tokens", "?")))
                _say(report.LINE)

            with worker.next_in_line():
                worker.submit(lambda: call("show the protein as sticks and colour the "
                                           "helices blue", reference_text, context_text),
                              done)
            return
    else:
        key = ai.api_key()
        state = "OFF (prospector ai on  switches it on)" if not ai.switched_on() else (
            "ON" if key else "waiting for a key (prospector ai key)")
        _say("STATUS       %s" % state)
        _say("KEY          %s" % ("%s, %s" % (ai.masked(key), "stored in the Prospector "
                                              "folder" if ai.key_source() == "stored"
                                              else "from the environment variable "
                                              "ANTHROPIC_API_KEY") if key else "none"))
        _say("MODEL        %s" % ai.model())
        _say("")
        _say(_ABOUT % ai.DEFAULT_MODEL)
    _say(report.LINE)
