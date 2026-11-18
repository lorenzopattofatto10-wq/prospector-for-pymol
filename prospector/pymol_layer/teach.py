"""Teach mode: show the raw PyMOL commands behind every Prospector action.

Commands call ``note("show sticks, organic", "draw ligands as sticks")`` while
they work; when the command finishes, the collected lines are printed if
teach mode is on.  The setting is remembered between PyMOL sessions.
"""
from __future__ import annotations

from ..core import store

_notes = []
_hinted = False


def enabled():
    return bool(store.get_setting("teach", False))


def set_enabled(value):
    store.set_setting("teach", bool(value))


def note(command, comment=""):
    _notes.append((command, comment))


def clear():
    del _notes[:]


def take():
    """Lines to print for the command that just finished ([] when there are none)."""
    global _hinted
    notes = list(_notes)
    clear()
    if not notes:
        return []
    if not enabled():
        if _hinted:
            return []
        _hinted = True
        return ["Tip:  prospector teach on  shows the PyMOL commands behind every step."]
    notes = [(c if len(c) <= 110 else c[:106] + " ...", comment) for c, comment in notes]
    width = min(52, max(len(command) for command, _ in notes))
    lines = ["", "Teach mode - the PyMOL commands behind this (you can type them yourself):"]
    for command, comment in notes:
        if comment:
            lines.append("  %-*s  # %s" % (width, command, comment))
        else:
            lines.append("  " + command)
    return lines
