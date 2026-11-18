"""Persistent search history and favourites (JSON files in ~/.prospector)."""
from __future__ import annotations

import json
import os
import time

from . import cache

HISTORY_LIMIT = 200


def _path(name):
    folder = cache.home()
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, name)


def _read(name):
    try:
        with open(_path(name), "r", encoding="utf-8") as handle:
            data = json.load(handle)
        return data if isinstance(data, list) else []
    except (OSError, ValueError):
        return []


def _write(name, data):
    tmp = _path(name) + ".tmp"
    try:
        with open(tmp, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=1)
        os.replace(tmp, _path(name))
        return True
    except OSError:
        return False


# ------------------------------------------------------------------ history
def add_history(command, query, outcome, top=""):
    """Remember one search.  ``outcome`` is the result kind, ``top`` the best hit."""
    if not query.strip():
        return
    items = _read("history.json")
    items.append({"t": time.time(), "command": command, "query": query.strip(),
                  "outcome": outcome, "top": top})
    _write("history.json", items[-HISTORY_LIMIT:])


def history(count=None):
    items = _read("history.json")
    return items[-count:] if count else items


def clear_history():
    count = len(_read("history.json"))
    _write("history.json", [])
    return count


# --------------------------------------------------------------- favourites
def favorites():
    return _read("favorites.json")


def add_favorite(identifier, title="", note=""):
    """Returns True if added, False if it was already a favourite (note updated)."""
    items = favorites()
    for item in items:
        if item.get("id", "").upper() == identifier.upper():
            if note:
                item["note"] = note
            if title and not item.get("title"):
                item["title"] = title
            _write("favorites.json", items)
            return False
    items.append({"id": identifier, "title": title, "note": note, "t": time.time()})
    _write("favorites.json", items)
    return True


def remove_favorite(identifier):
    items = favorites()
    kept = [i for i in items if i.get("id", "").upper() != identifier.upper()]
    if len(kept) == len(items):
        return False
    _write("favorites.json", kept)
    return True


def when(timestamp):
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(timestamp or 0))


# ----------------------------------------------------------------- settings
def get_setting(name, default=None):
    try:
        with open(_path("settings.json"), "r", encoding="utf-8") as handle:
            data = json.load(handle)
        return data.get(name, default) if isinstance(data, dict) else default
    except (OSError, ValueError):
        return default


def set_setting(name, value):
    try:
        with open(_path("settings.json"), "r", encoding="utf-8") as handle:
            data = json.load(handle)
        if not isinstance(data, dict):
            data = {}
    except (OSError, ValueError):
        data = {}
    data[name] = value
    try:
        with open(_path("settings.json"), "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=1)
    except OSError:
        pass
