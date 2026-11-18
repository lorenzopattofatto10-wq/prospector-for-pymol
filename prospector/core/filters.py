"""Parse the text after a command: free-text query plus filters.

    hemoglobin organism=human method=xray res<2.0 ligand=HEM year>2015 limit=5
"""
from __future__ import annotations

import difflib
import re

from . import organisms

METHODS = {
    "xray": "X-RAY DIFFRACTION", "x-ray": "X-RAY DIFFRACTION", "xrd": "X-RAY DIFFRACTION",
    "crystal": "X-RAY DIFFRACTION",
    "em": "ELECTRON MICROSCOPY", "cryoem": "ELECTRON MICROSCOPY",
    "cryo-em": "ELECTRON MICROSCOPY",
    "nmr": "SOLUTION NMR",
    "neutron": "NEUTRON DIFFRACTION",
}
KEYS = ["organism", "method", "res", "ligand", "year", "limit", "mode", "pick", "kind"]
_ALIASES = {"resolution": "res", "org": "organism", "species": "organism", "lig": "ligand",
            "n": "limit", "top": "limit"}

_FILTER = re.compile(r"""^(?P<key>[A-Za-z_]+)\s*(?P<op><=|>=|=|<|>)\s*(?P<value>.+)$""")


class FilterError(ValueError):
    """A filter the user typed could not be understood."""


class Query(object):
    def __init__(self):
        self.text = ""            # the protein name part
        self.organism = None      # organisms.Organism
        self.organism_from_text = False
        self.method = None        # key of METHODS: xray | em | nmr | neutron
        self.res_max = None
        self.res_min = None
        self.ligands = []         # chemical component ids, upper case
        self.year_min = None      # inclusive
        self.year_max = None      # inclusive
        self.limit = 10
        self.mode = "auto"        # auto | text
        self.pick = 1             # which row of the ranked list to take (grab)
        self.kind = "auto"        # auto | protein | molecule
        self.raw = ""

    def filter_words(self):
        """The filters, written back the way the user would type them."""
        out = []
        if self.organism and not self.organism_from_text:
            typed = self.organism.typed
            out.append("organism=%s" % (typed.replace(" ", "_")))
        if self.method:
            out.append("method=%s" % self.method)
        if self.res_max is not None:
            out.append("res<%g" % self.res_max)
        if self.res_min is not None:
            out.append("res>%g" % self.res_min)
        for ligand in self.ligands:
            out.append("ligand=%s" % ligand)
        if self.year_min is not None:
            out.append("year>%d" % (self.year_min - 1))
        if self.year_max is not None:
            out.append("year<%d" % (self.year_max + 1))
        return out

    def has_filters(self):
        return bool(self.filter_words()) or (self.organism is not None)


def _split(raw):
    """Split on spaces, keeping "quoted strings" together."""
    return [m.group(0) for m in re.finditer(r'[^\s"\']*(?:"[^"]*"|\'[^\']*\')[^\s]*|\S+', raw)]


def _number(value, what, example):
    try:
        return float(value)
    except ValueError:
        raise FilterError("%s needs a number, e.g. %s (you wrote '%s')." % (what, example, value))


def parse(raw):
    query = Query()
    query.raw = raw.strip()
    words = []
    for token in _split(query.raw):
        match = _FILTER.match(token)
        if not match:
            words.append(token)
            continue
        key = match.group("key").lower()
        key = _ALIASES.get(key, key)
        op = match.group("op")
        value = match.group("value").strip().strip("\"'")
        if key not in KEYS:
            close = difflib.get_close_matches(key, KEYS, n=1, cutoff=0.6)
            hint = " Did you mean '%s'?" % close[0] if close else ""
            raise FilterError("Unknown filter '%s'.%s Available filters: organism=human  "
                              "method=xray|em|nmr  res<2.0  ligand=ATP  year>2015  limit=10  "
                              "pick=2"
                              % (match.group("key"), hint))
        if not value:
            raise FilterError("The filter '%s' has no value." % token)

        if key == "organism":
            if op != "=":
                raise FilterError("Write the organism filter as organism=human.")
            query.organism = organisms.lookup(value)
        elif key == "method":
            method = value.lower()
            if op != "=" or method not in METHODS:
                raise FilterError("method must be one of xray, em, nmr (you wrote '%s')."
                                  % value)
            canonical = METHODS[method]
            query.method = [k for k in ("xray", "em", "nmr", "neutron")
                            if METHODS[k] == canonical][0]
        elif key == "res":
            number = _number(value.rstrip("AaÅ"), "res", "res<2.0")
            if op in ("<", "<=", "="):
                query.res_max = number
            else:
                query.res_min = number
        elif key == "ligand":
            if op != "=":
                raise FilterError("Write the ligand filter as ligand=ATP "
                                  "(the 3-letter PDB ligand code).")
            for code in value.split("+"):
                code = code.strip().upper()
                if not re.match(r"^[A-Z0-9]{1,5}$", code):
                    raise FilterError(
                        "ligand= needs a PDB ligand code such as ATP, HEM or STI "
                        "(you wrote '%s'). Look codes up at rcsb.org/ligand." % value)
                query.ligands.append(code)
        elif key == "year":
            year = int(_number(value, "year", "year>2015"))
            if op == ">":
                query.year_min = year + 1
            elif op == ">=":
                query.year_min = year
            elif op == "<":
                query.year_max = year - 1
            elif op == "<=":
                query.year_max = year
            else:
                query.year_min = query.year_max = year
        elif key == "limit":
            limit = int(_number(value, "limit", "limit=10"))
            if op != "=" or not 1 <= limit <= 50:
                raise FilterError("limit must be between 1 and 50, e.g. limit=10.")
            query.limit = limit
        elif key == "pick":
            pick = int(_number(value, "pick", "pick=2"))
            if op != "=" or not 1 <= pick <= 50:
                raise FilterError("pick chooses a row of the list, e.g. pick=2 for the "
                                  "second-best structure.")
            query.pick = pick
        elif key == "kind":
            kinds = {"protein": "protein", "proteins": "protein", "structure": "protein",
                     "molecule": "molecule", "mol": "molecule", "ligand": "molecule",
                     "small": "molecule", "auto": "auto"}
            if op != "=" or value.lower() not in kinds:
                raise FilterError("kind must be 'protein' or 'molecule', e.g. "
                                  "find cholesterol kind=protein")
            query.kind = kinds[value.lower()]
        elif key == "mode":
            if value.lower() not in ("auto", "text"):
                raise FilterError("mode must be 'auto' or 'text'.")
            query.mode = value.lower()

    query.limit = max(query.limit, query.pick)
    query.text = " ".join(words).strip().strip("\"'")
    return query
