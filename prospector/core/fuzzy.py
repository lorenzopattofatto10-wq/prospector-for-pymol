"""Spelling help: British/American variants and typo suggestions.

Works fully offline, using a bundled vocabulary (core/data/names.txt) built
from UniProt protein names and human gene symbols.
"""
from __future__ import annotations

import difflib
import os
import re

# British / older spellings -> the spelling UniProt and the PDB use.
_VARIANTS = [
    (r"haem", "hem"),            # haemoglobin, haemagglutinin, haem
    (r"oestr", "estr"),          # oestrogen
    (r"foet", "fet"),            # foetal
    (r"tumour", "tumor"),
    (r"sulph", "sulf"),          # sulphatase
    (r"leuco", "leuko"),         # leucocyte -> leukocyte (leucine is untouched)
    (r"caesium", "cesium"),
    (r"aluminium", "aluminum"),
    (r"signalling", "signaling"),
    (r"fibre", "fiber"),
    (r"coeliac", "celiac"),
    (r"anaemi", "anemi"),
    (r"glycaem", "glycem"),
    (r"paed", "ped"),
]

_vocab = None


def _load():
    global _vocab
    if _vocab is None:
        words, genes = set(), set()
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "names.txt")
        section = None
        try:
            with open(path, "r", encoding="utf-8") as handle:
                for line in handle:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    if line.startswith("["):
                        section = line
                    elif section == "[words]":
                        words.add(line.lower())
                    elif section == "[genes]":
                        genes.add(line.lower())
        except OSError:
            pass
        _vocab = (words, genes, sorted(words))
    return _vocab


def normalize_spelling(text):
    """Return (text_in_US_spelling, list_of_notes)."""
    result = text
    for pattern, replacement in _VARIANTS:
        result = re.sub(pattern, replacement, result, flags=re.IGNORECASE)
    notes = []
    if result.lower() != text.lower():
        notes.append("Reading '%s' as '%s' (the spelling used by UniProt and the PDB)."
                     % (text, result))
    return result, notes


def known(word):
    words, genes, _ = _load()
    return word.lower() in words or word.lower() in genes


def suggest(text):
    """Suggest a corrected query for a probable typo, or return None.

    Only words that are NOT in the vocabulary are changed, and only when a
    close match exists - so 'asdfqwerty' gets no suggestion at all."""
    words, genes, word_list = _load()
    changed = False
    out = []
    for token in text.split():
        low = token.lower()
        if len(low) < 4 or low in words or low in genes or not low.isalpha():
            out.append(token)
            continue
        close = difflib.get_close_matches(low, word_list, n=1, cutoff=0.8)
        if close and close[0] != low:
            out.append(close[0])
            changed = True
        else:
            out.append(token)
    return " ".join(out) if changed else None
