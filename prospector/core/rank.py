"""Scoring: why one structure ranks above another.

Each hit gets 0-100 points from six components (weights below), plus a small
bonus if it is on Prospector's hand-curated list of classic structures.
Change WEIGHTS to change the ranking - nothing else depends on the numbers.
"""
from __future__ import annotations

import json
import os
import re
import textwrap
import time

from . import fuzzy, organisms

WEIGHTS = {
    "relevance": 30,    # right protein, title matches, natural (unmutated) sequence
    "coverage": 20,     # fraction of the UniProt sequence present in the structure
    "resolution": 20,   # 1.0 A or better = full marks, 4.0 A or worse = none
    "organism": 15,     # your organism filter; otherwise human > model organisms > other
    "method": 10,       # X-ray / cryo-EM above NMR (unless you filter by method)
    "recency": 5,       # deliberately small, so classics are not buried
}
CLASSIC_BONUS = 8
ORDER = ["relevance", "coverage", "resolution", "organism", "method", "recency"]

_classics = None


def classics():
    """{PDB id: short reason}: Prospector's curated list of textbook structures."""
    global _classics
    if _classics is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "classics.json")
        try:
            with open(path, "r", encoding="utf-8") as handle:
                _classics = json.load(handle)
        except (OSError, ValueError):
            _classics = {}
    return _classics


def resolution_score(resolution):
    if resolution is None:
        return 0.5              # NMR etc.: neutral rather than punished
    return max(0.0, min(1.0, (4.0 - resolution) / 3.0))


def method_score(methods):
    best = 0.5
    for method in methods:
        upper = (method or "").upper()
        if "X-RAY" in upper or upper == "X-RAY":
            best = max(best, 1.0)
        elif "MICROSCOPY" in upper or upper == "EM" or "NEUTRON" in upper:
            best = max(best, 0.9)
        elif "NMR" in upper:
            best = max(best, 0.6)
    return best


def recency_score(year, now=None):
    now = now or time.localtime().tm_year
    if not year:
        return 0.0
    return max(0.0, min(1.0, (year - 1990) / float(max(1, now - 1990))))


def _tokens(text):
    text, _ = fuzzy.normalize_spelling(text)
    return set(re.findall(r"[a-z0-9]+", text.lower()))


class Hit(object):
    def __init__(self, entry):
        self.entry = entry
        self.parts = {}            # component -> points earned
        self.score = 0.0
        self.coverage = None       # 0..1 or None when unknown
        self.coverage_of = ""      # gene/accession the coverage refers to
        self.accession = ""        # UniProt accession of the chain that was scored
        self.mutations = None      # engineered mutations in that chain (None = unknown)
        self.partners = 0          # other, unrelated protein chains in the entry
        self.textbook = False      # a classic structure of exactly what was searched
        self.classic = ""          # reason text if this is a classic structure

    @property
    def pdb_id(self):
        return self.entry.pdb_id


class Context(object):
    """Everything the scorer needs to know about the search."""

    def __init__(self, query_text, targets=None, lengths=None, text_scores=None,
                 organism_filter=False, extra_words=()):
        self.words = [w for w in _tokens(query_text) if w not in ("protein", "the", "of")]
        self.extra_words = set()
        for text in extra_words:
            self.extra_words |= _tokens(text)
        self.targets = targets or {}        # {accession: uniprot.Protein}
        self.lengths = lengths or {}        # {accession: sequence length}
        self.text_scores = text_scores or {}  # {pdb id: RCSB relevance 0..1}
        self.organism_filter = organism_filter


def score(entry, ctx):
    hit = Hit(entry)
    text_words = _tokens(entry.title + " " + " ".join(entry.descriptions))
    if ctx.words:
        matched = sum(1 for w in ctx.words if w in text_words) / float(len(ctx.words))
    else:
        matched = 0.0
    if ctx.extra_words and matched < 1.0 and (ctx.extra_words & text_words):
        matched = max(matched, 0.8)   # title uses the gene / official name instead

    # coverage: best over the proteins we are looking for (or any, in text mode)
    best_cov, best_acc = None, ""
    accessions = ctx.targets if ctx.targets else entry.coverage
    for accession in accessions:
        covered = entry.coverage.get(accession)
        length = (ctx.targets[accession].length if accession in ctx.targets
                  else ctx.lengths.get(accession))
        if covered and length:
            # the deposited sequence can be much longer than what was actually
            # built into the model, so scale by the fraction that was modelled
            fraction = min(1.0, covered / float(length)) * entry.modelled_fraction
            if best_cov is None or fraction > best_cov:
                best_cov, best_acc = fraction, accession
    hit.coverage = best_cov
    if best_acc:
        target = ctx.targets.get(best_acc)
        hit.coverage_of = (target.gene or best_acc) if target else best_acc
        hit.accession = best_acc
        hit.mutations = entry.mutations.get(best_acc)
    natural = 1.0 if hit.mutations == 0 else (0.5 if hit.mutations is None else 0.0)

    if ctx.targets:
        hit.partners = sum(1 for accessions in entry.protein_entities
                           if not any(a in ctx.targets for a in accessions))
        alone = 1.0 if hit.partners == 0 else 0.0
        # identity is confirmed by the UniProt mapping (60%); the rest rewards a
        # title that names the protein, a sequence without engineered mutations,
        # and a structure of the protein on its own (no antibodies, inhibitor
        # proteins or other partners)
        relevance = 0.6 + 0.15 * matched + 0.15 * natural + 0.10 * alone
    else:
        relevance = (0.3 * ctx.text_scores.get(entry.pdb_id, 0.5) + 0.5 * matched
                     + 0.2 * natural)

    if ctx.organism_filter:
        organism = 1.0
    elif ctx.targets and best_acc in ctx.targets:
        organism = organisms.preference(ctx.targets[best_acc].taxid)
    else:
        organism = max([organisms.preference(tax, name) for name, tax in entry.organisms]
                       or [0.4])

    values = {
        "relevance": relevance,
        "coverage": 0.5 if best_cov is None else best_cov,
        "resolution": resolution_score(entry.resolution),
        "organism": organism,
        "method": method_score(entry.methods),
        "recency": recency_score(entry.year),
    }
    hit.parts = dict((k, WEIGHTS[k] * v) for k, v in values.items())
    total = sum(hit.parts.values())
    reason = classics().get(entry.pdb_id)
    if reason:
        hit.classic = reason
        hit.parts["classic"] = CLASSIC_BONUS
        total += CLASSIC_BONUS
    hit.score = min(100.0, total)
    return hit


def rank(entries, ctx):
    hits = [score(entry, ctx) for entry in entries]
    hits.sort(key=lambda h: (-h.score, h.entry.resolution or 99, h.pdb_id))
    return hits


def _quality(entry):
    """Plain words for the experimental quality of an entry."""
    method = entry.method_short()
    if entry.resolution is None:
        return "%s structure (no resolution value)" % method
    res = entry.resolution
    word = ("very high" if res < 1.5 else "high" if res < 2.2 else
            "medium" if res < 3.0 else "low")
    return "%s resolution (%.2f A, %s)" % (word, res, method)


def explain(hits):
    """Why #1 won: plain reasons first, then the points."""
    if not hits:
        return ""
    top = hits[0]
    entry = top.entry
    reasons = ["from %s" % entry.organism_text(top.accession), _quality(entry)]
    if top.mutations == 0:
        reasons.append("natural sequence (no engineered mutations)")
    elif top.mutations:
        reasons.append("%d engineered mutation%s" % (top.mutations,
                                                     "" if top.mutations == 1 else "s"))
    if top.coverage is not None:
        reasons.append("%.0f%% of %s is actually modelled" % (100 * top.coverage,
                                                              top.coverage_of or "the protein"))
    if top.partners:
        reasons.append("in complex with %d other protein%s"
                       % (top.partners, "" if top.partners == 1 else "s"))
    elif top.accession and top.coverage_of:
        reasons.append("the protein on its own")
    if top.classic:
        reasons.append("a classic textbook structure (%s)" % top.classic)
    if top.textbook:
        sentence = "It is the textbook structure of what you searched for (%s)." % top.classic
        lines = ["WHY %s IS #1" % top.pdb_id,
                 textwrap.fill(sentence, 78, initial_indent="   ", subsequent_indent="   "),
                 "   Textbook structures are listed first; the rest follow by score."]
        return "\n".join(lines)
    def wrap(text):
        return textwrap.fill(text, 78, initial_indent="   ", subsequent_indent="   ")

    sentence = "; ".join(reasons) + "."
    lines = ["WHY %s IS #1" % top.pdb_id, wrap(sentence[:1].upper() + sentence[1:])]
    if len(hits) > 1:
        second = hits[1]
        gaps = [(top.parts.get(k, 0) - second.parts.get(k, 0), k)
                for k in ORDER + ["classic"]]
        gap, component = max(gaps)
        margin = top.score - second.score
        if margin < 0.05 and gap <= 0.05:
            lines.append("   It ties with #2 (%s); the better resolution decides."
                         % second.pdb_id)
        else:
            label = "being a classic structure" if component == "classic" else component
            lines.append("   It is %.1f points ahead of #2 (%s), mainly on %s."
                         % (margin, second.pdb_id, label))
    points = " | ".join("%s %.0f/%d" % (k, top.parts[k], WEIGHTS[k]) for k in ORDER)
    if top.classic:
        points += " | classic +%d" % CLASSIC_BONUS
    lines.append(textwrap.fill("Points: %s = %.0f/100" % (points, top.score), 78,
                               initial_indent="   ", subsequent_indent="           "))
    return "\n".join(lines)
