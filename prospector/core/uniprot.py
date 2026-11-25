"""UniProt REST client (https://rest.uniprot.org).

Endpoint and field names checked against the live service:
  GET /uniprotkb/search?query=...&fields=...&size=...
  query fields : gene_exact, protein_name, organism_id, organism_name,
                 reviewed, structure_3d, accession
  return fields: accession, id, gene_primary, gene_synonym, protein_name,
                 organism_name, organism_id, length, reviewed, xref_pdb,
                 annotation_score
"""
from __future__ import annotations

import re
import urllib.parse

from . import cache, http

BASE = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = ("accession,id,gene_primary,gene_synonym,protein_name,organism_name,"
          "organism_id,length,reviewed,xref_pdb,annotation_score")
ACCESSION_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$", re.I)


class PdbRef(object):
    """One PDB cross-reference from a UniProt entry."""
    __slots__ = ("pdb_id", "method", "resolution", "ranges")

    def __init__(self, pdb_id, method, resolution, ranges):
        self.pdb_id = pdb_id
        self.method = method          # "X-ray", "EM", "NMR", ...
        self.resolution = resolution  # float or None
        self.ranges = ranges          # list of (start, end) in UniProt numbering

    def covered(self):
        """Number of UniProt residues covered by this entry (union of ranges)."""
        residues = 0
        last = 0
        for start, end in sorted(self.ranges):
            start = max(start, last + 1)
            if end >= start:
                residues += end - start + 1
                last = end
        return residues


class Protein(object):
    """The parts of a UniProt entry Prospector needs."""

    def __init__(self, record):
        self.accession = record.get("primaryAccession", "")
        self.entry_name = record.get("uniProtkbId", "")
        self.reviewed = "unreviewed" not in record.get("entryType", "").lower()
        self.annotation_score = record.get("annotationScore", 0) or 0
        organism = record.get("organism", {})
        self.organism = organism.get("scientificName", "")
        self.organism_common = organism.get("commonName", "")
        self.taxid = organism.get("taxonId")
        self.length = (record.get("sequence") or {}).get("length", 0) or 0

        self.gene = ""
        self.gene_synonyms = []
        for gene in record.get("genes", []) or []:
            if "geneName" in gene and not self.gene:
                self.gene = gene["geneName"].get("value", "")
            for synonym in gene.get("synonyms", []) or []:
                self.gene_synonyms.append(synonym.get("value", ""))
            if "geneName" in gene and gene["geneName"].get("value") != self.gene:
                self.gene_synonyms.append(gene["geneName"].get("value", ""))

        self.name = ""
        self.other_names = []
        description = record.get("proteinDescription", {}) or {}
        self._read_names(description, top=True)
        for part in (description.get("contains", []) or []) + (description.get("includes", []) or []):
            self._read_names(part, top=False)
        if not self.name:
            self.name = self.entry_name
        # only present when asked for (FAMILY_FIELDS): what UniProt files it under
        self.keywords = [k.get("name", "") for k in record.get("keywords", []) or []]
        self.families = []
        for comment in record.get("comments", []) or []:
            if comment.get("commentType") == "SIMILARITY":
                for text in comment.get("texts", []) or []:
                    self.families.extend(part.strip() for part in
                                         re.split(r"\.\s*", text.get("value", ""))
                                         if part.strip())

        self.pdb = {}
        for xref in record.get("uniProtKBCrossReferences", []) or []:
            if xref.get("database") != "PDB":
                continue
            props = dict((p.get("key"), p.get("value")) for p in xref.get("properties", []))
            resolution = None
            match = re.match(r"([0-9.]+)", props.get("Resolution", "") or "")
            if match:
                try:
                    resolution = float(match.group(1))
                except ValueError:
                    pass
            ranges = [(int(a), int(b)) for a, b in
                      re.findall(r"=(\d+)-(\d+)", props.get("Chains", "") or "")]
            pdb_id = xref.get("id", "").upper()
            self.pdb[pdb_id] = PdbRef(pdb_id, props.get("Method", ""), resolution, ranges)

    def _read_names(self, block, top):
        def value(item):
            return ((item or {}).get("fullName") or {}).get("value", "")

        recommended = value(block.get("recommendedName"))
        if recommended:
            if top and not self.name:
                self.name = recommended
            else:
                self.other_names.append(recommended)
        for short in (block.get("recommendedName") or {}).get("shortNames", []) or []:
            self.other_names.append(short.get("value", ""))
        for alt in block.get("alternativeNames", []) or []:
            if value(alt):
                self.other_names.append(value(alt))
            for short in alt.get("shortNames", []) or []:
                self.other_names.append(short.get("value", ""))
        for sub in block.get("submissionNames", []) or []:
            if top and not self.name:
                self.name = value(sub)

    def all_names(self):
        return [self.name] + [n for n in self.other_names if n]

    def label(self):
        gene = " (%s)" % self.gene if self.gene else ""
        return "%s%s" % (self.name, gene)

    def __repr__(self):
        return "Protein(%s %s %s pdb=%d)" % (self.accession, self.gene, self.organism,
                                             len(self.pdb))


def _escape(text):
    """Make user text safe inside a UniProt (Lucene) query."""
    return re.sub(r'(["\\])', r"\\\1", text)


class Found(list):
    """A list of Protein objects plus ``total``: how many UniProt matched overall."""
    total = 0


def _run(query, size):
    url = BASE + "?" + urllib.parse.urlencode({"query": query, "fields": FIELDS, "size": size})
    data, total = http.request_json(url, ttl=7 * cache.DAY, header="x-total-results")
    found = Found(Protein(record) for record in (data or {}).get("results", []))
    try:
        found.total = int(total)
    except (TypeError, ValueError):
        found.total = len(found)
    return found


def organism_clause(organism):
    if organism is None:
        return ""
    if organism.taxid:
        return " AND (organism_id:%d)" % organism.taxid
    return ' AND (organism_name:"%s")' % _escape(organism.name)


def search_name(name, organism=None, extra="", size=25, reviewed=True):
    """Proteins whose gene symbol equals, or whose protein name contains, ``name``."""
    safe = _escape(name)
    parts = ['(protein_name:"%s")' % safe]
    if " " not in name.strip():
        parts.insert(0, '(gene_exact:"%s")' % safe)
    query = "(%s)" % " OR ".join(parts)
    if reviewed:
        query += " AND (reviewed:true)"
    query += organism_clause(organism) + extra
    return _run(query, size)


def search_words(words, organism=None, size=25):
    """Looser search: every word must appear in the protein name."""
    query = " AND ".join("(protein_name:%s)" % _escape(w) for w in words)
    query = "(%s) AND (reviewed:true)%s" % (query, organism_clause(organism))
    return _run(query, size)


FAMILY_FIELDS = FIELDS + ",keyword,protein_families"


def search_family(phrases, organism=None, size=500, pages=2):
    """Reviewed proteins with a 3D structure whose name, UniProt keyword or
    protein family mentions any of ``phrases`` ("kinase", "G protein-coupled
    receptor").  Loose on purpose: resolve.family_member decides who belongs.
    Up to ``pages`` pages of ``size`` proteins."""
    parts = []
    for phrase in phrases:
        safe = _escape(phrase)
        parts.extend(['(protein_name:"%s")' % safe, '(keyword:"%s")' % safe,
                      '(family:"%s")' % safe])
    query = "(%s) AND (reviewed:true) AND (structure_3d:true)%s" % (
        " OR ".join(parts), organism_clause(organism))
    url = BASE + "?" + urllib.parse.urlencode({"query": query, "fields": FAMILY_FIELDS,
                                               "size": size})
    found = Found()
    for _ in range(pages):
        data, link = http.request_json(url, ttl=7 * cache.DAY, header="link")
        found.extend(Protein(record) for record in (data or {}).get("results", []))
        following = re.search(r'<([^>]+)>;\s*rel="next"', link or "")
        if not following:
            break
        url = following.group(1)
    found.total = len(found)
    return found


def by_accession(accession):
    found = _run("(accession:%s)" % _escape(accession.upper()), 1)
    return found[0] if found else None


DETAIL_FIELDS = ("accession,id,gene_primary,gene_synonym,protein_name,organism_name,"
                 "organism_id,length,sequence,xref_pdb,cc_function,ft_act_site,ft_binding,ft_site,"
                 "ft_domain,ft_dna_bind,ft_zn_fing,ft_disulfid,ft_mod_res,ft_carbohyd,"
                 "ft_lipid,ft_crosslnk,ft_variant")


def details(accession):
    """One UniProt entry with function text, sequence and feature annotations.
    Returns (Protein, raw record) or (None, None)."""
    url = BASE + "?" + urllib.parse.urlencode(
        {"query": "(accession:%s)" % _escape(accession.upper()), "fields": DETAIL_FIELDS,
         "size": 1})
    data = http.request_json(url, ttl=7 * cache.DAY) or {}
    records = data.get("results") or []
    if not records:
        return None, None
    return Protein(records[0]), records[0]


def function_text(record):
    """The 'Function' comment of a raw UniProt record, as plain text."""
    for comment in record.get("comments") or []:
        if comment.get("commentType") == "FUNCTION":
            texts = [t.get("value", "") for t in comment.get("texts") or []]
            text = " ".join(t for t in texts if t)
            # drop evidence tags such as (PubMed:12345) for readability
            return re.sub(r"\s*\((?:PubMed|Ref\.|By similarity)[^)]*\)", "", text).strip()
    return ""


def lengths(accessions):
    """Sequence lengths for many accessions in one request: {accession: length}."""
    result = {}
    accessions = sorted(set(a for a in accessions if ACCESSION_RE.match(a)))
    for start in range(0, len(accessions), 80):
        chunk = accessions[start:start + 80]
        query = " OR ".join("(accession:%s)" % a for a in chunk)
        url = BASE + "?" + urllib.parse.urlencode(
            {"query": query, "fields": "accession,length", "size": len(chunk)})
        data = http.request_json(url, ttl=30 * cache.DAY) or {}
        for record in data.get("results", []):
            result[record.get("primaryAccession")] = (record.get("sequence") or {}).get("length", 0)
    return result
