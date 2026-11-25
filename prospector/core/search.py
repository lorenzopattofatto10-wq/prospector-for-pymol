"""The `find` pipeline: text in, ranked answer out.  No PyMOL imports here."""
from __future__ import annotations

import re

from . import alphafold, chemcomp, filters, fuzzy, http, organisms, pubchem, rank, rcsb
from . import resolve, uniprot

CONCEPT_RATIO = 0.01
DETAIL_POOL = 40   # how many entries get full metadata before the final ranking


class Result(object):
    """What `find` produced.  ``kind`` is one of:
    results | ambiguous | suggest | none | nofilter | alphafold | molecule | error"""

    def __init__(self, kind, query=None):
        self.kind = kind
        self.query = query
        self.raw = query.raw if query is not None else ""
        self.notes = []          # lines explaining how the query was read
        self.hits = []           # rank.Hit, best first (already cut to the limit)
        self.total = 0           # entries that matched before the limit
        self.mode = ""           # "protein" or "text"
        self.resolution = None   # resolve.Resolution
        self.suggestion = ""     # corrected query for kind == "suggest"
        self.model = None        # alphafold.Model for kind == "alphafold"
        self.protein = None      # uniprot.Protein the AlphaFold model belongs to
        self.molecule = None     # chemcomp.Molecule for kind == "molecule"
        self.message = ""        # for error / none / nofilter
        self.offline = False

    @property
    def top(self):
        """The structure this search stands for: the best one, or the row the
        user asked for with pick=N."""
        if not self.hits:
            return None
        pick = self.query.pick if self.query is not None else 1
        return self.hits[min(pick, len(self.hits)) - 1]

    def command(self):
        """The exact PyMOL command that acts on this result."""
        if self.kind == "results" and self.hits:
            return "fetch %s" % self.top.pdb_id
        if self.kind == "molecule" and self.molecule:
            return "molecule %s" % (self.molecule.comp_id or self.molecule.name)
        if self.kind == "alphafold" and self.model:
            return "alphafold %s" % self.model.accession
        if self.kind == "suggest":
            return "find %s" % " ".join([self.suggestion] + self.query.filter_words())
        return ""


def find(raw):
    """Run a search.  Never raises for bad input or network trouble: those come
    back as Result(kind="error")."""
    http.take_notes()
    query = None
    try:
        query = filters.parse(raw)
        result = _find(query)
    except filters.FilterError as err:
        result = Result("error", query)
        result.message = str(err)
    except http.NetworkError as err:
        result = Result("error", query)
        result.offline = err.offline
        if err.offline:
            result.message = ("You seem to be offline: %s Prospector needs the internet for "
                              "new searches; searches you ran before still work from the "
                              "cache." % err)
        else:
            result.message = ("Network problem: %s I tried %d times. Please run the same "
                              "command again in a moment - whatever was already fetched "
                              "is kept, so the next try is quicker."
                              % (err, http.ATTEMPTS))
    result.raw = raw.strip()
    result.notes.extend(http.take_notes())
    return result


def _find(query):
    if not query.text:
        result = Result("error", query)
        result.message = ("Type a protein name, gene or PDB ID after the command - for "
                          "example:  grab hemoglobin   or   find p53 organism=human "
                          "res<2.0")
        return result

    notes = []
    name, spelling_notes = fuzzy.normalize_spelling(query.text)
    notes.extend(spelling_notes)

    # 0. the user asked for the small molecule of that name
    if query.kind == "molecule":
        molecule = chemcomp.lookup(name)
        if molecule is None or not molecule.exact:
            # not in the PDB's dictionary under that name: ask PubChem
            other = pubchem.lookup(name)
            if other is not None and (other.has_shape() or molecule is None):
                molecule = other
        return _molecule_result(query, name, notes, molecule)

    # 1. a PDB ID typed directly
    if " " not in name and rcsb.looks_like_pdb_id(name):
        found = rcsb.entries([name])
        if found:
            result = Result("results", query)
            result.mode = "id"
            result.notes = notes + ["Recognised %s as a PDB ID." % name.upper()]
            entry = list(found.values())[0]
            lengths = uniprot.lengths(entry.coverage)
            result.hits = rank.rank([entry], rank.Context("", lengths=lengths))
            result.total = 1
            return result

    # 2. an organism mentioned in the text ("spike protein sars-cov-2")
    if query.organism is None:
        rest, organism = organisms.extract(name)
        if organism is not None:
            name = rest
            query.organism = organism
            query.organism_from_text = True
            notes.append("Reading '%s' as the organism: %s." % (organism.typed, organism.name))

    # "a protein from e. coli": no name at all, only an organism
    # the whole text is an organism ("e. coli"): any protein from it
    if query.organism is None and re.search(r"[ .\-]", name.strip()):
        whole = organisms.lookup(name)
        if whole.taxid is not None:
            query.organism = whole
            query.organism_from_text = True
            notes.append("Reading '%s' as the organism: %s." % (name, whole.name))
            name = ""
    generic = " ".join(w for w in name.lower().split()
                       if w not in ("a", "an", "any", "some", "the", "from", "of", "in"))
    if generic in ("", "protein", "proteins", "structure", "structures") and (
            query.organism is not None or query.has_filters()):
        name = "protein"
        query.mode = "text"

    # 3. which protein is meant?
    if query.mode == "text":
        resolution = resolve.Resolution("broad", name)
    else:
        resolution = resolve.resolve(name, query.organism)
    early_notes = list(notes)            # how the query was read, before any search
    first = resolution
    notes.extend(resolution.notes)

    key = " ".join(re.findall(r"[a-z0-9]+", name.lower()))
    if key in resolve.COMPLEXES and resolution.kind == "ambiguous" and \
            not query.has_filters():
        # "nucleosome": one machine built from many proteins - its structures are
        # found by the word, not by choosing one of its parts
        notes.append("'%s' is a complex of several proteins, so I searched for "
                     "structures of the whole complex." % name)
        resolution = resolve.Resolution("broad", name)

    # A loose name match with very few structures, for a word that appears in
    # thousands of PDB entries, is a concept ("ribosome"), not one protein.
    if resolution.kind in ("protein", "ambiguous") and \
            resolution.top_tier <= resolve.TIER_PREFIX and resolution.leader_count > 0:
        scope = filters.Query()
        scope.organism = query.organism      # compare like with like
        mentions, _ = rcsb.search_text(name, scope, rows=1)
        if resolution.leader_count < CONCEPT_RATIO * mentions and \
                resolution.kind == "ambiguous" and not query.has_filters() and \
                _family_members(name, resolution.choices) >= 2:
            # "aquaporin": a family whose members are literally called Aquaporin-1,
            # Aquaporin Z ...  Offer them as a list instead of picking one.
            notes.append("'%s' names a family of proteins, so I won't pick one for you: "
                         "these are its members with the most structures." % name)
        elif resolution.leader_count < CONCEPT_RATIO * mentions:
            notes.append("'%s' is mentioned in %d PDB entries, but no single protein by "
                         "that name has more than %d structures, so I searched the text "
                         "of every PDB entry instead."
                         % (name, mentions, resolution.leader_count))
            resolution = resolve.Resolution("broad", name)

    # A name that no single protein owns may be a small molecule: cholesterol,
    # ATP, caffeine.  Only a molecule really called that takes over.
    decisive = (resolution.kind == "protein" and resolution.top_tier >= resolve.TIER_NAME
                and resolution.leader_count > 0)
    if query.kind == "auto" and query.mode != "text" and not decisive \
            and not query.has_filters():
        molecule = chemcomp.lookup(name, allow_close=False, allow_id=False)
        unknown_protein = resolution.kind in ("none", "broad")
        if molecule is not None and not molecule.strong and (
                not unknown_protein or len(name) <= 4):
            # only a lesser-known synonym of some molecule ("ras"), while
            # proteins of that name exist: the proteins keep the name
            molecule = None
        if molecule is None and resolution.kind == "none":
            # no protein and no PDB molecule is called this: PubChem as a last
            # source, before concluding that nothing exists (or that it is a typo)
            try:
                other = pubchem.lookup(name)
            except http.NetworkError:
                other = None           # PubChem busy or unreachable: carry on
            if other is not None and other.has_shape():
                molecule = other
        if molecule is not None:
            return _molecule_result(query, name, spelling_notes, molecule, resolution)

    # A word that names a whole family ("kinase", "GPCR", "globin"): list its
    # best-studied members and let the user choose, instead of loading one.
    family_word = " ".join(re.findall(r"[a-z0-9]+", name.lower())) in resolve.FAMILY_TERMS
    if query.kind == "auto" and query.mode != "text" and not query.has_filters() and (
            family_word or resolution.kind == "broad" or first.kind == "broad" or (
                first.kind in ("protein", "ambiguous")
                and first.top_tier <= resolve.TIER_PREFIX)):
        # a name proteins are called by ("hemoglobin subunit alpha"): only a
        # human family counts; a bare word ("globin") may look further
        everywhere = family_word or not (first.kind == "protein"
                                         and first.top_tier >= resolve.TIER_PREFIX)
        members = resolve.family(name, query.organism, everywhere)
        if members:
            family = resolve.Resolution(
                "ambiguous", name, choices=members[:resolve.FAMILY_SIZE],
                hidden=len(members) - resolve.FAMILY_SIZE)
            family.top_tier = resolve.TIER_WORDS
            family.leader_count = members[0].count
            family.family = True
            scope = (" These are human members." if members.human else
                     "" if query.organism is not None else
                     " Few human members have structures, so all organisms are included.")
            result = Result("ambiguous", query)
            result.resolution = family
            result.notes = early_notes + [
                "Nothing was loaded: '%s' names a whole family of proteins.%s"
                % (name, scope)]
            return result

    if resolution.kind == "ambiguous":
        result = Result("ambiguous", query)
        result.resolution = resolution
        result.notes = notes
        return result

    if resolution.kind == "protein":
        return _protein_search(query, name, resolution, notes)

    if resolution.kind == "none":
        suggestion = fuzzy.suggest(name)
        if suggestion:
            result = Result("suggest", query)
            result.suggestion = suggestion
            result.notes = notes
            result.message = "No structure found for '%s' - did you mean '%s'?" % (
                query.text, suggestion)
            return result

    return _text_search(query, name, notes)


def _family_members(name, groups):
    """How many choices are named after the word itself plus a short tag:
    'Aquaporin-1', 'Aquaporin Z', 'Myosin-7' (not 'Ribosome biogenesis protein')."""
    pattern = re.compile(r"^%s(?:[\s-]?[a-z0-9]{1,3})?$" % re.escape(name.strip().lower()))
    count = 0
    for group in groups:
        names = {p.name.lower() for p in group.members if p.name}
        if any(pattern.match(n) for n in names):
            count += 1
    return count


def _molecule_result(query, name, notes, molecule, resolution=None):
    if molecule is None:
        result = Result("none", query)
        result.notes = notes
        result.message = ("No small molecule called '%s' was found, neither in the PDB's "
                          "dictionary of molecules nor in PubChem. Check the spelling or "
                          "try a more specific name." % name)
        return result
    if molecule.source == "PubChem" and not molecule.has_shape():
        result = Result("none", query)
        result.notes = notes
        result.message = ("PubChem knows the small molecule '%s' (%s, %s g/mol), but has "
                          "no 3D shape for it, so there is nothing to show. That is usual "
                          "for very large or floppy molecules and for metal complexes."
                          % (molecule.title(), molecule.formula or "formula unknown",
                             "%.0f" % molecule.weight if molecule.weight else "?"))
        return result
    result = Result("molecule", query)
    result.molecule = molecule
    result.resolution = resolution
    result.notes = notes
    if molecule.source == "PDB":
        molecule.entries = chemcomp.count_entries(molecule.comp_id)
    return result


def _protein_search(query, name, resolution, notes):
    group = resolution.group
    lead = group.lead
    targets = dict((p.accession, p) for p in group.members if p.pdb)
    target_ids = group.targets()

    total, ids = (0, [])
    if target_ids:
        total, ids = rcsb.search_accessions(target_ids, query)

    if not ids:
        if target_ids and query.has_filters():
            result = Result("nofilter", query)
            result.resolution = resolution
            result.notes = notes
            result.message = (
                "%s has %d experimental structure(s), but none pass your filters (%s). "
                "Try relaxing one of them." % (group.label(), group.count,
                                               _describe_filters(query)))
            return result
        return _alphafold_fallback(query, resolution, notes)

    # Cheap first pass using what UniProt already told us (method, resolution,
    # residue range), so only the most promising entries need full metadata.
    def preliminary(pdb_id):
        best = None
        for accession in target_ids:
            protein = targets[accession]
            ref = protein.pdb.get(pdb_id)
            if ref is None:
                continue
            coverage = min(1.0, ref.covered() / float(protein.length)) if protein.length else 0.5
            value = (rank.WEIGHTS["coverage"] * coverage
                     + rank.WEIGHTS["resolution"] * rank.resolution_score(ref.resolution)
                     + rank.WEIGHTS["method"] * rank.method_score([ref.method])
                     + rank.WEIGHTS["organism"] * (1.0 if query.organism else
                                                   organisms.preference(protein.taxid)))
            best = value if best is None else max(best, value)
        if best is None:      # released after UniProt's last update: middle of the road
            best = 0.6 * (rank.WEIGHTS["coverage"] + rank.WEIGHTS["resolution"]
                          + rank.WEIGHTS["method"] + rank.WEIGHTS["organism"])
        if pdb_id in rank.classics():
            best += rank.CLASSIC_BONUS + rank.WEIGHTS["recency"]
        return best

    pool_size = max(DETAIL_POOL, 4 * query.limit)
    shortlist = sorted(ids, key=lambda i: (-preliminary(i), i))[:pool_size]
    details = rcsb.entries(shortlist)
    extra = []
    for protein in group.members:
        extra.extend([protein.gene, protein.name])
    ctx = rank.Context(name, targets=targets, organism_filter=query.organism is not None,
                       extra_words=extra)
    hits = rank.rank(details.values(), ctx)

    result = Result("results", query)
    result.mode = "protein"
    result.resolution = resolution
    result.total = total
    if query.pick > len(hits):
        failed = Result("error", query)
        failed.message = ("pick=%d asks for row %d, but only %d structure%s matched."
                          % (query.pick, query.pick, len(hits),
                             "" if len(hits) == 1 else "s"))
        return failed
    result.hits = hits[:query.limit]
    result.protein = lead
    result.notes = notes
    return result


def _text_search(query, name, notes):
    total, scored = rcsb.search_text(name, query, rows=max(60, 4 * query.limit),
                                     titles_first=True)
    if not scored:
        result = Result("none", query)
        result.notes = notes
        extra = ""
        if query.has_filters():
            extra = " with your filters (%s)" % _describe_filters(query)
        result.message = "No structure found for '%s'%s." % (query.text, extra)
        return result
    # Text relevance from the PDB is coarse, so classic structures whose
    # description mentions every query word are added to the candidates
    # (they must still pass the filters).
    words = [w for w in name.lower().split() if w not in resolve.STOPWORDS]
    wanted = [pdb_id for pdb_id, _ in scored]
    classic_ids = [pdb_id for pdb_id, why in rank.classics().items()
                   if words and all(w in why.lower() for w in words) and pdb_id not in wanted]
    details = rcsb.entries(wanted + classic_ids)
    for pdb_id in classic_ids:
        if pdb_id in details and not _passes(details[pdb_id], query):
            del details[pdb_id]
    accessions = set()
    for entry in details.values():
        accessions.update(entry.coverage)
    lengths = uniprot.lengths(accessions)
    ctx = rank.Context(name, lengths=lengths, text_scores=dict(scored),
                       organism_filter=query.organism is not None)
    hits = rank.rank(details.values(), ctx)
    # a textbook structure of exactly this thing goes first: that is what a
    # student searching for "ubiquitin" or "ribosome" is looking for
    for hit in hits:
        hit.textbook = hit.pdb_id in classic_ids or (
            bool(words) and bool(hit.classic)
            and all(w in hit.classic.lower() for w in words))
    hits.sort(key=lambda h: (not h.textbook, -h.score))
    result = Result("results", query)
    result.mode = "text"
    result.total = total
    result.hits = hits[:query.limit]
    result.notes = notes
    return result


def _alphafold_fallback(query, resolution, notes):
    group = resolution.group
    lead = group.lead
    model = alphafold.prediction(lead.accession)
    result = Result("alphafold" if model else "none", query)
    result.resolution = resolution
    result.protein = lead
    result.model = model
    result.notes = notes
    if model is None:
        result.message = ("%s (%s, %s) has no experimental structure in the PDB, and "
                          "AlphaFold DB has no model for it either."
                          % (lead.label(), lead.accession, lead.organism))
    return result


def _passes(entry, query):
    """Check an entry against the filters locally (same rules as the RCSB search)."""
    if query.organism is not None and query.organism.taxid and \
            query.organism.taxid not in [tax for _, tax in entry.organisms]:
        return False
    if query.method and filters.METHODS[query.method] not in [m.upper() for m in entry.methods]:
        return False
    if query.res_max is not None and (entry.resolution is None
                                      or entry.resolution > query.res_max):
        return False
    if query.res_min is not None and (entry.resolution is None
                                      or entry.resolution < query.res_min):
        return False
    if any(ligand not in entry.ligands for ligand in query.ligands):
        return False
    if query.year_min is not None and (entry.year or 0) < query.year_min:
        return False
    if query.year_max is not None and (entry.year or 9999) > query.year_max:
        return False
    return True


def _describe_filters(query):
    words = query.filter_words()
    if query.organism is not None and query.organism_from_text:
        words.insert(0, "organism=%s" % query.organism.typed)
    return " ".join(words)
