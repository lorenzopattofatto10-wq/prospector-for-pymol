"""Name resolution: free text -> UniProt protein(s), or an honest "I'm not sure".

The outcome is one of
    protein    one protein (or one complex / set of orthologs) was identified
    ambiguous  several clearly different proteins fit - the user must choose
    broad      the word names a whole class (e.g. "kinase") - use text search
    none       UniProt knows nothing by that name
"""
from __future__ import annotations

import re

from . import organisms, uniprot

STOPWORDS = {"protein", "proteins", "the", "of", "a", "an", "and", "structure", "enzyme",
             "domain", "complex", "from"}

# How well a UniProt entry matches the typed name.
TIER_GENE = 4      # equals the gene symbol or a gene synonym (p53 -> TP53)
TIER_NAME = 3      # equals a protein name (insulin -> "Insulin")
TIER_PREFIX = 2    # a protein name starts with it (hemoglobin -> "Hemoglobin subunit alpha")
TIER_WORDS = 1     # all words appear in a protein name
TIER_NONE = 0

SUBSTANTIAL = 0.10   # a rival needs >= 10% of the leader's structures to count
BROAD_TOTAL = 100    # more UniProt matches than this (and no exact match) = a class word
MAX_CHOICES = 6
MAX_TARGETS = 10


def _norm(text):
    return re.sub(r"[^a-z0-9]", "", text.lower())


def _words(text):
    return re.findall(r"[a-z0-9]+", text.lower())


def match_tier(protein, name):
    target = _norm(name)
    query_words = _words(name)
    for gene in [protein.gene] + protein.gene_synonyms:
        if gene and _norm(gene) == target:
            return TIER_GENE
    best = TIER_NONE
    for full in protein.all_names():
        words = _words(full)
        if words[:1] == ["protein"] and len(words) > 1:
            words = words[1:]            # "Protein myomaker" == "myomaker"
        if "".join(words) == target:
            return TIER_NAME
        if words[:len(query_words)] == query_words:
            best = max(best, TIER_PREFIX)
        elif all(w in words for w in query_words):
            best = max(best, TIER_WORDS)
    return best


class Group(object):
    """Proteins that count as 'the same answer': orthologs with the same name,
    or partners that are nearly always in the same PDB entries (a complex)."""

    def __init__(self, members, tiers):
        def order(protein):
            return (-organisms.preference(protein.taxid), -len(protein.pdb),
                    -protein.annotation_score, protein.accession)
        self.members = sorted(members, key=order)
        self.tier = max(tiers[p.accession] for p in members)
        self.pdb_ids = set()
        for protein in members:
            self.pdb_ids.update(protein.pdb)

    @property
    def count(self):
        return len(self.pdb_ids)

    @property
    def lead(self):
        """The member that represents the group (most structures, human first)."""
        with_structures = [p for p in self.members if p.pdb]
        pool = with_structures or self.members
        return max(pool, key=lambda p: (len(p.pdb), organisms.preference(p.taxid)))

    def distinct_names(self):
        seen, names = set(), []
        for protein in sorted(self.members, key=lambda p: -len(p.pdb)):
            if protein.name.lower() not in seen:
                seen.add(protein.name.lower())
                names.append(protein)
        return names

    def label(self, max_names=2):
        named = self.distinct_names()
        text = " + ".join(p.label() for p in named[:max_names])
        if len(named) > max_names:
            text += " + %d related" % (len(named) - max_names)
        return text

    def organisms_text(self, max_items=3):
        seen = []
        for protein in self.members:
            if not protein.pdb and self.count:
                continue
            name = protein.organism_common or protein.organism
            name = name.split("(")[0].strip()
            if name and name not in seen:
                seen.append(name)
        text = ", ".join(seen[:max_items])
        if len(seen) > max_items:
            text += " +%d more" % (len(seen) - max_items)
        return text

    def targets(self):
        """Accessions to search the PDB with."""
        return [p.accession for p in self.members if p.pdb][:MAX_TARGETS]

    def suggested_query(self, others=()):
        """What to type to get exactly this group: the gene symbol when that is
        unambiguous among the choices, otherwise the UniProt accession."""
        lead = self.lead
        rivals = [g.lead.gene.lower() for g in others if g is not self]
        if len(lead.gene) >= 3 and lead.gene.lower() not in rivals:
            return lead.gene
        return lead.accession


class Resolution(object):
    def __init__(self, kind, name, group=None, choices=None, hidden=0, notes=None):
        self.kind = kind            # protein | ambiguous | broad | none
        self.name = name
        self.group = group
        self.choices = choices or []
        self.hidden = hidden        # number of smaller matches not listed
        self.top_tier = TIER_GENE   # how well the best candidate matched the name
        self.leader_count = 0       # structures of the best-populated candidate
        self.notes = notes or []
        self.family = False         # the choices are members of a protein family


def _group(candidates, tiers):
    """Union-find over the candidates."""
    parent = dict((p.accession, p.accession) for p in candidates)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        parent[find(a)] = find(b)

    for i, first in enumerate(candidates):
        for second in candidates[i + 1:]:
            same_name = _norm(first.name) == _norm(second.name)
            # orthologs usually share the gene symbol (cas9, TP53/Tp53, recA)
            same_gene = (len(first.gene) >= 3
                         and first.gene.lower() == second.gene.lower())
            # Partners of an obligate complex (hemoglobin alpha + beta) share nearly
            # all their PDB entries.  Insulin and its receptor share only some.
            shared = len(set(first.pdb) & set(second.pdb))
            either = len(set(first.pdb) | set(second.pdb))
            complex_partners = either > 0 and shared >= 0.5 * either
            if same_name or same_gene or complex_partners:
                union(first.accession, second.accession)

    buckets = {}
    for protein in candidates:
        buckets.setdefault(find(protein.accession), []).append(protein)
    return [Group(members, tiers) for members in buckets.values()]


CLASS_WORDS = {"kinase", "protein", "enzyme", "protease", "phosphatase", "polymerase"}


def resolve(name, organism=None):
    """Resolve a protein name.  If the full phrase leads nowhere useful (nothing
    found, or a protein with no experimental structure) and the phrase is a
    name plus a class word - "abl kinase" - the name alone is tried as well."""
    first = _resolve(name, organism)
    words = name.split()
    dead_end = first.kind == "none" or (first.kind == "protein" and first.group.count == 0)
    if dead_end and len(words) > 1:
        core = [w for w in words if w.lower() not in CLASS_WORDS]
        if len(core) == 1 and len(core) < len(words):
            second = _resolve(core[0], organism)
            if second.kind == "protein" and second.group.count > 0:
                second.name = name
                second.notes.append("No structure is filed under '%s' as a phrase, so I "
                                    "used '%s' on its own." % (name, core[0]))
                return second
    return first


def _resolve(name, organism=None):
    """Resolve a protein name.  ``organism`` is an organisms.Organism or None."""
    name = name.strip()
    notes = []

    if uniprot.ACCESSION_RE.match(name):
        protein = uniprot.by_accession(name)
        if protein:
            return Resolution("protein", name, Group([protein], {protein.accession: TIER_GENE}),
                              notes=["Recognised %s as a UniProt accession." % name.upper()])

    found = []
    total = 0
    if organism is not None:
        hits = uniprot.search_name(name, organism)
        total = hits.total
        found.extend(hits)
    else:
        human = uniprot.search_name(name, organisms.Organism("Homo sapiens", organisms.HUMAN))
        total = human.total
        found.extend(human)
        found.extend(uniprot.search_name(name, extra=" AND (structure_3d:true)"))
        if not found:
            found.extend(uniprot.search_name(name))

    words = [w for w in name.split() if w.lower() not in STOPWORDS]
    if not found and words and (len(words) > 1 or len(words) != len(name.split())):
        hits = uniprot.search_words(words, organism)
        total = max(total, hits.total)
        found.extend(hits)

    # de-duplicate, score how well each candidate matches the typed name
    candidates, tiers = [], {}
    match_name = " ".join(words) if words else name
    for protein in found:
        if protein.accession in tiers:
            continue
        tier = max(match_tier(protein, name), match_tier(protein, match_name))
        if tier == TIER_NONE:
            continue
        tiers[protein.accession] = tier
        candidates.append(protein)
    if not candidates:
        return Resolution("none", name, notes=notes)

    top_tier = max(tiers.values())
    if (top_tier <= TIER_PREFIX and total > BROAD_TOTAL) or \
            (top_tier <= TIER_WORDS and total > 25):
        return Resolution("broad", name, notes=[
            "'%s' is a whole family of proteins (%d in UniProt), not one protein, so I "
            "searched the text of every PDB entry instead." % (name, total)])

    if top_tier == TIER_GENE:
        floor = TIER_GENE          # a gene symbol is an identifier: decisive
    elif top_tier == TIER_NAME:
        floor = TIER_PREFIX        # "insulin" vs "insulin receptor": keep both in play
    else:
        floor = top_tier
    pool = [p for p in candidates if tiers[p.accession] >= floor]
    groups = _group(pool, tiers)
    groups.sort(key=lambda g: (-g.count, -g.tier, g.lead.accession))

    leader_count = groups[0].count
    if leader_count == 0:
        # nobody has an experimental structure: take the best name match
        best = max(groups, key=lambda g: (g.tier, organisms.preference(g.lead.taxid),
                                          g.lead.annotation_score))
        return Resolution("protein", name, best, notes=notes)

    def finish(resolution):
        resolution.top_tier = top_tier
        resolution.leader_count = leader_count
        return resolution

    exact = [g for g in groups if g.tier == top_tier and top_tier >= TIER_NAME and g.count]
    best_exact = exact[0] if exact else None   # groups are sorted by size

    exact_leads = top_tier == TIER_NAME and groups[0].tier == TIER_NAME

    def substantial(group):
        share = SUBSTANTIAL
        if exact_leads and group.tier < TIER_NAME:
            # you typed a protein's exact name and it is also the best-studied
            # match: a protein that merely starts with those words ("insulin
            # receptor substrate") only counts as a rival if it is a big one
            share = 0.25
        if group.count >= max(3, share * leader_count):
            return True
        # an exact name match stays in play even if a looser match has more structures
        return group is best_exact and groups[0].tier < top_tier

    strong = [g for g in groups if substantial(g)] or [groups[0]]
    if len(strong) == 1:
        return finish(Resolution("protein", name, strong[0], notes=notes))

    strong.sort(key=lambda g: (-g.tier, -g.count))
    return finish(Resolution("ambiguous", name, choices=strong[:MAX_CHOICES],
                             hidden=len(groups) - min(len(strong), MAX_CHOICES),
                             notes=notes))


# ------------------------------------------------------------------ families
# A word such as "kinase" or "GPCR" names a whole family.  `family` lists its
# best-studied members - human first - so the user can choose one.
COMPLEXES = {"ribosome", "ribosomes", "proteasome", "proteasomes", "nucleosome",
             "nucleosomes", "spliceosome", "spliceosomes", "photosystem", "photosystems",
             "photosystem i", "photosystem ii", "exosome", "apoptosome", "inflammasome",
             "replisome", "chromatin"}
# one protein made from several identical genes (UBB, UBC, RPS27A, UBA52 all give
# ubiquitin): UniProt files them as a family, but there is nothing to choose
ONE_PROTEIN = {"ubiquitin", "polyubiquitin", "calmodulin"}
# parts found in thousands of unrelated proteins: no family to choose from
DOMAINS = {"zinc finger", "zinc fingers", "leucine zipper", "homeodomain", "homeobox",
           "bromodomain", "sh2 domain", "sh3 domain", "pdz domain", "helix turn helix"}
# other ways UniProt says the same thing (names, keywords, family sentences)
FAMILY_TERMS = {
    "gpcr": ["G protein-coupled receptor"],
    "gpcrs": ["G protein-coupled receptor"],
    "g protein coupled receptor": ["G protein-coupled receptor"],
    "g protein coupled receptors": ["G protein-coupled receptor"],
    "nuclear receptor": ["nuclear receptor", "nuclear hormone receptor"],
    "nuclear receptors": ["nuclear receptor", "nuclear hormone receptor"],
    "tyrosine kinase": ["tyrosine kinase", "tyrosine-protein kinase", "Tyr protein kinase"],
    "tyrosine kinases": ["tyrosine kinase", "tyrosine-protein kinase", "Tyr protein kinase"],
    "serine protease": ["serine protease", "peptidase S1"],
    "ubiquitin ligase": ["ubiquitin ligase", "E3 ubiquitin-protein ligase",
                         "ubiquitin-protein ligase"],
    "e3 ligase": ["E3 ubiquitin-protein ligase"],
    "heat shock protein": ["heat shock protein", "heat shock"],
}
FAMILY_SIZE = 6       # members listed
FAMILY_FLOOR = 10     # the best-studied member of a real family has at least this many
_GREEK = set("alpha beta gamma delta epsilon zeta eta theta iota kappa lambda mu nu xi "
             "omicron pi rho sigma tau upsilon phi chi psi omega".split())
_TAGS = set("type subunit chain member isoform class form homolog kda like isozyme variant "
            "family superfamily subfamily".split())


def _flat(text):
    """'G-protein coupled receptor 1 (TC 1.A.8)' -> 'g protein coupled receptor 1'."""
    text = re.sub(r"\([^()]*\)|\[[^\[\]]*\]", " ", text.lower())
    return " ".join(re.findall(r"[a-z0-9.,']+", text.replace("-", " ").replace("/", " ")))


def _heads(text, phrase):
    """Is ``phrase`` what ``text`` is: 'Caspase-3', 'Tubulin alpha-1B chain',
    'Actin, cytoplasmic 1' for caspase / tubulin / actin - but not 'Cyclin-dependent
    kinase 2' for cyclin, nor 'Histone deacetylase 1' for histone."""
    flat, words = _flat(text), _flat(phrase)
    for match in re.finditer(r"(?<![a-z0-9])%s(?![a-z0-9])" % re.escape(words), flat):
        rest = flat[match.end():]
        if rest.strip() in ("", ".") or rest.startswith(", "):
            return True
        # "4,6-dehydratase": a comma inside a chemical name does not end the name
        tags = re.findall(r"[a-z0-9.,']+", re.split(r",\s", rest)[0])
        if all(re.search(r"\d", t) or len(t) <= 3 or t in _GREEK or t in _TAGS
               for t in tags):
            return True
    return False


def family_member(protein, phrases):
    """How ``protein`` belongs to the family: "name" (the word heads its name),
    "family" (UniProt files it in that family), "keyword" (a UniProt keyword
    says so) - or None."""
    if any(_heads(protein.name, p) for p in phrases):
        return "name"
    if any(_heads(sentence.replace("Belongs to the ", ""), p)
           for sentence in getattr(protein, "families", []) for p in phrases):
        return "family"
    keys = set(_norm(k) for k in getattr(protein, "keywords", []))
    if any(_norm(p) in keys for p in phrases):
        return "keyword"
    return None


def family_phrases(name):
    key = " ".join(re.findall(r"[a-z0-9]+", name.lower()))
    return FAMILY_TERMS.get(key, [name])


def family(name, organism=None, everywhere=True):
    """The best-studied members of the family ``name`` names, as Groups (most
    structures first), or [] when ``name`` is not a family: one protein, a
    complex, or a word with fewer than two members of any weight.  Human
    proteins are tried first; ``everywhere``: then all organisms ("globin":
    few human globins have structures, sperm whale myoglobin has hundreds)."""
    key = " ".join(re.findall(r"[a-z0-9]+", name.lower()))
    if key in COMPLEXES or key in DOMAINS or key in ONE_PROTEIN or not key:
        return []
    phrases = family_phrases(name)
    places = [organism] if organism is not None else \
        [organisms.Organism("Homo sapiens", organisms.HUMAN)] + ([None] if everywhere else [])
    for place in places:
        found = uniprot.search_family(phrases, place)
        how = {}
        for protein in found:
            if protein.pdb and protein.accession not in how:
                kind = family_member(protein, phrases)
                if kind:
                    how[protein.accession] = (kind, protein)

        def groups_of(kinds):
            members = [p for k, p in how.values() if k in kinds]
            groups = _group(members, dict((p.accession, TIER_WORDS) for p in members))
            groups.sort(key=lambda g: (-g.count, g.lead.accession))
            return groups

        def strength(groups):
            return sum(g.count for g in groups[:FAMILY_SIZE])

        # a name or a family sentence says exactly what a protein is; a keyword is
        # looser ("Protease" is also given to precursors of protein-cutting enzymes),
        # so keyword members count only when they make a far better list
        # (protease, GPCR, ion channel, lectin: words UniProt rarely puts in names)
        precise = groups_of(("name", "family"))
        loose = groups_of(("name", "family", "keyword"))
        groups = loose if strength(loose) > 2 * strength(precise) else precise
        if not groups:
            continue
        top = groups[0].count
        weighty = [g for g in groups if g.count >= max(3, 0.1 * top)]
        if len(weighty) >= 2 and top >= FAMILY_FLOOR:
            listed = FamilyList(groups)
            listed.human = place is not None and place.taxid == organisms.HUMAN
            return listed
    return []


class FamilyList(list):
    """Groups, plus whether only human proteins were looked at."""
    human = False
