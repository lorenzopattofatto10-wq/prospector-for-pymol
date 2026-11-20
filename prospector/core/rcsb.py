"""RCSB PDB clients: Search API v2 and Data API (GraphQL).

Checked against the live schemas:
  POST https://search.rcsb.org/rcsbsearch/v2/query     (HTTP 204 = no hits)
  POST https://data.rcsb.org/graphql
"""
from __future__ import annotations

import re

from . import cache, http
from .filters import METHODS

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
DOWNLOAD_URL = "https://files.rcsb.org/download/%s.cif"

# Searchable attributes (names verified in the Search API metadata schema).
A_ACCESSION = ("rcsb_polymer_entity_container_identifiers."
               "reference_sequence_identifiers.database_accession")
A_METHOD = "exptl.method"
A_RESOLUTION = "rcsb_entry_info.resolution_combined"
A_RELEASE = "rcsb_accession_info.initial_release_date"
A_LIGAND = "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id"
A_TAX_ID = "rcsb_entity_source_organism.taxonomy_lineage.id"
A_TAX_NAME = "rcsb_entity_source_organism.taxonomy_lineage.name"

# Buffer components, cryo-protectants and common ions: listed after real ligands.
ADDITIVES = {"HOH", "DOD", "EDO", "GOL", "PEG", "PGE", "PG4", "1PE", "P6G", "MPD", "DMS",
             "ACT", "ACY", "FMT", "SO4", "PO4", "CL", "NA", "K", "BR", "IOD", "NO3", "TRS",
             "EPE", "MES", "BME", "DTT", "IPA", "EOH", "CIT", "TLA", "IMD", "NH4", "SCN",
             "BOG", "LDA", "UNX", "UNL"}

# Detergents (and lipidic-cubic-phase lipids) used to pull membrane proteins out of
# the membrane and crystallise them: part of the experiment, not of the protein.
DETERGENTS = {
    "BNG": "beta-nonylglucoside", "BOG": "beta-octylglucoside", "LDA": "LDAO",
    "LMT": "dodecyl maltoside (DDM)", "DMU": "decyl maltoside", "UMQ": "undecyl maltoside",
    "HTG": "heptyl thioglucoside", "C8E": "C8E4", "SDS": "sodium dodecyl sulfate",
    "CPS": "CHAPS", "OLC": "monoolein (from crystallisation in lipid)",
    "OLA": "oleic acid", "Y01": "cholesterol hemisuccinate", "LMN": "LMNG",
}
ADDITIVES |= set(DETERGENTS)

# Simple ions: sometimes part of how the protein works (K+ in a potassium channel),
# sometimes just salt from the solution - never called "not real".
IONS = {"NA", "K", "CL", "BR", "IOD", "NH4", "LI", "RB", "CS", "F"}

# Membrane lipids and sterols that come with membrane proteins (from the cell
# membrane, a nanodisc or the lipid used for crystallisation).
LIPIDS = {
    "CLR": "cholesterol", "PLM": "palmitic acid", "STE": "stearic acid",
    "MYR": "myristic acid", "OLA": "oleic acid", "POV": "POPC", "PCW": "DOPC",
    "PEE": "DOPE", "PC1": "phosphatidylcholine", "PTY": "phosphatidylethanolamine",
    "3PE": "phosphatidylethanolamine", "CDL": "cardiolipin", "LHG": "phosphatidylglycerol",
    "PGW": "phosphatidylglycerol", "PEF": "phosphatidylethanolamine", "LPP": "phosphatidic acid",
    "D10": "decane", "HEX": "hexane", "OCT": "octane", "R16": "hexadecane",
    "LMG": "galactolipid", "DGD": "galactolipid", "SQD": "sulfolipid",
    "PCF": "phosphatidylcholine", "P5S": "phosphatidylserine", "PIO": "phosphatidylinositol",
}
# Sugars that form the sugar chains attached to many proteins (glycosylation).
GLYCANS = {"NAG", "NDG", "MAN", "BMA", "FUC", "FUL", "GAL", "GLA", "SIA", "BGC", "GLC",
           "XYS", "A2G", "NGA"}

PDB_ID_RE = re.compile(r"^[1-9][A-Za-z0-9]{3}$")


def looks_like_pdb_id(text):
    return bool(PDB_ID_RE.match(text)) and any(c.isalpha() for c in text)


def _text(attribute, operator, value):
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attribute, "operator": operator, "value": value}}


def filter_nodes(query):
    """Translate a filters.Query into Search API nodes."""
    nodes = []
    if query.organism is not None:
        if query.organism.taxid:
            nodes.append(_text(A_TAX_ID, "exact_match", str(query.organism.taxid)))
        else:
            nodes.append(_text(A_TAX_NAME, "exact_match", query.organism.name))
    if query.method:
        nodes.append(_text(A_METHOD, "exact_match", METHODS[query.method]))
    if query.res_max is not None:
        nodes.append(_text(A_RESOLUTION, "less_or_equal", query.res_max))
    if query.res_min is not None:
        nodes.append(_text(A_RESOLUTION, "greater_or_equal", query.res_min))
    for ligand in query.ligands:
        nodes.append(_text(A_LIGAND, "exact_match", ligand))
    if query.year_min is not None:
        nodes.append(_text(A_RELEASE, "greater_or_equal", "%d-01-01" % query.year_min))
    if query.year_max is not None:
        nodes.append(_text(A_RELEASE, "less_or_equal", "%d-12-31" % query.year_max))
    return nodes


def _search(nodes, options):
    node = nodes[0] if len(nodes) == 1 else {
        "type": "group", "logical_operator": "and", "nodes": nodes}
    payload = {"query": node, "return_type": "entry", "request_options": options}
    data = http.request_json(SEARCH_URL, payload, ttl=cache.DAY)
    if not data:   # HTTP 204: nothing found
        return 0, []
    hits = []
    for item in data.get("result_set", []):
        if isinstance(item, dict):
            hits.append((item.get("identifier", "").upper(), item.get("score", 0.0)))
        else:
            hits.append((str(item).upper(), 1.0))
    return data.get("total_count", len(hits)), hits


def search_accessions(accessions, query):
    """All entries that contain any of the UniProt accessions and pass the filters.
    Returns (total, [pdb ids])."""
    nodes = [_text(A_ACCESSION, "in", list(accessions))] + filter_nodes(query)
    total, hits = _search(nodes, {"return_all_hits": True, "results_verbosity": "compact"})
    return total, [pdb_id for pdb_id, _ in hits]


def search_text(text, query, rows=60, titles_first=False):
    """Text search plus filters.  Returns (total, [(pdb id, relevance 0..1)]).

    With ``titles_first`` the entry TITLES are searched first: an entry whose
    title names the thing is far more likely to be about it than one that
    merely mentions the word somewhere.  If fewer than five titles match, the
    full text of every entry is searched instead."""
    options = {"paginate": {"start": 0, "rows": rows},
               "sort": [{"sort_by": "score", "direction": "desc"}]}
    if titles_first:
        nodes = [_text("struct.title", "contains_words", text)] + filter_nodes(query)
        total, hits = _search(nodes, options)
        if total >= 5:
            return total, hits
    nodes = [{"type": "terminal", "service": "full_text", "parameters": {"value": text}}]
    nodes += filter_nodes(query)
    return _search(nodes, options)


_ENTRY_QUERY = """query($ids:[String!]!){entries(entry_ids:$ids){
 rcsb_id
 struct{title}
 exptl{method}
 rcsb_entry_info{resolution_combined deposited_polymer_entity_instance_count
   deposited_modeled_polymer_monomer_count deposited_polymer_monomer_count}
 rcsb_accession_info{initial_release_date}
 polymer_entities{
   rcsb_polymer_entity{pdbx_description}
   entity_poly{rcsb_mutation_count rcsb_entity_polymer_type}
   rcsb_entity_source_organism{ncbi_scientific_name ncbi_taxonomy_id}
   rcsb_polymer_entity_container_identifiers{auth_asym_ids uniprot_ids}
   rcsb_polymer_entity_align{reference_database_name reference_database_accession
     aligned_regions{ref_beg_seq_id length}}
 }
 nonpolymer_entities{rcsb_nonpolymer_entity_container_identifiers{nonpolymer_comp_id}}
}}"""
_ENTRY_QUERY = " ".join(_ENTRY_QUERY.split())


_DETAIL_QUERY = """query($id:String!){entry(entry_id:$id){
 rcsb_id struct{title} exptl{method}
 refine{ls_R_factor_R_free ls_R_factor_R_work}
 rcsb_entry_info{resolution_combined}
 rcsb_accession_info{initial_release_date deposit_date}
 rcsb_primary_citation{title rcsb_journal_abbrev year pdbx_database_id_DOI
   pdbx_database_id_PubMed rcsb_authors journal_volume page_first}
 polymer_entities{
   rcsb_polymer_entity{pdbx_description}
   entity_poly{rcsb_sample_sequence_length rcsb_entity_polymer_type rcsb_mutation_count}
   rcsb_entity_source_organism{ncbi_scientific_name}
   rcsb_polymer_entity_container_identifiers{entity_id auth_asym_ids uniprot_ids}
   rcsb_polymer_entity_align{reference_database_name reference_database_accession
     aligned_regions{entity_beg_seq_id ref_beg_seq_id length}}
   polymer_entity_instances{rcsb_polymer_entity_instance_container_identifiers{
     auth_asym_id asym_id auth_to_entity_poly_seq_mapping}}
 }
 nonpolymer_entities{nonpolymer_comp{chem_comp{id name formula}}}
}}"""
_DETAIL_QUERY = " ".join(_DETAIL_QUERY.split())


def entry_detail(pdb_id):
    """Everything about one entry (raw GraphQL dict), or None if it does not exist.
    Used by `about`, `features` and `pocket`."""
    data = http.request_json(GRAPHQL_URL, {"query": _DETAIL_QUERY,
                                           "variables": {"id": pdb_id.upper()}},
                             ttl=7 * cache.DAY, timeout=20) or {}
    return (data.get("data") or {}).get("entry")


def ligand_names(detail):
    """{ligand code: (chemical name, formula)} from an entry_detail() answer."""
    names = {}
    for item in (detail or {}).get("nonpolymer_entities") or []:
        comp = ((item.get("nonpolymer_comp") or {}).get("chem_comp") or {})
        if comp.get("id"):
            names[comp["id"]] = (comp.get("name") or "", comp.get("formula") or "")
    return names


class Entry(object):
    """Metadata of one PDB entry, flattened from the GraphQL answer."""

    def __init__(self, record):
        info = record.get("rcsb_entry_info") or {}
        self.pdb_id = (record.get("rcsb_id") or "").upper()
        self.title = ((record.get("struct") or {}).get("title") or "").strip()
        self.methods = [m.get("method", "") for m in (record.get("exptl") or [])]
        resolutions = info.get("resolution_combined") or []
        self.resolution = min(resolutions) if resolutions else None
        ligands = []
        for item in record.get("nonpolymer_entities") or []:
            code = ((item.get("rcsb_nonpolymer_entity_container_identifiers") or {})
                    .get("nonpolymer_comp_id"))
            if code and code not in ligands:
                ligands.append(code)
        # interesting ligands first, crystallisation additives and ions last
        self.ligands = sorted(ligands, key=lambda c: (c in ADDITIVES, ligands.index(c)))
        self.chain_count = info.get("deposited_polymer_entity_instance_count") or 0
        # how much of the deposited sequence was actually built into the model
        deposited = info.get("deposited_polymer_monomer_count") or 0
        modelled = info.get("deposited_modeled_polymer_monomer_count") or 0
        self.modelled_fraction = min(1.0, modelled / float(deposited)) if deposited else 1.0
        self.protein_entities = []   # one list of UniProt accessions per protein entity
        date = (record.get("rcsb_accession_info") or {}).get("initial_release_date") or ""
        self.year = int(date[:4]) if date[:4].isdigit() else None

        self.organisms = []     # [(scientific name, taxid)]
        self.descriptions = []  # one per polymer entity
        self.chains = []        # author chain ids
        self.coverage = {}      # {uniprot accession: residues covered}
        regions = {}            # {accession: [(first residue, last residue)]}
        largest = {}            # {accession: size of its biggest chain so far}
        self.mutations = {}     # {uniprot accession: engineered mutations in that chain}
        self.organism_of = {}   # {uniprot accession: scientific name of that chain}
        for entity in record.get("polymer_entities") or []:
            description = ((entity.get("rcsb_polymer_entity") or {})
                           .get("pdbx_description") or "")
            if description:
                self.descriptions.append(description)
            entity_organism = ""
            for organism in entity.get("rcsb_entity_source_organism") or []:
                pair = (organism.get("ncbi_scientific_name") or "",
                        organism.get("ncbi_taxonomy_id"))
                entity_organism = entity_organism or pair[0]
                if pair[0] and pair not in self.organisms:
                    self.organisms.append(pair)
            mutation_count = (entity.get("entity_poly") or {}).get("rcsb_mutation_count") or 0
            ids = entity.get("rcsb_polymer_entity_container_identifiers") or {}
            self.chains.extend(ids.get("auth_asym_ids") or [])
            if (entity.get("entity_poly") or {}).get("rcsb_entity_polymer_type") == "Protein":
                self.protein_entities.append(list(ids.get("uniprot_ids") or []))
            for align in entity.get("rcsb_polymer_entity_align") or []:
                if (align.get("reference_database_name") or "").lower() != "uniprot":
                    continue
                accession = align.get("reference_database_accession")
                covered = 0
                for region in align.get("aligned_regions") or []:
                    start, length = region.get("ref_beg_seq_id"), region.get("length") or 0
                    if start and length:
                        regions.setdefault(accession, []).append((start, start + length - 1))
                        covered += length
                if covered >= largest.get(accession, 0):
                    largest[accession] = covered
                    self.organism_of[accession] = entity_organism
                self.mutations[accession] = max(self.mutations.get(accession, 0),
                                                mutation_count)
        # union of the covered residues, so two chains cut from one precursor
        # (insulin A + B chain) add up, but identical copies do not count twice
        for accession, spans in regions.items():
            total, last = 0, 0
            for start, end in sorted(spans):
                start = max(start, last + 1)
                if end >= start:
                    total += end - start + 1
                    last = end
            self.coverage[accession] = total

    def method_short(self):
        short = []
        for method in self.methods:
            upper = method.upper()
            if "X-RAY" in upper:
                short.append("X-ray")
            elif "ELECTRON MICROSCOPY" in upper:
                short.append("EM")
            elif "NMR" in upper:
                short.append("NMR")
            elif "NEUTRON" in upper:
                short.append("Neutron")
            elif "ELECTRON CRYSTALLOGRAPHY" in upper:
                short.append("ED")
            else:
                short.append(method.title())
        return "+".join(dict.fromkeys(short)) or "?"

    def organism_text(self, accession=None):
        """Organism of the chain carrying ``accession`` (or of the first natural chain)."""
        if accession and self.organism_of.get(accession):
            return self.organism_of[accession]
        natural = [name for name, _ in self.organisms if "synthetic" not in name.lower()]
        if natural:
            return natural[0]
        return self.organisms[0][0] if self.organisms else "?"


def entries(pdb_ids):
    """Metadata for many entries (batched GraphQL).  Returns {pdb id: Entry}."""
    result = {}
    ids = [i.upper() for i in pdb_ids]
    for start in range(0, len(ids), 50):
        chunk = ids[start:start + 50]
        data = http.request_json(GRAPHQL_URL, {"query": _ENTRY_QUERY, "variables": {"ids": chunk}},
                                 ttl=7 * cache.DAY, timeout=20) or {}
        for record in ((data.get("data") or {}).get("entries") or []):
            if record:
                entry = Entry(record)
                result[entry.pdb_id] = entry
    return result


# ----------------------------------------------------------------- Tier 3
_SEQUENCE_QUERY = ("query($id:String!){entry(entry_id:$id){polymer_entities{"
                   "entity_poly{pdbx_seq_one_letter_code_can rcsb_entity_polymer_type} "
                   "rcsb_polymer_entity_container_identifiers{auth_asym_ids}}}}")


def chain_sequences(pdb_id):
    """{author chain id: full deposited protein sequence} for an entry."""
    data = http.request_json(GRAPHQL_URL, {"query": _SEQUENCE_QUERY,
                                           "variables": {"id": pdb_id.upper()}},
                             ttl=30 * cache.DAY, timeout=20) or {}
    result = {}
    entry = (data.get("data") or {}).get("entry") or {}
    for entity in entry.get("polymer_entities") or []:
        poly = entity.get("entity_poly") or {}
        if poly.get("rcsb_entity_polymer_type") != "Protein":
            continue
        sequence = (poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "")
        chains = (entity.get("rcsb_polymer_entity_container_identifiers") or {}).get(
            "auth_asym_ids") or []
        for chain in chains:
            result[chain] = sequence
    return result


def search_sequence(sequence, identity=0.3, rows=100):
    """Sequence-similarity search (MMseqs2 at RCSB).  Hits are grouped into
    clusters of >= 95% identical sequences and one representative per cluster
    is returned, so 300 structures of the same protein do not hide its
    relatives.  Returns (total matching chains, [(entity id like '104M_1',
    identity 0..1, e-value)]), best first."""
    payload = {
        "query": {"type": "terminal", "service": "sequence",
                  "parameters": {"evalue_cutoff": 0.1, "identity_cutoff": identity,
                                 "sequence_type": "protein", "value": sequence}},
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": rows},
                            "results_verbosity": "verbose",
                            "scoring_strategy": "sequence",
                            "group_by": {"aggregation_method": "sequence_identity",
                                         "similarity_cutoff": 95,
                                         "ranking_criteria_type": {"sort_by": "score",
                                                                   "direction": "desc"}},
                            "group_by_return_type": "representatives"},
    }
    data = http.request_json(SEARCH_URL, payload, ttl=cache.DAY, timeout=25)
    if not data:
        return 0, []
    hits = []
    for item in data.get("result_set") or []:
        best_identity, evalue = 0.0, None
        for service in item.get("services") or []:
            for node in service.get("nodes") or []:
                for context in node.get("match_context") or []:
                    if (context.get("sequence_identity") or 0) >= best_identity:
                        best_identity = context.get("sequence_identity") or 0
                        evalue = context.get("evalue")
        hits.append((item.get("identifier", ""), best_identity, evalue))
    return data.get("total_count", len(hits)), hits


_ENTITY_QUERY = """query($ids:[String!]!){polymer_entities(entity_ids:$ids){
 rcsb_id
 rcsb_polymer_entity{pdbx_description}
 rcsb_entity_source_organism{ncbi_scientific_name}
 rcsb_polymer_entity_container_identifiers{entry_id auth_asym_ids uniprot_ids}
 entry{struct{title} exptl{method} rcsb_entry_info{resolution_combined}
       rcsb_accession_info{initial_release_date}}
}}"""
_ENTITY_QUERY = " ".join(_ENTITY_QUERY.split())


def entities(entity_ids):
    """Metadata for polymer entities ('4HHB_1'): {id: dict}."""
    result = {}
    ids = list(entity_ids)
    for start in range(0, len(ids), 100):
        chunk = ids[start:start + 100]
        data = http.request_json(GRAPHQL_URL, {"query": _ENTITY_QUERY,
                                               "variables": {"ids": chunk}},
                                 ttl=7 * cache.DAY, timeout=25) or {}
        for record in (data.get("data") or {}).get("polymer_entities") or []:
            if not record:
                continue
            ids_block = record.get("rcsb_polymer_entity_container_identifiers") or {}
            entry = record.get("entry") or {}
            resolutions = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
            organisms = [o.get("ncbi_scientific_name") for o in
                         record.get("rcsb_entity_source_organism") or [] if o]
            methods = [m.get("method", "") for m in entry.get("exptl") or []]
            result[record.get("rcsb_id")] = {
                "pdb_id": (ids_block.get("entry_id") or "").upper(),
                "chains": ids_block.get("auth_asym_ids") or [],
                "accessions": ids_block.get("uniprot_ids") or [],
                "description": (record.get("rcsb_polymer_entity") or {}).get(
                    "pdbx_description") or "",
                "organism": (organisms[0] if organisms else None) or "?",
                "title": ((entry.get("struct") or {}).get("title") or "").strip(),
                "resolution": min(resolutions) if resolutions else None,
                "method": methods[0] if methods else "",
            }
    return result
