"""The `about` info card (pure Python): labelled sections, aligned fields."""
from __future__ import annotations

import re
import textwrap

from . import rcsb, uniprot

WIDTH = 78
RULE = "=" * WIDTH
LABEL = 13          # width of the field-name column


def _section(title):
    return ["", " " + title, " " + "-" * (WIDTH - 2)]


def _field(label, value):
    """'   Label        value', wrapped so continuation lines stay aligned."""
    first = "   %-*s" % (LABEL, label)
    return textwrap.fill(str(value), width=WIDTH, initial_indent=first,
                         subsequent_indent=" " * len(first)) or first.rstrip()


def _cut(text, width):
    text = " ".join(str(text).split())
    return text if len(text) <= width else text[:width - 2] + ".."


def _short_function(text, limit=260):
    """The first sentence or two of a UniProt function text."""
    text = " ".join(text.split())
    if len(text) <= limit:
        return text
    sentences = re.split(r"(?<=[.!?])\s+", text)
    out = ""
    for sentence in sentences:
        if out and len(out) + len(sentence) + 1 > limit:
            break
        out = (out + " " + sentence).strip()
    if not out or len(out) > limit + 60:
        out = text[:limit].rsplit(" ", 1)[0] + " ..."
    return out


def _method(detail):
    methods = [m.get("method", "") for m in detail.get("exptl") or []]
    return " + ".join(m.capitalize() for m in methods) or "unknown"


def _protein_block(out, protein, record, heading):
    out.extend(_section(heading))
    out.append(_field("Name", protein.name))
    out.append(_field("Gene", protein.gene or "-"))
    out.append(_field("Organism", protein.organism))
    out.append(_field("UniProt", "%s  (%d amino acids)" % (protein.accession, protein.length)))
    function = uniprot.function_text(record)
    if function:
        out.append(_field("Function", _short_function(function)))


def entry_card(pdb_id):
    """Text card for a PDB entry, or None if the entry does not exist."""
    detail = rcsb.entry_detail(pdb_id)
    if not detail:
        return None
    out = [RULE]
    title = " ".join(((detail.get("struct") or {}).get("title") or "").split())
    out.append(textwrap.fill(title, width=WIDTH,
                             initial_indent=" %s  |  " % detail.get("rcsb_id", pdb_id.upper()),
                             subsequent_indent=" " * 10))
    out.append(RULE)

    # ---- experiment
    info = detail.get("rcsb_entry_info") or {}
    resolutions = info.get("resolution_combined") or []
    refine = (detail.get("refine") or [{}])[0] or {}
    r_work, r_free = refine.get("ls_R_factor_R_work"), refine.get("ls_R_factor_R_free")
    released = ((detail.get("rcsb_accession_info") or {}).get("initial_release_date")
                or "")[:4] or "?"
    out.extend(_section("EXPERIMENT"))
    out.append(_field("Method", _method(detail)))
    out.append(_field("Resolution", "%.2f A" % min(resolutions) if resolutions
                      else "not applicable"))
    if detail.get("refine"):
        out.append(_field("R-work", "%.3f" % r_work if r_work is not None
                          else "not reported"))
        out.append(_field("R-free", "%.3f" % r_free if r_free is not None
                          else "not reported"))
    out.append(_field("Released", released))

    # ---- proteins (via UniProt)
    accessions = []
    for entity in detail.get("polymer_entities") or []:
        ids = entity.get("rcsb_polymer_entity_container_identifiers") or {}
        for accession in ids.get("uniprot_ids") or []:
            if accession not in accessions:
                accessions.append(accession)
    shown = accessions[:4]
    for number, accession in enumerate(shown, 1):
        protein, record = uniprot.details(accession)
        if protein is None:
            continue
        heading = "PROTEIN" if len(shown) == 1 else "PROTEIN %d of %d" % (number,
                                                                          len(accessions))
        _protein_block(out, protein, record, heading)
    if len(accessions) > 4:
        out.append("   (+%d more proteins in this entry)" % (len(accessions) - 4))
    if not accessions:
        out.extend(_section("PROTEIN"))
        out.append("   No chain of this entry is linked to UniProt (DNA/RNA, a designed "
                   "protein\n   or an antibody).")

    # ---- chains
    out.extend(_section("CHAINS"))
    row = "   %-7s %-29s %-7s %-12s %s"
    out.append(row % ("Chain", "Molecule", "Length", "UniProt res.", "Organism"))
    for entity in detail.get("polymer_entities") or []:
        ids = entity.get("rcsb_polymer_entity_container_identifiers") or {}
        poly = entity.get("entity_poly") or {}
        chains = ",".join(ids.get("auth_asym_ids") or []) or "?"
        description = (entity.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "?"
        organisms = [o.get("ncbi_scientific_name") for o in
                     entity.get("rcsb_entity_source_organism") or [] if o]
        kind = poly.get("rcsb_entity_polymer_type") or "Polymer"
        length = poly.get("rcsb_sample_sequence_length") or 0
        unit = "aa" if kind == "Protein" else "nt"
        span = "-"
        for align in entity.get("rcsb_polymer_entity_align") or []:
            if (align.get("reference_database_name") or "").lower() != "uniprot":
                continue
            regions = align.get("aligned_regions") or []
            if regions:
                first = min(r["ref_beg_seq_id"] for r in regions)
                last = max(r["ref_beg_seq_id"] + r["length"] - 1 for r in regions)
                span = "%d-%d" % (first, last)
        mutations = poly.get("rcsb_mutation_count") or 0
        if kind != "Protein" and description.upper().startswith(("DNA (", "RNA (")):
            description = kind          # the deposited name is just the base sequence
        molecule = _cut(description, 19 if mutations else 29)
        if mutations:
            molecule += " (%d mut.)" % mutations
        out.append(row % (_cut(chains, 7), molecule, "%d %s" % (length, unit), span,
                          _cut(organisms[0] if organisms else "-", 17)))

    # ---- ligands
    out.extend(_section("LIGANDS"))
    names = rcsb.ligand_names(detail)
    if names:
        out.append("   %-6s %-44s %s" % ("Code", "Name", "Formula"))
        for code in sorted(names, key=lambda c: (c in rcsb.ADDITIVES, c)):
            name, formula = names[code]
            name = name.capitalize() if name.isupper() else name
            if code in rcsb.ADDITIVES:
                name = _cut(name, 30) + " [additive]"
            out.append("   %-6s %-44s %s" % (code, _cut(name, 44), formula))
    else:
        out.append("   none")

    # ---- citation
    out.extend(_section("CITATION"))
    citation = detail.get("rcsb_primary_citation") or {}
    if citation:
        authors = citation.get("rcsb_authors") or []
        who = authors[0] + (" et al." if len(authors) > 1 else "") if authors else "?"
        journal = citation.get("rcsb_journal_abbrev") or ""
        if citation.get("journal_volume"):
            journal += " %s" % citation["journal_volume"]
        if citation.get("page_first"):
            journal += ":%s" % citation["page_first"]
        out.append(_field("Authors", who))
        out.append(_field("Year", citation.get("year") or "?"))
        out.append(_field("Title", (citation.get("title") or "?").rstrip(".")))
        out.append(_field("Journal", journal or "-"))
        doi = citation.get("pdbx_database_id_DOI")
        out.append(_field("DOI", "https://doi.org/" + doi if doi else "none recorded"))
        if citation.get("pdbx_database_id_PubMed"):
            out.append(_field("PubMed", citation["pdbx_database_id_PubMed"]))
    else:
        out.append("   none recorded (the structure may be unpublished)")
    out.append(RULE)
    return "\n".join(out)


def molecule_card(molecule):
    """Text card for a small molecule from the PDB's chemical dictionary."""
    out = [RULE, " %s  |  %s  |  small molecule" % (
        molecule.comp_id or "PubChem %d" % molecule.cid, molecule.title()), RULE]
    out.extend(_section("MOLECULE"))
    out.append(_field("Name", molecule.name))
    out.append(_field("Kind", molecule.kind()))
    out.append(_field("Formula", molecule.formula or "?"))
    if molecule.weight:
        out.append(_field("Mass", "%.1f g/mol" % molecule.weight))
    out.append(_field("Atoms", "%d (%d without hydrogens)"
                      % (molecule.atoms, molecule.heavy_atoms) if molecule.atoms
                      else "%d without hydrogens" % molecule.heavy_atoms))
    charge = molecule.charge
    out.append(_field("Charge", "none (neutral)" if not charge else "%+d" % charge))
    if molecule.chiral_centres:
        out.append(_field("Chirality", "%d chiral centre%s (atoms whose mirror image is a "
                          "different molecule)" % (molecule.chiral_centres,
                                                   "" if molecule.chiral_centres == 1
                                                   else "s")))
    out.append(_field("Rings", "has aromatic (flat, benzene-like) rings"
                      if molecule.aromatic_bonds else "no aromatic rings"))
    names = [n for n in molecule.synonyms if n.lower() != molecule.name.lower()][:5]
    if names:
        out.append(_field("Also called", "; ".join(names)))
    if molecule.description:
        out.extend(_section("WHAT IT IS"))
        out.extend(textwrap.wrap(molecule.description[:900], WIDTH - 3, initial_indent="   ", subsequent_indent="   "))
    if molecule.categories:
        out.extend(_section("CLASSES"))
        out.extend(textwrap.wrap("; ".join(molecule.categories[:10]), WIDTH - 3, initial_indent="   ", subsequent_indent="   "))
    if molecule.source == "PubChem":
        out.extend(_section("PROPERTIES"))
        if getattr(molecule, "iupac", ""):
            out.append(_field("IUPAC name", molecule.iupac))
        if molecule.xlogp is not None:
            out.append(_field("Likes water?", "XLogP %s  (below 0 = dissolves in water, "
                              "above 3 = oily, prefers membranes and fat)" % molecule.xlogp))
        if molecule.donors is not None and molecule.acceptors is not None:
            out.append(_field("H-bonds", "can give %s and accept %s hydrogen bonds"
                              % (molecule.donors, molecule.acceptors)))
        out.extend(_section("SOURCE"))
        out.append(_field("Database", "PubChem (NIH), compound %d%s"
                          % (molecule.cid, "; description from %s"
                             % molecule.description_source
                             if getattr(molecule, "description_source", "") else "")))
        out.append(_field("In the PDB", "not in the PDB's dictionary: no structure has "
                          "this molecule bound"))
        out.append(_field("3D shape", "a computed, energy-minimised shape from PubChem, "
                          "not a measurement"))
        out.append(RULE)
        return "\n".join(out)
    out.extend(_section("IN THE PDB"))
    if molecule.entries is not None:
        out.append(_field("Found in", "%d structure%s, bound to a protein or nucleic acid"
                          % (molecule.entries, "" if molecule.entries == 1 else "s")))
    out.append(_field("See one", "grab protein ligand=%s" % molecule.comp_id))
    out.append(_field("Source", "PDB Chemical Component Dictionary; description from "
                      "DrugBank where available"))
    out.append(_field("3D shape", "an idealised geometry computed for the dictionary, not "
                      "a measurement of this molecule on its own"))
    out.append(RULE)
    return "\n".join(out)


def model_card(model):
    """Text card for an AlphaFold model."""
    out = [RULE, " %s  |  AlphaFold predicted model" % model.model_id, RULE]
    out.extend(_section("MODEL"))
    out.append(_field("Type", "a computed prediction, not an experiment"))
    out.append(_field("Version", model.version or "?"))
    out.append(_field("Residues", "%s-%s, one chain (A), numbered like UniProt"
                      % (model.start, model.end)))
    if model.mean_plddt is not None:
        out.append(_field("Confidence", "mean pLDDT %.0f out of 100" % model.mean_plddt))
    for label, value in model.fractions:
        if value is not None:
            out.append(_field("", "%3.0f%%  %s" % (100 * value, label)))
    protein, record = uniprot.details(model.accession)
    if protein is not None:
        _protein_block(out, protein, record, "PROTEIN")
        known = len(protein.pdb)
        out.append(_field("PDB entries", "%d experimental structure%s of this protein"
                          % (known, "" if known == 1 else "s") if known
                          else "none - no experimental structure exists"))
    out.extend(_section("CITATION"))
    out.append(_field("Method", "Jumper et al. (2021) Nature 596:583"))
    out.append(_field("Database", "Varadi et al. (2024) Nucleic Acids Res 52:D368"))
    out.append(_field("DOI", "https://doi.org/10.1038/s41586-021-03819-2"))
    out.append(RULE)
    return "\n".join(out)
