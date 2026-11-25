"""Turn a search.Result into the text PyMOL prints.  Plain ASCII on purpose,
so it displays correctly in every PyMOL console on every operating system."""
from __future__ import annotations

import re
import textwrap

from . import rank

LINE = "-" * 78


def _cut(text, width):
    text = " ".join(str(text).split())
    return text if len(text) <= width else text[:width - 2] + ".."


def _ligands(entry, width=14):
    if not entry.ligands:
        return "-"
    shown = ",".join(entry.ligands[:3])
    if len(entry.ligands) > 3:
        shown += "+%d" % (len(entry.ligands) - 3)
    if len(shown) > width:          # too wide: show fewer codes instead of cutting one
        shown = "%s+%d" % (",".join(entry.ligands[:2]), len(entry.ligands) - 2)
    return _cut(shown, width)


def table(hits):
    """Two lines per structure, so everything fits in an 80-column window."""
    rows = [" #  PDB    Score  Method  Res.     Year  Organism            Chains  Ligands"]
    for index, hit in enumerate(hits, 1):
        entry = hit.entry
        resolution = "%.2f A" % entry.resolution if entry.resolution is not None else "  -   "
        star = "*" if hit.classic else " "
        rows.append("%2d  %-4s%s  %5.1f  %-6s  %-7s  %4s  %-18s  %6s  %s" % (
            index, entry.pdb_id, star, hit.score, entry.method_short()[:6], resolution,
            entry.year or "?", _cut(entry.organism_text(hit.accession), 18),
            entry.chain_count or "?", _ligands(entry, 11)))
        rows.append("           %s" % _cut(entry.title, 67))
    if any(hit.classic for hit in hits):
        rows.append("    * = classic / textbook structure (Prospector's curated list)")
    rows.append("    Score = 0-100, higher is better (see below for what it is made of).")
    return "\n".join(rows)


def resolved_line(result):
    """How the typed name was understood, in plain words."""
    resolution = result.resolution
    if resolution is None or resolution.group is None:
        return ""
    group = resolution.group
    lead = group.lead
    typed = resolution.name
    if lead.gene and lead.gene.lower() != typed.lower():
        what = "the gene %s = UniProt %s" % (lead.gene, lead.accession)
    else:
        what = "UniProt %s" % lead.accession
    line = textwrap.fill("'%s' is %s (%s, %s)." % (typed, what, lead.name, lead.organism),
                         78, initial_indent="   ", subsequent_indent="   ")
    others = [p for p in group.members if p.accession != lead.accession and p.pdb]
    if others:
        named = ", ".join("%s %s" % (p.gene or p.name, p.accession) for p in others[:3])
        more = " +%d more" % (len(others) - 3) if len(others) > 3 else ""
        line += ("\n   Also counted: the same protein in other species, or its partner "
                 "chains\n   (%s%s)." % (named, more))
    return line


def facts(hit, result=None):
    """Two or three plain lines describing a structure: what it is, how it was solved."""
    entry = hit.entry
    lines = []
    protein = result.protein if result is not None else None
    if protein is not None and result.mode == "protein":
        lines.append("%s%s, %s" % (protein.name,
                                   " (gene %s)" % protein.gene if protein.gene else "",
                                   entry.organism_text(hit.accession)))
    else:
        lines.append(entry.organism_text(hit.accession))
    resolution = "%.2f A" % entry.resolution if entry.resolution is not None \
        else "no resolution value"
    chains = entry.chain_count or "?"
    ligands = "ligands: %s" % ", ".join(entry.ligands[:5]) if entry.ligands else "no ligands"
    lines.append("%s, %s, %s | %s chain%s | %s" % (
        entry.method_short(), resolution, entry.year or "?", chains,
        "" if chains == 1 else "s", ligands))
    return lines


def one_line(hit):
    """A one-line summary of a structure."""
    entry = hit.entry
    resolution = "%.2f A" % entry.resolution if entry.resolution is not None else "no resolution"
    ligands = ", ligands: %s" % ",".join(entry.ligands[:5]) if entry.ligands else ""
    return "%s: %s | %s | %s, %s, %s | %s chain(s)%s" % (
        entry.pdb_id, _cut(entry.title, 70), entry.organism_text(hit.accession),
        entry.method_short(),
        resolution, entry.year or "?", entry.chain_count or "?", ligands)


def molecule_lines(molecule, label):
    """The answer-first lines for a small molecule."""
    if molecule.source == "PubChem":
        out = ["%s%s  -  %s, from PubChem (compound %d)"
               % (label, molecule.title(), molecule.kind(), molecule.cid)]
    else:
        out = ["%s%s  -  %s, PDB code %s" % (label, molecule.title(), molecule.kind(),
                                              molecule.comp_id)]
    if molecule.nickname:
        out.append("             official name: %s" % molecule.name)
    weight = "%.1f g/mol" % molecule.weight if molecule.weight else "mass unknown"
    atoms = ("%d atoms (%d without hydrogens)" % (molecule.atoms, molecule.heavy_atoms)
             if molecule.atoms else "%d atoms without hydrogens" % molecule.heavy_atoms)
    out.append("             %s | %s | %s" % (molecule.formula or "formula unknown",
                                               weight, atoms))
    if molecule.how and not molecule.exact:
        out.append("NOTE         %s." % molecule.how)
    if molecule.description:
        first = re.split(r"(?<=[a-z0-9\)])\.\s", molecule.description)[0]
        out.append("WHAT IT IS   %s." % _cut(first.rstrip("."), 300))
    return [textwrap.fill(line, 78, subsequent_indent=" " * 13) for line in out]


def choices(result):
    """The commands behind the numbers of a numbered list, in order: typing 1
    in PyMOL runs the first one.  Empty when nothing was numbered."""
    query = result.query
    filters = " ".join(query.filter_words()) if query is not None else ""
    if result.kind == "ambiguous" and result.resolution is not None:
        groups = result.resolution.choices
        # the UniProt ID: exactly the protein on that line (MYH7 alone could
        # load the cow version when the line says human)
        return [" ".join(("grab %s %s" % (g.lead.accession, filters)).split())
                for g in groups]
    if result.kind == "results" and len(result.hits) > 1:
        return ["grab %s" % hit.entry.pdb_id for hit in result.hits]
    return []


def render(result, command_name="find"):
    out = []
    query = result.query
    typed = result.raw
    out.append("")
    out.append("Prospector: %s %s" % (command_name, typed))

    if result.kind != "results":
        for note in result.notes:
            out.append("Note: " + note)

    if result.kind == "error":
        out.append("")
        out.append(textwrap.fill(result.message, 78))

    elif result.kind == "results":
        top = result.top
        # 1. the answer
        out.append("")
        label = "THIS IS      " if result.mode == "id" else "BEST MATCH   "
        if query.pick > 1:
            label = "YOUR PICK    #%d of the list below:  " % query.pick
        out.append("%s%s - %s" % (label, top.pdb_id, _cut(top.entry.title, 58)))
        for line in facts(top, result):
            out.append("             " + line)
        if result.mode == "id":
            out.append("TO LOAD IT   grab %s" % top.pdb_id)
        else:
            again = " ".join([query.text] + query.filter_words())
            out.append("TO LOAD IT   grab %s" % again)
            out.append("             (PyMOL's own command for this entry:  fetch %s)"
                       % top.pdb_id)
        if result.mode == "id":
            out.append("")
            out.append("Run this:   %s" % result.command())
            out.append(LINE)
            return "\n".join(out)

        # 2. how the search was understood
        out.append("")
        out.append("HOW I READ YOUR SEARCH")
        for note in result.notes:
            out.append(textwrap.fill(note, 78, initial_indent="   ",
                                     subsequent_indent="   "))
        if result.mode == "protein":
            out.append(resolved_line(result))
            out.append("   %d structure%s of it %s in the PDB%s."
                       % (result.total, "" if result.total == 1 else "s",
                          "passes your filters" if query.has_filters() and result.total == 1
                          else "pass your filters" if query.has_filters()
                          else "exists" if result.total == 1 else "exist",
                          "" if result.total <= len(result.hits)
                          else " - these are the %d best" % len(result.hits)))
        else:
            if not result.notes:
                out.append("   No single protein has this name, so I searched the text of "
                           "every PDB entry.")
            out.append("   A more specific name gives a better answer.")
            out.append("   %d entr%s mention%s it%s."
                       % (result.total, "y" if result.total == 1 else "ies",
                          "s" if result.total == 1 else "",
                          "" if result.total <= len(result.hits)
                          else " - these are the best %d" % len(result.hits)))
        # 3. the alternatives
        out.append("")
        out.append(table(result.hits))
        out.append("")
        # 4. the reasoning
        out.append(rank.explain(result.hits))
        if result.mode == "protein":
            out.append("")
            out.append("Not the protein you meant?  Search all PDB text instead with:")
            out.append("   find %s mode=text" % query.text)
        if len(result.hits) > 1:
            other = 2 if query.pick == 1 else 1
            out.append("")
            out.append("Want a different row?  Just type its number, e.g.  %d" % other)
        out.append("")
        out.append("Run this:   %s" % result.command())

    elif result.kind == "ambiguous":
        resolution = result.resolution
        out.append("")
        if getattr(resolution, "family", False):
            out.append("'%s' is a family of many proteins. Its best-studied members - "
                       "pick one\nand run the command on its right:" % resolution.name)
        else:
            out.append("'%s' could mean several different proteins. I won't guess - "
                       "pick one\nand run the command on its right:" % resolution.name)
        out.append("")
        filters_text = " ".join(query.filter_words())
        for index, group in enumerate(resolution.choices, 1):
            command = "%s %s" % ("grab" if command_name == "grab" else "find",
                                 group.suggested_query(resolution.choices))
            if filters_text:
                command += " " + filters_text
            out.append("  %d. %s  -  %d structures" % (index, _cut(group.label(), 52),
                                                       group.count))
            out.append("     %-34s ->  %s" % (_cut(group.organisms_text(), 34), command))
        if resolution.hidden > 0:
            out.append("  (+%d smaller matches not shown)" % resolution.hidden)
        out.append("")
        out.append("TO LOAD ONE  just type its number:  1   (or 2, 3 ...)")
        out.append("Or search every PDB entry that mentions the word:  find %s mode=text"
                   % query.text)

    elif result.kind == "suggest":
        out.append("")
        out.append(result.message)
        out.append("")
        out.append("Run this:   %s" % result.command())

    elif result.kind == "nofilter":
        out.append("")
        out.append(textwrap.fill(result.message, 78))
        out.append("")
        out.append("How I read your search:")
        out.append(resolved_line(result))

    elif result.kind == "alphafold":
        protein, model = result.protein, result.model
        out.append("")
        out.append("NO EXPERIMENTAL STRUCTURE of %s exists in the PDB." % protein.label())
        confidence = ""
        if model.mean_plddt is not None:
            confidence = ", mean confidence %.0f out of 100" % model.mean_plddt
        out.append("BUT a computer prediction is available: AlphaFold model %s"
                   % model.model_id)
        out.append("             (residues %s-%s%s)." % (model.start, model.end, confidence))
        out.append("             A prediction is a good guide to the fold, but it is not "
                   "an experiment.")
        out.append("TO LOAD IT   alphafold %s" % model.accession)
        out.append("")
        out.append("How I read your search:")
        out.append(resolved_line(result))
        out.append("")
        out.append("Run this:   %s" % result.command())

    elif result.kind == "molecule":
        out.append("")
        out.extend(molecule_lines(result.molecule, "BEST MATCH   "))
        out.append("TO LOAD IT   grab %s" % typed)
        molecule = result.molecule
        if molecule.entries:
            out.append("IN PROTEINS  bound in %d PDB structure%s; to see one:  grab protein "
                       "ligand=%s" % (molecule.entries, "" if molecule.entries == 1 else "s",
                                      molecule.comp_id))
        if molecule.source == "PubChem":
            out.append("NOTE         not in the PDB: no structure has it bound. The shape "
                       "is a computed one from PubChem.")
        out.append("NOT THIS?    proteins named after it:  find %s kind=protein"
                   % (query.text if query else typed))
        for other in molecule.alternatives:
            out.append("             or the molecule %s (%s):  grab %s kind=molecule"
                       % (_cut(other.title(), 34), other.comp_id, other.comp_id))

    elif result.kind == "none":
        out.append("")
        out.append(result.message)
        if "small molecule" not in result.message:
            out.append("Tips: check the spelling; try the gene name (TP53) or a UniProt ID "
                       "(P04637);\n      remove filters; or use fewer words.")

    if result.kind in ("results", "ambiguous") and result.mode != "id" and query \
            and query.kind == "auto" and not query.has_filters() \
            and (result.resolution is None or result.resolution.top_tier < 3) \
            and not getattr(result.resolution, "family", False) \
            and " " not in query.text.strip():
        # a name no protein owns outright could also be a small molecule
        out.append("MOLECULE?    if you meant the small molecule of that name:  "
                   "molecule %s" % query.text)
    out.append(LINE)
    return "\n".join(line for line in out if line is not None)
