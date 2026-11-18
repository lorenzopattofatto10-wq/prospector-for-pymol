---
title: 'Prospector for PyMOL: finding, showing and explaining protein structures by name'
tags:
  - Python
  - PyMOL
  - structural biology
  - protein structure
  - education
authors:
  - name: Lorenzo Pattofatto
    orcid: 0009-0008-7977-2084
    affiliation: 1
affiliations:
  - name: Imperial College London, United Kingdom
    index: 1
date: 18 November 2026
bibliography: paper.bib
---

# Summary

Prospector for PyMOL is an open-source plugin for the molecular viewer PyMOL
[@pymol] that lets students work with protein structures by name instead of by
code. A student types `grab haemoglobin`; Prospector finds the protein in
UniProt [@uniprot2023], ranks its experimental structures in the Protein Data
Bank [@berman2000pdb], loads the best one and says why it was chosen. It then
explains what is on the screen: functional sites, bound molecules, known
disease variants and the meaning of each colour. Further commands recolour and
redraw the structure, measure distances and interpret them, show a mutation
next to the original protein, compare two structures, and build idealised
models of DNA, RNA, peptides, lipid membranes and micelles, alongside real
biological assemblies rebuilt from the PDB. Requests can also be written as
plain sentences (`ask where does oxygen bind`), and every feature is available
both as a command and as a button in a panel docked inside PyMOL. Because
Prospector can print the plain PyMOL commands behind each step, it doubles as
a way to learn PyMOL itself.

# Statement of need

PyMOL is one of the most widely used programs for viewing biomolecules, and it
is common in undergraduate biology and biochemistry teaching. For a beginner,
however, it raises three barriers before any biology is learned. First, a
structure must be requested by its four-character PDB identifier, so students
have to leave the program, search a database and choose among dozens or
hundreds of entries without knowing which are suitable. Second, almost every
action requires PyMOL's command and selection syntax. Third, the program shows
a structure but does not explain it: which residues bind the ligand, whether a
distance means a bond, or what a mutation changes.

Web resources address parts of this. Proteopedia [@hodis2008proteopedia]
pairs structures with explanatory text, and the Mol\* viewer
[@sehnal2021molstar] makes structures easy to open in a browser. Neither brings
name-based searching, explanation and experimentation into the desktop viewer
that many courses and research groups already use. Prospector is aimed at
students in their first years of biology, biochemistry or biotechnology, and
at teachers who want a lab session to be about proteins rather than about
software syntax.

# Features and design

Prospector adds 41 commands to PyMOL. The main groups are:

- **Find and load:** `find` and `grab` resolve a protein name, gene symbol,
  UniProt accession or PDB ID, score candidate structures on how well they
  match the name, sequence coverage, resolution, organism, method and age
  (with a small bonus for textbook structures), and explain the ranking. Ambiguous names (such as
  "myosin") produce a numbered list instead of a guess. Proteins without an
  experimental structure can be loaded as predictions from the AlphaFold
  Protein Structure Database [@jumper2021alphafold; @varadi2024afdb], clearly
  labelled as such.
- **Understand:** `about`, `tour`, `explain`, `features` and `pocket` combine
  the structure with UniProt annotations, mapped to the structure's residues
  through SIFTS [@dana2019sifts], to describe sites, domains and disease
  variants in plain words.
- **Look and select:** `look`, `viewas`, `scout` and `choose` change colours
  and representations, or select atoms described in words, showing the
  equivalent PyMOL selection.
- **Experiment:** `whatif` swaps a residue and shows the original and changed
  protein side by side with the changes in charge, size, contacts and hydrogen
  bonds (\autoref{fig:whatif}); `measure` interprets distances; `compare`
  superimposes two structures and reports where they differ.
- **Build:** `build` creates idealised nucleic acids, peptides, membranes and
  micelles, each marked as an idealised model, and real biological assemblies
  from the PDB. Small molecules not in the PDB can be loaded from PubChem
  [@kim2023pubchem].

![The sickle-cell mutation shown by `whatif haemoglobin B/6 V`: the original
glutamate (left, green) and the modelled valine (right, magenta) in the beta
chain of deoxyhaemoglobin (PDB 1A3N). Prospector reports that a negative
charge is lost and an oily side chain is placed on the protein
surface.\label{fig:whatif}](figure_whatif.png)

The code is written in Python using only the standard library, so the plugin
installs through PyMOL's own Plugin Manager with no further dependencies. A
pure-Python core (database clients, name resolution, ranking and caching) is
kept separate from the PyMOL layer (commands, background worker and Qt
panel). All network requests run in a background thread with timeouts, so
PyMOL never freezes, and answers are cached so that a prepared class keeps
working without a network connection. Plain-sentence requests are first
interpreted by built-in rules; only sentences the rules cannot read are sent,
if the user has supplied their own API key, to an optional large language
model, whose suggested commands are checked against a fixed list before they
run and are always shown to the user. The test suite (177 tests) runs offline
on recorded database answers, and a further set of 385 example sentences
checks how requests are understood.

# Use in teaching

The repository includes a one-hour laboratory worksheet for first-year
students, with teacher notes and an answer key, built around haemoglobin. In
four parts, students load deoxyhaemoglobin by name [@tame2000deoxy], locate
the heme pocket and measure that the proximal histidine is bonded to the iron
(2.2 Å) while the distal histidine is not (4.2 Å); model the Glu6Val change of
sickle-cell anaemia and relate the new surface hydrophobic patch to fibre
formation; and compare the deoxy and oxy forms [@shaanan1983oxy], finding that
each alpha–beta pair barely changes (RMSD 0.5 Å) while the whole tetramer does
(2.3 Å), the quaternary change underlying cooperative oxygen binding. Every
value in the answer key was produced with the plugin itself. The worksheet is
designed to be adapted to other proteins, and further materials are welcome
through the repository's contribution guidelines.

# AI usage disclosure

Generative AI was used extensively in this project. The plugin's code, tests,
documentation, the teaching worksheet and a first draft of this paper were
written by Claude (Anthropic), a large language model, working from the
author's instructions over many sessions. The author decided what Prospector
should do and how it should behave, tested each version in PyMOL, reported
problems found in use (for example, a protein name being misread as an amino
acid), requested changes, and reviewed and edited the documentation, the
worksheet and this paper. The author takes responsibility for the content of
the software and the paper.

# Acknowledgements

Prospector relies on the open data and services of the RCSB Protein Data Bank,
UniProt, the AlphaFold Protein Structure Database, PDBe and PubChem.

# References
