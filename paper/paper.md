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
disease variants and the meaning of each colour. Students can recolour and
redraw the structure, measure and interpret distances, see a mutation next to
the original protein, compare structures, and build models of DNA, RNA,
peptides, membranes and micelles. Requests can also be typed as plain
sentences, and every feature is a button in a panel docked inside PyMOL. The
plugin comes with a one-hour laboratory worksheet on haemoglobin.

# Statement of need

PyMOL is one of the most widely used programs for viewing biomolecules and is
common in undergraduate teaching. For a beginner, however, it raises three
barriers before any biology is learned. A structure must be requested by its
four-character PDB identifier, so students leave the program to search a
database and choose among dozens or hundreds of entries. Almost every action
needs PyMOL's command and selection syntax. And the program shows a structure
without explaining it: which residues bind the ligand, whether a distance
means a bond, or what a mutation changes.

Web resources address parts of this. Proteopedia [@hodis2008proteopedia]
pairs structures with explanatory text, and the Mol\* viewer
[@sehnal2021molstar] makes structures easy to open in a browser. Neither brings
name-based searching, explanation and experimentation into the desktop viewer
that many courses and research groups already use. Prospector is aimed at
students in the first years of biology, biochemistry or biotechnology, and at
teachers who want a lab session to be about proteins rather than syntax.

# Learning objectives

With Prospector, students can:

1. find and load a protein structure from its name, and judge whether it is a
   suitable one (method, resolution, coverage, organism);
2. locate functional sites and bound molecules, and relate them to the
   protein's function;
3. interpret measured distances (bonds, hydrogen bonds, contacts);
4. connect a change in sequence to its structural and medical consequences;
5. distinguish changes within a chain from changes between chains;
6. gradually learn the PyMOL commands behind each step.

# Instructional design

Every command answers in the same plain-language layout: what is now on
screen, what it means, and what to try next. Ambiguous names (such as
"myosin") give a numbered list instead of a guess, and models that are not
experimental structures, such as AlphaFold predictions
[@jumper2021alphafold; @varadi2024afdb] or idealised membranes, are labelled as
such. Experimenting is safe: `whatif` changes a residue in a copy of the
protein and shows the original and the changed version side by side, with the
differences in charge, size, contacts and hydrogen bonds
(\autoref{fig:whatif}), and `back` undoes any view. A teaching mode prints the
plain PyMOL commands used at each step, so students who start with Prospector
can move on to PyMOL itself. Annotations come from UniProt, mapped to the
structure through SIFTS [@dana2019sifts], and small molecules not in the PDB
from PubChem [@kim2023pubchem].

![The sickle-cell mutation shown by `whatif haemoglobin B/6 V`: the original
glutamate (left, green) and the modelled valine (right, magenta) in the beta
chain of deoxyhaemoglobin (PDB 1A3N). Prospector reports that a negative
charge is lost and an oily side chain is placed on the protein
surface.\label{fig:whatif}](figure_whatif.png)

The repository includes a one-hour worksheet for first-year students, with
teacher notes and an answer key. Students load deoxyhaemoglobin by name
[@tame2000deoxy]; locate the heme pocket and measure that the proximal
histidine is bonded to the iron (2.2 Å) while the distal histidine is not
(4.2 Å; \autoref{fig:pocket}); model the Glu6Val change of sickle-cell anaemia and relate the new
surface hydrophobic patch to fibre formation; and compare the deoxy and oxy
forms [@shaanan1983oxy], finding that each alpha–beta pair barely changes
(RMSD 0.5 Å) while the whole tetramer does (2.3 Å), the quaternary change
behind cooperative oxygen binding (\autoref{fig:compare}). Every value in the
answer key was produced with the plugin itself.

![The heme pocket of the alpha chain of deoxyhaemoglobin (PDB 1A3N), as
measured with `measure A/87 to the iron` and `measure A/58 to the iron`. The
proximal His87 is bonded to the iron (2.2 Å); the distal His58 is not
(4.2 Å), leaving the site where oxygen binds.\label{fig:pocket}](figure_pocket.png)

![Deoxy (T state, PDB 1A3N, grey) and oxy (R state, PDB 1HHO, coloured)
haemoglobin after superimposing one alpha–beta pair (right); colour shows how
far each residue moved. The superimposed pair barely changes (mean 0.7 Å)
while the other pair is displaced (mean 4.9 Å): the change lies in how the
pairs sit against each other. The oxy tetramer was rebuilt with `build
assembly 1HHO`; the image was rendered in PyMOL.\label{fig:compare}](figure_compare.png)

# Software

Prospector is written in Python using only the standard library, so it
installs through PyMOL's own Plugin Manager with no further dependencies; it
is designed for both the open-source and the commercial PyMOL. A pure-Python core (database
clients, name resolution, ranking and caching) is kept separate from the
PyMOL layer (41 commands, a background worker and a Qt panel;
\autoref{fig:panel}). Network
requests run in the background with timeouts, so PyMOL never freezes, and
answers are cached so that a prepared class still works without a network.
Plain sentences are first interpreted by built-in rules; only those the rules
cannot read are sent, if the user supplies their own API key, to an optional
large language model, whose suggested commands are checked against a fixed
list and always shown before they run. An offline test suite of 177 tests and
385 example sentences checks the behaviour.

![Prospector's panel docked inside PyMOL during the sickle-cell experiment of
\autoref{fig:whatif}. (A) The 3D view shows the original haemoglobin (left)
and the changed copy (right) with the same camera. (B) Requests can be typed
as plain sentences. (C) Five tabs follow the order of a lab session: find,
understand, look, experiment, keep. (D) The student picks the residue (B/6)
and what it becomes (valine) from a list; the change runs at once. (E) The
result is explained in plain words: a negative charge is lost and a
water-avoiding side chain is placed on the protein
surface.\label{fig:panel}](figure_panel.png)

# Reuse and modification

The worksheet is written to be adapted to other proteins: the same sequence of
commands (`grab`, `pocket`, `measure`, `explain`, `whatif`, `compare`) works for
any structure in the PDB. Prospector, the worksheet and its teacher notes are
released under the MIT licence, and new teaching materials are welcome
through the repository's contribution guidelines.

# Story of the project

Prospector began when the author, a first-year biotechnology student, wanted
to find and understand protein structures in PyMOL without first having to
learn its syntax. It was developed step by step: each version was used in
PyMOL, and what was confusing to a new user was changed. For example,
distance measurement became a deliberate two-click mode after clicks started
measurements by accident, and the choices offered for an ambiguous name can
now be loaded by typing their number.

# Software availability

- **Source code and releases:** https://github.com/lorenzopattofatto10-wq/prospector-for-pymol
- **Archive and DOI (all versions):** https://doi.org/10.5281/zenodo.24355535
- **Licence:** MIT (software)
- **Requirements:** PyMOL 2 or later with Python 3; no further dependencies
- **Teaching material:** `docs/worksheets/` in the repository (student sheet, teacher notes)

# Acknowledgements

Prospector relies on the open data and services of the RCSB Protein Data Bank,
UniProt, the AlphaFold Protein Structure Database, PDBe and PubChem.

# References
