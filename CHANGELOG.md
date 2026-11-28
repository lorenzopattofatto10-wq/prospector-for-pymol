# Changelog - Prospector for PyMOL

## 2.9.8

Three problems found by recording a real PyMOL session (the worksheet's order:
pocket, measure, whatif, whatif off):

- **`whatif` shows two panels again after `pocket` or `measure`.** The dashed lines
  those commands leave each took a panel of their own, so students saw three or
  four squeezed panels. Everything other than the two structures is now hidden
  while the panels are shown.
- **`whatif off` puts the original back as it was:** its colours, solid ribbon
  and the dashes and other objects `whatif` hid. Before, the protein stayed grey
  and see-through, so `pocket` and `measure` afterwards looked washed out.
- **`viewas surface` draws the surface in PyMOL's window** (and `look
  hydrophobicity` its coloured surface). PyMOL left a surface undrawn when it was
  asked for right after the snapshot that `back` uses; it is now redrawn a
  moment later. Tests without a window could not see this.

## 2.9.7

- **New command `piece`: part of a structure as a separate object.** `piece A 19-40`,
  `piece chymotrypsinogen A/19-20`, `piece residues 19 to 40 of chain A`,
  `piece chain B`, several ranges (`piece A 19-20 30-40`), PyMOL's own
  `/4CHA//B/40-60`, and an optional colour (`piece A 19-40 red`). The new object
  (e.g. `4CHA_B19-40`) appears in PyMOL's object list and can be coloured, hidden
  or shown on its own; the original keeps every atom (the residues are only hidden
  there, so they are not drawn twice). Short pieces are drawn as sticks. A missing
  chain or residue range is named together with the ones that exist; numbers that
  exist only in UniProt's numbering are converted, and a difference between the
  file's numbers and UniProt's is pointed out. Also in the panel (Look tab,
  **Separate part**) and in `ask` ("make residues 19 to 40 of chain A a separate
  object"). `back` removes it.
- **Every big family gives a choice.** `grab kinase`, `protease`, `phosphatase`,
  `GPCR`, `ion channel`, `potassium channel`, `transporter`, `ABC transporter`,
  `ATPase`, `nuclear receptor`, `collagen`, `keratin`, `tubulin`, `actin`,
  `histone`, `immunoglobulin`, `heat shock protein`, `chaperone`, `globin`,
  `interleukin`, `interferon`, `ubiquitin ligase`, `dehydrogenase`, `synthase`,
  `reductase` and others used to load one structure from a text search (or, for
  GPCR, a plant protein that happens to share the abbreviation). They now list
  the best-studied members - human ones first, most structures first - and a
  number loads one. Members are proteins whose name the word heads, proteins
  UniProt files in that family, or (where that makes a far better list) proteins
  with that UniProt keyword. Lists that already existed lost proteins that only
  mention the word (Pannexin-1 under caspase, CD47 under integrin).
- Complexes load a structure of the whole complex: `nucleosome` and `spliceosome`
  no longer offer a list of their parts. New textbook entries: the yeast
  spliceosome (3JB9) and photosystems I and II (1JB0, 3WU2).
- Unchanged on purpose: single proteins (`haemoglobin`, `p53`, `GFP`, `ubiquitin`)
  load at once, `insulin` and `lysozyme` still ask which one, searches with a filter
  (`grab kinase ligand=STI`) load the best match, and `antibody` and `zinc
  finger` (parts shared by thousands of unrelated proteins) still search the text.
- `info` fits the console again (one line was wrapping).

## 2.9.6

Found with `grab aquaporin`, then tested on 25 structures of every kind (enzymes,
channels, receptors, glycoproteins, a virus shell, ferritin, chaperonins):

- **The working molecule, always.** When a PDB file holds only part of the working
  molecule (one chain of the aquaporin tetramer, one of the 60 of a virus shell),
  or extra copies that only sit side by side in the crystal (two copies of a
  kinase that works alone), every command that loads a structure now loads the
  PDB's assembly instead and says so.
- **New command `unit whole|one|file`** (also in the panel's Look tab and in `ask`:
  "show the monomer", "back to the tetramer"): switch between the working
  molecule, a single chain, and exactly what the PDB file holds. The first copy keeps its chain names; the copies get the next free
  letters and share its UniProt annotations (features, explain, tour). PyMOL's own
  `fetch` still loads a file exactly as it is.
- **Everything that is not protein is explained**, for any structure: water,
  buffer and cryo-protectant, simple ions (which may be functional, like K+ in a
  potassium channel), detergents, membrane lipids, sugar chains (only when
  chemically bonded to the protein), and molecules that do not touch the protein.
  A file that holds extra copies packed in the crystal is pointed out too.
- **Family names give a list.** "aquaporin" now lists AQP1, AqpZ, AQP3 ...
  instead of quietly picking one structure; concepts such as "ribosome" still use
  the text search. A protein picked from a list is called by its gene name.
- **`build` tolerates typos:** `build asembly` works, and typing 8Y80 (zero) while
  8Y8O (letter O) is on screen uses the one on screen. Its help is laid out
  cleanly.
- Worksheet question 16 uses `fetch 1HHO, oxy_half` for the half-molecule file.

## 2.9.5

No change to how the plugin works. This version archives:

- the paper (`paper/`), with four figures, references and an AI usage disclosure;
- printable PDF versions of the haemoglobin worksheet and teacher notes;
- the author's ORCID and affiliation in the citation file;
- a CC BY 4.0 licence for the text, worksheets and figures (`LICENSE-docs.md`);
  the software stays under the MIT licence.

## 2.9.4

Found while testing the haemoglobin worksheet:

- `whatif haemoglobin B/6 V` (an address, then the new amino acid) now works,
  also as `B/6 to valine`.
- `explain` suggests a change in the structure's own numbering (`whatif 1A3N
  B/E6V`), not UniProt's, which can differ (haemoglobin's Glu6 is UniProt's 7).
- `measure A/87 to the iron` finds a metal atom inside a molecule, such as the
  iron of heme; the closest one is measured.
- `pocket oxygen` (and `ask where does oxygen bind`) with a protein on screen
  uses the molecule in that structure instead of downloading a lone O2; in a
  heme protein without bound O2 it shows the heme pocket and says why.
- `build assembly` names the chains A, B, C, D ... copy by copy, like a normal
  entry, so `compare` pairs the right chains (deoxy vs oxy haemoglobin now gives
  the expected quaternary change).

Added: a haemoglobin worksheet with teacher notes (`docs/worksheets/`),
`CONTRIBUTING.md`, and a paper template with checked references (`paper/`).

## 2.9.3

- Archived on Zenodo, so this version has a DOI and can be cited.
- Citation file now includes a summary of what Prospector does.
- No changes to the plugin's behaviour.

## 2.9.2 - first public release

What Prospector does at this version:

- **Find and load by name**: `find` / `grab` turn a protein name, gene, UniProt ID
  or PDB ID into the best structure, rank the alternatives and say why. Ambiguous
  names give a numbered list; type the number (`1`, `2` ...) to load one.
- **Understand**: `about`, `tour`, `explain`, `features`, `pocket`, `sequence`.
- **Look**: `viewas`, `look`, `scout`, `choose` (selections in plain words, with
  the PyMOL syntax shown).
- **Experiment**: `whatif` / `mutate` (before/after side by side), `measure`,
  `compare`.
- **Build**: DNA, RNA, peptides, helices, sheets, biological assemblies,
  membranes and micelles (`build`), clearly marked as idealised models or PDB data.
- **Plain sentences**: `ask` understands many requests by rules and hands only
  the unclear ones to an optional AI model (your own Anthropic API key), which
  then remembers the answer.
- **Panel**: a docked Qt panel with tabs for every feature; picking an option
  runs it at once.
- **Safe to use**: network work runs in the background with timeouts, answers are
  cached in `~/.prospector`, and everything works offline from the cache.
