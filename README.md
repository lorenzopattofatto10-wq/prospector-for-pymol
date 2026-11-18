# Prospector for PyMOL

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24355143.svg)](https://doi.org/10.5281/zenodo.24355143)

A PyMOL plugin for working with protein structures - and small molecules - by name. You type what you
want in plain words or short commands; Prospector finds the structure, shows
it, and explains what is on screen and why.

```
PyMOL> grab haemoglobin
PyMOL> about haemoglobin
PyMOL> ask where does the heme bind
```

**Version 2.9.** It is built for students who know some biology but are new to
PyMOL. It uses fixed rules and public databases: no account or key is needed
and the same input always gives the same answer. An AI model can optionally be
added for the sentences the rules do not understand (see below).

## Install

1. In PyMOL: **Plugin → Plugin Manager → Installed Plugins**, and uninstall any
   older `prospector`.
2. **Install New Plugin → Choose file…** and select `dist/prospector.zip`.
3. Close PyMOL completely and open it again.
4. The log should show `Prospector 2.9.4 loaded.` If it shows an older number,
   the old copy is still running.

Requirements: PyMOL 2.x or 3.x with Python 3.7+ (open-source or Incentive) and
an internet connection for the first search of each protein. No extra Python
packages are needed. If you see a certificate error on macOS, installing the
optional `certifi` package into PyMOL's Python fixes it.

To rebuild the zip after editing the code: `python build_zip.py`.

## Start here

| Command | What it does |
|---|---|
| `start` | First steps for a new user, and how to use the mouse. |
| `ask <anything>` | Say it in your own words. Prospector shows which command that is and runs it: `ask where does the drug bind in p53`, `ask load haemoglobin and show it as sticks`. Begin with `how do I` to see the command without running it. When it is unsure it offers choices and does nothing. |
| `status` | What is loaded, how it is drawn, which names you can use, what you could do next. |
| `glossary <word>` | One or two sentences on a term: RMSD, resolution, pLDDT, chain, any amino acid… |
| `info`, `info <command>`, `info filters` | The short command list, or details and examples for one command. |

## Names

Every command that works on a structure accepts a protein name (`p53`), a PDB
ID (`4HHB`), or nothing at all for the structure you worked on last.

- After `grab p53`, the name `p53` stands for what it loaded.
- A name that is not loaded yet is searched and loaded first:
  `tour haemoglobin` works on an empty session.
- Ambiguous names (`insulin`) and typos (`lysozime`) are explained, never
  guessed.

## Type it your way

Commands are forgiving about how you write them. These all do the same thing:

| You can type | Read as |
|---|---|
| `explain 248`, `explain R248`, `explain arginine 248`, `explain residue 248` | `explain 248` |
| `measure 176 238`, `measure 176 and 238`, `measure cys 176 to the zinc` | `measure 176 238` / `176 ZN` |
| `whatif R248Q`, `whatif 248 Q`, `whatif arg248gln`, `whatif R248 to Q` | `whatif R248Q` |
| `look rainbow p53`, `look p53 rainbow`, `look colour by chain` | the colour scheme, either order |
| `scout cysteines in p53`, `scout p53 cysteines` | `scout p53 cysteines` |
| `compare p53 with 1TUP`, `compare p53, 1TUP` | `compare p53, 1TUP` |

A word that cannot be a protein's name - a keyword, a number, a near-miss of
a keyword - is never looked up and loaded. You get "did you mean" instead.

## Commands

### Get a structure

| Command | What it does |
|---|---|
| `grab <name> [filters]` | Finds the best structure and loads it. |
| `find <name> [filters]` | Lists the candidates with the reasoning; loads nothing. |
| `alphafold <name>` | Loads the AlphaFold predicted model, coloured by confidence. Used automatically when no experimental structure exists. |
| `molecule <name>` | Loads a small molecule on its own (steroid, lipid, sugar, nucleotide, drug, cofactor): `molecule cholesterol`, `molecule ATP`. `grab` does the same when the name belongs to a molecule and to no protein. |
| `batch <file.txt>` | Loads a list of names or IDs, one per line. |
| `surprise` | Loads a random classic structure with its info card. |

Filters go after the name, separated by spaces:

| Filter | Example | Meaning |
|---|---|---|
| `organism=` | `organism=human`, `organism="e. coli"` | source organism |
| `method=` | `method=xray`, `em`, `nmr` | experimental method |
| `res<` | `res<2.0` | resolution better than 2.0 Å |
| `ligand=` | `ligand=ATP` | contains this bound molecule (PDB 3-letter code) |
| `year>` / `year<` | `year>2015` | release year |
| `limit=` | `limit=5` | rows shown by `find` (default 10) |
| `pick=` | `grab thrombin pick=2` | take row 2 of the list instead of the top one |

**Numbered lists:** after any numbered list (`grab myosin` when the name could
mean several proteins, or the ranked table of `find hemoglobin`), just type the
number - `1`, `2`, `3` ... - to load that line.
| `mode=text` | `find ribosome mode=text` | search the text of all entries |
| `kind=` | `find cholesterol kind=protein`, `kind=molecule` | insist on proteins, or on the small molecule |

### Small molecules

`grab cholesterol`, `grab caffeine`, `grab ATP`, `grab testosterone` load the
molecule itself as balls and sticks and say what it is. `about` prints its
card (formula, mass, charge, chirality, what it is used for). `viewas`,
`measure`, `choose`, `snap` and `back` work on it; commands made for protein
chains (`tour`, `whatif`, `sequence`, `scout` ...) say so instead of guessing.

- The first source is the PDB's own dictionary of molecules: about 45,000
  molecules that were seen bound in some structure, each with an idealised 3D
  shape.
- The second source is PubChem (NIH), for molecules the PDB does not have:
  `molecule limonene`, `molecule geosmin`. Its shapes are computed. PubChem has
  no 3D shape for very large or floppy molecules and for metal complexes;
  Prospector then says so. `grab <name>` asks PubChem only when no protein and
  no PDB molecule goes by that name; `molecule <name>` always may.
- A protein keeps its name: `grab insulin` and `grab p53` are still proteins.
  The molecule is chosen only when the name is one of its own names and no
  single protein owns it. `kind=protein` or `kind=molecule` overrides this.
- `grab protein ligand=CLR` loads a protein that has the molecule bound; the
  code (CLR) is printed whenever a molecule is found.
- Membranes and bilayers are not single molecules and cannot be loaded.

### Build a model

`build` makes structures that are not simply downloaded. Every result says
plainly what it is: an **IDEALISED MODEL** (built from geometry rules, nothing
measured) or **FROM EXPERIMENT** (real PDB data, completed with its symmetry).

| Command | What it builds |
|---|---|
| `build dna ATGCGTACGTTA [form=A]` | A DNA double helix of your sequence, B-form (as in cells) or A-form, with the partner strand added and paired. `single` for one strand. |
| `build rna GGCUAGCUAC` | One RNA strand, twisted as in an A-form helix. PyMOL cannot build double-stranded RNA; `build dna ... form=A` has nearly its shape. |
| `build helix AKLEELLKKLAE` | A peptide as an alpha-helix, with its i to i+4 hydrogen bonds drawn. |
| `build strand KVTVKVE`, `build peptide ACDEFGH` | An extended beta-strand; a fully stretched chain. |
| `build sheet KVTVKVE TVKVTVE [...]` | Two to six strands laid antiparallel side by side, slid to the position with the most hydrogen bonds. |
| `build membrane [popc|dopc|dope] [size nm] [cholesterol=30]` | A lipid bilayer patch, 2 to 20 nm wide, optionally with cholesterol. |
| `build micelle [sds|ldao] [number]` | A ball of detergent molecules, heads out and tails in. |
| `build assembly <name or PDB ID>` | The complete biological unit - e.g. the 24-chain ferritin cage from a one-chain entry. Assemblies over 300,000 atoms (virus capsids) are shown coarse: one ball per few amino acids, built from the PDB's symmetry operators, never downloaded whole. |

What to know about the models:

- **Membranes and micelles are built from straightened lipids.** Each lipid's
  tails are set to the all-trans zigzag (a cis double bond keeps its kink) and
  the shapes are packed without overlaps. The bilayer is therefore thicker
  (about 50 A between the phosphates) than a living cell's fluid membrane
  (38-40 A), and a micelle's centre is emptier than a real one. Both are said on
  screen.
- **Helices and sheets have their shape imposed.** Whether a real sequence
  would fold that way is not predicted.
- The shapes of the lipids and detergents ship with the plugin
  (`prospector/core/data/shapes`); how they were made is described in
  `README_prep.py.txt` there.
- In `ask`: `ask build a dna double helix of ATGCGTACGT`, `ask make a membrane
  with 30% cholesterol`, `ask show the whole ferritin cage`.

### Understand it

| Command | What it does |
|---|---|
| `about <name>` | Info card (protein, function, method, resolution, chains, ligands, citation) and the list of key places. Loads the structure first if it is not loaded. |
| `tour <name>` | Guided walk: fold → domains → functional sites → ligand → disease variants. `next`, `prev`, `tour 3`. |
| `focus <place>` | Zooms to a key place and says what is on screen: `focus 1`, `focus active site`, `focus heme`, `focus chain B`, `focus H/57`, `focus all`. `focus` alone lists the places. |
| `explain <A/57>` | One residue: the amino acid and what it is good for, helix/strand/loop, buried or exposed, what it touches, annotations and disease variants. |
| `sequence <name> [chain=A]` | The sequence with helices, strands and features marked underneath. |
| `features <name> [types=…]` | Colours UniProt annotations and makes a selection for each (`active`, `binding`, `domain`, `disulfide`, `ptm`, `variant`, and `site` on request). Residue numbers are converted with SIFTS and checked against the atoms. |
| `pocket <name> [ligand] [cutoff=4.0]` | The bound molecule, the residues around it and its possible hydrogen bonds. |
| `choose <description> [as <name>]` | Selects atoms from plain words and prints PyMOL's own selection line with every piece explained: `choose histidines within 5 of the zinc`, `choose charged residues in chain A`, `choose residues 100 to 120 except glycines`. `choose` alone lists the words. |
| `scout <name> <what>` | Spots things: any amino acid, `hydrophobic`, `charged`, `aromatic`, `helices`, `strands`, `disulfides`, `salt bridges`, `surface`, `buried`, `interface`, `flexible`, `ligands`, `metals`, `sugars`, `dna`, `modified`, `glycosylation sites`, `motif CxxC`. `scout` alone lists everything. |

### Experiment

| Command | What it does |
|---|---|
| `whatif <name> <R248Q> [more…]` | Changes up to five residues on a copy and shows original and changed version side by side in two panels, explaining charge, size, hydrogen bonds, contacts, DNA contact, clashes and known disease variants. `whatif off` leaves the two-panel view. |
| `mutate <name> <R248Q>` | The same change in a single view. |
| `measure <A> <B>` | Distance between two residues (`H/57`), atoms (`H/57/NE2`) or ligands (`HEM`), drawn and explained. |
| `compare <a> <b> [method=…]` | Superposes two structures, reports the RMSD and colours both by how far each residue moved. Add `:A` for one chain. |
| `similar <name> [identity=30]` | Related proteins by sequence, one row per protein. |

A change is written as original letter, position, new letter: `R248Q` means
arginine 248 becomes glutamine. Put a chain in front if needed: `B/E6V`.
Shorter forms are completed and the completion is stated: `248W` (original
read from the structure), `248`, `first`, `last` (become alanine).
`whatif` and `mutate` are sketches: only the side chain is swapped, and the
rest of the protein does not adjust.

### Change the picture

| Command | What it does |
|---|---|
| `viewas <style>` | How it is drawn: `lines`, `sticks`, `ballstick`, `spheres`, `cartoon`, `ribbon`, `tube`, `surface`, `mesh`, `dots`. Combine with `+`. |
| `look <preset>` | Colour schemes: `publication`, `chains`, `rainbow`, `secondary` (helix/strand/loop), `types` (kind of amino acid), `elements`, `hydrophobicity`, `putty`, `ligand`. |
| `back`, `back 3` | Goes back one step, or several (up to 10): restores the picture and removes what the step added (a loaded structure, a measurement, a changed copy, a selection). |
| `clear` | Empties the 3D window. `back` right after it restores everything. `clear text` wipes only the text output. |
| `snap <file> [width=2400] [dpi=300]` | Saves a ray-traced picture; never overwrites. |

### Keep your work

| Command | What it does |
|---|---|
| `session save <name>`, `session open <name>`, `session` | The whole PyMOL session and the names you used. |
| `report <name>` | One HTML page with the info card, key places and three pictures. |
| `history`, `fav`, `fav add`, `fav rm` | Past searches and favourites. |

### Practise

| Command | What it does |
|---|---|
| `quiz`, `quiz hint`, `quiz guess <name>`, `quiz reveal` | A mystery structure with three questions. |
| `prospector teach on` / `off` | Prints the PyMOL commands behind every action, with comments. |
| `prospector panel` | The same things as buttons (also under **Plugin > Prospector for PyMOL**). It opens docked inside PyMOL's window, under the 3D view, so it never covers the structure; drag its title bar to move it to another edge or pull it out as its own window. An Ask box on top, the structure you are working on under it, then five tabs in working order: 1 Find, 2 Understand, 3 Look, 4 Experiment, 5 Keep and learn. Back, Clear, More and the command that was run stay visible at the bottom. With nothing loaded, the buttons load the name typed in tab 1. Click a residue in the 3D view and the **Clicked** row names it and offers Explain it, Zoom to it, Change it; the boxes in tabs 2 and 4 are filled in as well. Measuring is deliberate: press **Measure a distance...**, then click two points, and the distance between them is measured at once; other clicks never start a measurement. A box on the right, **What just happened**, shows the result of every button in plain words. **Bigger panel / Smaller panel** switches its size. In tab 4 the new amino acid is picked from a list, so no letter codes are needed. |
| `prospector short on` / `off` | Short mode: action commands print only the answer (notes are always kept). `more` prints the full text of the last command. |

## Optional: an AI for the sentences the rules miss

`ask` first uses its built-in rules: instant, free, offline, the same answer
every time. If you store an API key, a sentence the rules do not understand
is sent to an AI model (Claude, by Anthropic), which picks Prospector commands
for it. Without a key nothing changes and nothing is sent.

```
PyMOL> prospector ai key        (opens a box; paste the key there)
PyMOL> prospector ai test       (one small question to check it works)
PyMOL> ask fammi vedere dove si lega il farmaco e colora il resto in grigio
```

| Command | What it does |
|---|---|
| `askai <anything>` | Like `ask`, but skips the rules and goes straight to the AI. For other languages, follow-ups, or when `ask` picked the wrong command. |
| `prospector ai` | Status, what is sent, what it costs. |
| `prospector ai key` | Store your key. Create one at console.anthropic.com. Never type the key on the command line: PyMOL logs it. |
| `prospector ai on` / `off` | Use the AI, or never use it. |
| `prospector ai test` | Check that the key and the connection work. |
| `prospector ai model <name>` | Default `claude-opus-5-5`; `claude-haiku-5-5` is much cheaper. |
| `prospector ai forget` | Delete the stored key. |
| `prospector ai learned` / `unlearn` | List, or forget, the sentences `ask` learned from the AI. |

- **What is sent:** your sentence, the list of commands, and the names of what
  is loaded (structures, chains, bound molecules, key places, your last few
  commands). No files and no structures.
- **What it costs:** about 1 cent per question with the default model, on your
  own Anthropic account. Questions the rules answer cost nothing.
- **What it can do:** only propose command lines. Prospector runs them if every
  one is a Prospector command or a simple display command of PyMOL (`hide`,
  `show`, `color`, `set` for a fixed list of settings, `distance`, `zoom`,
  `orient`, `center`, `turn`, `rock`, `bg_color`). If one line is anything
  else, none of them run. The lines are always shown, marked `(AI)`.
- **It learns.** When the AI answers a sentence the rules did not understand
  and the commands run, the pair is saved in `~/.prospector/learned.json`. The
  next time, `ask` runs that sentence at once, with no request and no cost,
  marked `(learned)`. If the answer named the structure on screen, it is saved
  as "the structure being worked on" and follows whatever is loaded next time.
  Sentences that point back ("the same for chain B") are not saved. `askai` always asks afresh.
- **Where the key is:** `~/.prospector/ai_key`, readable by you only. The
  environment variable `ANTHROPIC_API_KEY` is used if no key is stored.
- The AI can misread you, and the same sentence can give a different answer on
  another day. `back` and `clear` undo what it did to the picture.

## Names that differ from PyMOL's

Some obvious names belong to PyMOL itself, so Prospector uses others and leaves
PyMOL's alone:

| You might expect | Prospector uses | PyMOL's own command does |
|---|---|---|
| `get` | `grab` | reads a setting |
| `view` | `viewas` | saves and recalls camera positions |
| `undo` | `back` | undoes atom edits |
| `save` | `session save` | writes a file |

Every Prospector command also exists with a `p_` prefix (`p_find`, …).

## How the ranking works

Each structure gets up to 100 points:

| Component | Points | Full marks for |
|---|---|---|
| Relevance | 30 | the right protein, named in the title, no engineered mutations, on its own |
| Coverage | 20 | the whole UniProt sequence actually modelled |
| Resolution | 20 | 1.0 Å or better (4.0 Å = 0; NMR gets half) |
| Organism | 15 | your organism filter; otherwise human, then model organisms |
| Method | 10 | X-ray; cryo-EM 9; NMR 6 |
| Recency | 5 | released this year |

Structures on a hand-curated list of about 80 textbook entries (marked `*`) get
8 bonus points. When a name is a whole class (`kinase`, `ribosome`) the
textbook structures of that class are listed first. The weights live in
`prospector/core/rank.py`.

## Limits worth knowing

- In `ask`, "this residue" or "the selected residue" means the one you last
  clicked in the 3D window: `ask why is this residue important`.
- `ask` only runs a command when it can account for every meaningful word of
  your sentence. If a word is left over, the sentence is in another language,
  or a "name" does not look like a name, it does not guess: the sentence goes
  to the AI when a key is stored, and is refused with examples otherwise.
- Without the optional AI, `ask` follows fixed rules. It copes with plurals, common
  typos and many ways of saying the same thing ("catalytic residues" =
  "active site", "sulfur bridges" = "disulfides"), but a sentence about
  something Prospector has no command for is refused, not guessed.
  `tests/ask_phrases.txt` lists about 380 sentences it is tested on; add a
  line there when you find one it should understand.
- `back` cannot bring back something you deleted yourself with PyMOL's
  `delete`; only `clear` keeps a copy.
- `choose` knows a fixed vocabulary. It never guesses at an unknown word; it
  names the word and selects nothing.
- For a few simple requests `ask` runs one of PyMOL's own commands (`hide`,
  `color`, `rock`, `set`, `distance`) and shows it to you.
- A name picks the best-ranked structure, which may not be the one in your
  textbook. `find <name>` shows the list; type a row's number (or use its PDB ID) to get another.
- `whatif`, `mutate` and `compare` verdicts are aids for looking, not
  predictions.
- `snap`, `report` and the surface styles freeze PyMOL for a few seconds while
  they draw.

## Files Prospector writes

Everything is under `~/.prospector` (set `PROSPECTOR_HOME` to move it):

- `cache` — answers from the databases (searches 1 day, metadata 7 days; used
  beyond that when offline)
- `structures` — downloaded structure files
- `sessions` — saved sessions
- `history.json`, `favorites.json`, `settings.json`
- `ai_key` — your AI key, if you stored one

## For developers

```
prospector/core/         pure Python: API clients, resolver, ranking, cache
prospector/pymol_layer/  PyMOL commands, worker thread, Qt panel
tests/                   offline tests on recorded API answers
```

```bash
python -m unittest discover -s tests -v
```

```bash
PROSPECTOR_LIVE=1 python -m unittest tests.test_live -v
```

```bash
PROSPECTOR_QT_TEST=1 QT_QPA_PLATFORM=offscreen python -m unittest tests.test_pymol -v
```

```bash
python -m tests.record_fixtures
```

The first runs the offline suite, the second the tests that use the network,
the third the background-thread test, and the last re-records the API answers.

The AI client is `prospector/core/ai.py` (standard library only, one request
per question). Its tests use a fake service, so they run without a key.

`tests/manual/` holds PyMOL scripts (`.pml`) for checking things by eye.

To build the file for PyMOL's Plugin Manager:

```bash
python build_zip.py
```

It writes `dist/prospector.zip` (with the licence and this README inside).

## Teaching materials

- [Haemoglobin: from a name to a disease](docs/worksheets/haemoglobin-student.md) -
  a one-hour lab for first-year students: the heme pocket, sickle-cell anaemia
  (Glu6Val) and the change between haemoglobin with and without oxygen.
  [Teacher notes and answer key](docs/worksheets/haemoglobin-teacher.md).
  Printable PDFs: [student sheet](docs/worksheets/haemoglobin-student.pdf),
  [teacher notes](docs/worksheets/haemoglobin-teacher.pdf).

Want to help, report a problem or add a worksheet? See [CONTRIBUTING.md](CONTRIBUTING.md).

## Data sources and credits

Prospector ships no structures itself. It asks these free public services, and
what it shows comes from them:

| Service | Used for | Terms |
|---|---|---|
| [RCSB PDB](https://www.rcsb.org) (Search API, Data API, files, Chemical Component Dictionary) | structures, assemblies, ligands | PDB data is free of all restrictions (CC0 1.0) |
| [UniProt](https://www.uniprot.org) | protein names, genes, functions, sites, variants | CC BY 4.0 - credit UniProt |
| [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk) | predicted models | CC BY 4.0 - credit AlphaFold DB (Jumper et al. 2021; Varadi et al. 2024) |
| [PDBe SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) | which PDB chains match which UniProt protein | EMBL-EBI terms of use |
| [PubChem](https://pubchem.ncbi.nlm.nih.gov) (PUG REST) | small molecules not in the PDB | free to use; cite PubChem (Kim et al., Nucleic Acids Res.) and respect its request limits |
| [Anthropic API](https://www.anthropic.com) (optional) | plain-sentence requests the rules cannot read | your own key and account; nothing is sent without one |

The lipid and detergent shapes used by `build membrane` and `build micelle`
(`prospector/core/data/shapes/`) were prepared from the ideal coordinates in the
PDB's Chemical Component Dictionary.

If you publish work that used structures found with Prospector, cite the
structures themselves (their PDB IDs and papers) and the databases above.

PyMOL is a trademark of Schrödinger, LLC. Prospector for PyMOL is an independent plugin
and is not made or endorsed by Schrödinger.

## Licence

MIT - see [LICENSE](LICENSE). You may use, change and share Prospector for PyMOL freely,
including in teaching and commercial work, as long as the licence notice stays
with it. To cite Prospector, use the DOI [10.5281/zenodo.24355143](https://doi.org/10.5281/zenodo.24355143), for example:

> Pattofatto, L. (2026). *Prospector for PyMOL* (Version 2.9.3) [Computer software].
> Zenodo. https://doi.org/10.5281/zenodo.24355143

GitHub's **Cite this repository** button (from [CITATION.cff](CITATION.cff)) gives the
same in other formats, such as BibTeX.
