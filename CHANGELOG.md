# Changelog - Prospector for PyMOL

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
