# Prospector - Tier 1 test script.
# Paste everything into the PyMOL command line (or: File > Run Script...).
# Each block prints between two dashed lines; nothing here should freeze PyMOL.

prospector help

# 1. classic human entry near the top, with the reasoning ("Why #1")
find hemoglobin

# 2. p53 -> gene TP53 -> UniProt P04637
find p53

# 3. ambiguous name: choices, not a guess
find insulin

# 4. organism read from the text + cryo-EM filter
find spike protein sars-cov-2 method=em

# 5. typo -> "did you mean 'lysozyme'?"
find lysozime

# 6. imatinib (ligand code STI) bound kinases
find kinase ligand=STI

# 7. no experimental structure -> AlphaFold fallback
find MYMK organism=human

# 8. nonsense -> clean "nothing found"
find asdfqwerty

# extras: British spelling, a filter that removes everything, a mistyped filter
find haemoglobin organism=horse limit=3
find p53 res<0.5
find p53 organsim=human

# grab = find + load.  Loads 1 structure, 1 AlphaFold model, and refuses to guess.
grab hemoglobin
grab MYMK organism=human
grab insulin

# PyMOL's own "get" must still work (Prospector uses "grab" instead)
get cartoon_transparency

# the panel (also in the menu: Plugin > Prospector)
prospector panel
