# Prospector - Tier 3 test script (compare / similar / batch / look / snap / history / fav).
# Paste into the PyMOL command line. Start with an empty session (File > Reinitialize).

# --- compare: oxy vs deoxy haemoglobin (same protein -> method "align", RMSD < 1 A)
compare 1HHO 4HHB
# --- one haemoglobin chain vs myoglobin (28% identical -> method "super")
compare 4HHB:A 1MBN
# --- unrelated proteins (-> "cealign", with a caution that the fold differs)
compare 1MBN 1UBQ
# --- names instead of IDs, separated by a comma
compare human lysozyme c, hen lysozyme c
# --- must refuse politely
compare hemoglobin
compare insulin 1MBN

# --- similar: other myoglobins, one row per protein
similar 1MBN identity=20 limit=8

# --- look: list the presets, then try them
look
look chains 4HHB
look rainbow 1MBN
look hydrophobicity 1MBN
look putty 4HHB
look ligand 4HHB
look publication

# --- snap: writes prospector_test.png into PyMOL's current folder (type "pwd" to see it)
snap prospector_test width=1200
snap prospector_test width=1200

# --- history and favourites (both survive restarting PyMOL)
find p53
fav add 4HHB classic deoxyhaemoglobin
fav
history
fav rm 4HHB

# --- batch: needs a text file. Create one with a name or PDB ID per line, e.g.
#       4HHB
#       myoglobin organism=horse
#       insulin
#     then run:  batch /full/path/to/that_file.txt
batch no_such_file.txt
