# Prospector - Tier 2 test script (about / alphafold / features / pocket).
# Paste into the PyMOL command line. Start with an empty session (File > Reinitialize).

# --- thrombin: the hard numbering case (offset + insertion codes like 60A) ---
grab 1PPB
about
# Expect: active site UniProt 406/462/568 -> H/57, H/102, H/195, and
# "Mapping check: 18 of 18 mapped residues ..."
features 1PPB
# Expect: ligand 0G6, pocket residues incl. HIS57 ASP189 SER195, polar contacts
pocket 1PPB

# --- hemoglobin: constant offset of 1, two proteins in one entry ---
grab 4HHB
features 4HHB types=binding,domain
pocket 4HHB HEM

# --- imatinib bound to Abl kinase ---
grab 1IEP
pocket 1IEP STI cutoff=4.5

# --- AlphaFold: by UniProt ID and by name; features work on models too ---
alphafold P04637
about AF_P04637_F1
features AF_P04637_F1 types=binding,domain

# --- things that should fail politely ---
features 1PPB types=colour
pocket 1PPB XYZ
about nothere
grab 1BNA
features 1BNA
