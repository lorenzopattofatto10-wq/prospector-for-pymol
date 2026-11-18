# Prospector - Tier 4 test script (teach mode / tour / quiz / surprise).
# Paste into the PyMOL command line. Start with an empty session (File > Reinitialize).

# --- tour: five stops on thrombin; watch the view change at each "next"
grab 1PPB
tour 1PPB
next
next
next
next
next
prev

# --- teach mode: the same actions now print their PyMOL commands
prospector teach on
pocket 1PPB
look chains 1PPB
find hemoglobin limit=3
tour 1PPB
next
prospector teach off

# --- quiz: a mystery structure (different each time)
quiz
quiz hint
quiz guess banana
quiz hint
quiz reveal

# --- surprise: a random classic with its info card
surprise

# --- must fail politely
tour nothere
quiz hint
