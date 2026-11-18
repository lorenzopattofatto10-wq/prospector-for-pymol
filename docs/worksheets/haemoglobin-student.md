# Haemoglobin: from a name to a disease

**A one-hour lab with PyMOL and Prospector for PyMOL**

You will load the real 3D structure of human haemoglobin, find where oxygen
binds, see how one changed amino acid causes sickle-cell anaemia, and compare
haemoglobin without oxygen and with it.

**You need:** PyMOL (open-source or Incentive) with Prospector 2.9.4 or later installed
([download](https://github.com/lorenzopattofatto10-wq/prospector-for-pymol/releases/latest)),
and an internet connection.

**How it works:** type each command after `PyMOL>` in the command line and
press Enter. Prospector prints what it did and what you are looking at. Write
your answers in the spaces. Rotate with the left mouse button, zoom with the
scroll wheel.

> Tip: if something goes wrong, `back` undoes the last view and `focus all`
> shows the whole structure again.

---

## Part 1 – Meet the molecule (10 min)

```
grab haemoglobin
about haemoglobin
```

Prospector picks the best experimental structure of haemoglobin and loads it.

1. Which PDB entry was loaded? How was the structure determined, and at what
   resolution?

   ______________________________________________________________

2. How many chains does it have? Which chains are *alpha* and which are *beta*?

   ______________________________________________________________

```
look chains
```

3. Describe the shape: how are the four chains arranged? Is this haemoglobin
   with or without oxygen bound (look at the title from step 1)?

   ______________________________________________________________

---

## Part 2 – Where does oxygen bind? (15 min)

```
ask where does oxygen bind
```

4. Oxygen is not bound in this structure. Which molecule does it bind to
   instead, and which metal atom sits at its centre?

   ______________________________________________________________

5. Find the list of POCKET RESIDUES. Which **histidines** are in the pocket of
   chain A?

   ______________________________________________________________

Now measure how far two of those histidines are from the iron:

```
measure A/87 to the iron
measure A/58 to the iron
```

6. Fill in the table:

   | Residue | Distance to the iron (Å) | Bonded to the iron? |
   |---|---|---|
   | His87 (chain A) | | |
   | His58 (chain A) | | |

7. One histidine holds the heme in place; the other is on the opposite side of
   the iron, with a gap between them. What do you think goes into that gap
   when haemoglobin picks up oxygen?

   ______________________________________________________________

```
explain A/87
```

8. What does UniProt say about this residue (see ANNOTATED)? Why is histidine a
   good choice for holding a metal (see AMINO ACID)?

   ______________________________________________________________

---

## Part 3 – Sickle-cell anaemia: one letter changed (15 min)

Sickle-cell anaemia is caused by a single change in the **beta** chain.

```
explain B/6
```

9. Which amino acid is at position 6 of chain B? Is it charged? Is it buried
   or on the surface?

   ______________________________________________________________

10. Look at the line IN UNIPROT. Why does UniProt call this position 7 when the
    structure file calls it 6? (Hint: the first amino acid made by the ribosome
    is removed from the finished protein.)

    ______________________________________________________________

11. Which disease variants are known at this position? Which one is
    haemoglobin S (sickle-cell)?

    ______________________________________________________________

Now see what the change does:

```
whatif haemoglobin B/6 V
```

The original is on the LEFT, the changed protein on the RIGHT.

12. Fill in from the table Prospector printed:

    | | Before | After |
    |---|---|---|
    | Amino acid | | |
    | Charge | | |
    | Loves or avoids water? | | |

13. The change puts a water-avoiding (oily) side chain on the **surface** of the
    protein. Oily patches stick to other oily patches. Use this to explain why
    haemoglobin S molecules stick together into long fibres when they have no
    oxygen bound, and why this changes the shape of red blood cells.

    ______________________________________________________________

    ______________________________________________________________

14. Read the CAUTION line. What does this model *not* show?

    ______________________________________________________________

```
whatif off
```

---

## Part 4 – Without oxygen vs with oxygen (15 min)

The structure you loaded has no oxygen (the "T" or tense state). Now build the
oxygen-bound form (the "R" or relaxed state) and compare them.

```
build assembly 1HHO
compare haemoglobin with 1HHO_assembly
```

`build assembly` rebuilds the whole four-chain molecule: the file 1HHO stores
only half of it.

15. What is the overall RMSD between the two forms? Which residues moved the
    most (BIGGEST DIFFERENCES)?

    ______________________________________________________________

16. Try comparing just **one** pair of chains (alpha + beta):

    ```
    compare haemoglobin with 1HHO
    ```

    (1HHO on its own contains only one alpha and one beta chain.) Is the RMSD
    larger or smaller than in question 15? What does that tell you about
    *where* the change between the two forms happens – inside each chain, or in
    how the chains sit against each other?

    ______________________________________________________________

    ______________________________________________________________

---

## Part 5 – Keep your work (5 min)

```
snap haemoglobin_figure
```

saves a picture of the current view. Bonus: run `prospector teach on` and repeat
one command from this lab – Prospector now also shows the plain PyMOL commands
it used, so you can learn to write them yourself.

---

## Going further

- `tour haemoglobin` – a guided walk through the structure
- `features haemoglobin` – colours the known sites and disease variants
- `scout haemoglobin histidines` – every histidine in the molecule
- `explain B/92` – the iron-holding histidine of the beta chain
