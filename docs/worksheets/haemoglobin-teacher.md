# Haemoglobin lab – teacher notes and answer key

For the student sheet [haemoglobin-student.md](haemoglobin-student.md).

**Level:** first-year biology, biochemistry or biotechnology.
**Time:** about 60 minutes. Parts 1–3 can be run alone as a 40-minute session.
**Learning goals:** students can

1. load a protein structure by name and read its basic facts (method,
   resolution, chains);
2. locate a binding site and relate metal coordination to measured distances;
3. explain how a single amino-acid substitution (Glu6Val) leads to disease;
4. distinguish changes within a chain from changes between chains (quaternary
   structure, T vs R state).

**Before the class:** install Prospector on the lab computers
(Plugin → Plugin Manager → Install New Plugin → `prospector.zip`) and run
`grab haemoglobin` once on each machine: the files are then cached in
`~/.prospector`, so the class still works if the network is slow. The AI
feature is not needed for this lab.

The values below were obtained with Prospector 2.9.3. A later PDB release could
make Prospector pick a different haemoglobin entry in Part 1; if so, have
students type `grab 1A3N` so the numbers match.

---

## Part 1

1. **1A3N** – *Deoxy human hemoglobin*; X-ray diffraction, **1.80 Å**.
2. **Four** chains: **A and C = alpha** (141 aa), **B and D = beta** (146 aa).
3. A compact tetramer: two alpha–beta pairs related by symmetry, each chain
   holding one heme. The title says **deoxy**: no oxygen bound.

## Part 2

4. **Heme** (HEM, protoporphyrin IX containing Fe), with an **iron** atom at its
   centre. Prospector notes that O₂ is not bound in this entry and shows the
   heme pocket instead.
5. **His45, His58, His87** (chain A pocket residues include MET32 … HIS45 PHE46
   HIS58 … HIS87 …).
6. | Residue | Distance (Å) | Bonded? |
   |---|---|---|
   | His87 (A) | **2.24** (NE2–Fe) | **yes** – "close enough to be bonded to the metal ion" |
   | His58 (A) | **4.19** (CE1–Fe) | **no** – in contact, not bonded |

   (Beta chain equivalents: His92 2.26 Å, His63 4.22 Å.)
7. **Oxygen (O₂)** binds to the iron in the gap on the His58 side. His87 is the
   *proximal* histidine (F8), His58 the *distal* histidine (E7), which
   stabilises the bound O₂ by a hydrogen bond.
8. UniProt: *"Binding site: heme b – proximal binding residue"*. Histidine's
   imidazole ring has a nitrogen with a free electron pair that can coordinate
   metals such as iron, zinc and copper.

## Part 3

9. **Glutamate (Glu6)**, negatively charged, **on the surface** (56% of it can
   touch water), in an alpha-helix.
10. UniProt counts from the initiator **methionine**, which is removed from
    mature haemoglobin; the structure (and the traditional name "Glu6Val")
    counts the mature protein. HGVS notation is therefore p.Glu7Val.
11. E>A, E>K, E>Q, E>V. **E>V is Hb S** (sickle cell; heterozygotes are
    protected against malaria; dbSNP rs334). E>K is Hb C.
12. | | Before | After |
    |---|---|---|
    | Amino acid | glutamate | valine |
    | Charge | negative | none (a charge is lost) |
    | Water | water-loving | water-avoiding ("an oily patch on the surface") |

    Prospector also reports: one hydrogen bond lost (to Thr4), no clashes.
13. Expected reasoning: the surface Val6 of one Hb S molecule fits into a
    hydrophobic pocket on the beta chain of a neighbouring molecule (formed by
    Phe85 and Leu88). This contact exists only in the deoxy (T) shape, so when
    oxygen is released the molecules polymerise into long fibres, which distort
    red blood cells into the sickle shape, make them rigid and shorten their
    life.
14. Only the side chain was swapped, in its most likely shape; the rest of the
    protein was not allowed to adjust. It is a sketch of one molecule, not a
    prediction – and it cannot show the contact with a *second* molecule.

## Part 4

15. Whole tetramer: **RMSD 2.32 Å** ("same fold with local differences"). The
    biggest differences are at **His146 of both beta chains (8.1 and 7.7 Å)**,
    Asp94β, Lys95β and Tyr145β – the C-terminus of the beta chain.
16. One alpha–beta pair only (1A3N vs 1HHO): **RMSD 0.54 Å**, "essentially the
    same conformation". So each chain – and each alpha–beta pair – hardly
    changes; the T→R transition is mainly a **rotation of one alpha–beta pair
    relative to the other** (about 15°): a change in **quaternary** structure.
    The salt bridge between **His146β and Asp94β**, present in the T state, is
    broken in the R state; its protonation is part of the **Bohr effect** (lower
    pH stabilises the T state and promotes O₂ release in tissues).

---

## Common problems

| Problem | Fix |
|---|---|
| "Nothing is loaded yet" | run `grab haemoglobin` first |
| A different PDB ID was loaded in Part 1 | `grab 1A3N` |
| The view is cluttered | `focus all`, or `back` |
| `whatif` left two panels on screen | `whatif off` |
| No internet in the lab | run all `grab`/`build` commands once beforehand; the files are cached |
| Students want plain PyMOL | `prospector teach on` shows the PyMOL commands behind every step |

## Assessment ideas

- Explain in three sentences why a single surface mutation can cause a disease.
- Predict what `whatif haemoglobin B/6 K` (Hb C) changes, then check it.
- Use `measure` to test whether His92 of chain B is bonded to its iron.
