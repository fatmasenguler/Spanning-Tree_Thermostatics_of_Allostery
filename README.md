# Spanning-Tree Thermostatistics of Protein Allostery

**An Exact Kirchhoff Framework with Application to Oncogenic KRAS**

**Fatma Senguler Ciftci and Burak Erman**
Department of Chemical and Biological Engineering, Koç University, Istanbul, Turkey

Published in *Physical Biology* (2026), gold open access under CC BY 4.0.

[![DOI](https://img.shields.io/badge/DOI-10.1088%2F1478--3975%2Fae9e79-blue)](https://doi.org/10.1088/1478-3975/ae9e79)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey)](https://creativecommons.org/licenses/by/4.0/)

---

## Overview

This repository contains the computational notebooks accompanying the paper:

> **Spanning-Tree Thermostatistics of Protein Allostery: An Exact Kirchhoff Framework with Application to Oncogenic KRAS**

The framework treats protein allostery as a statistical-mechanical problem on the spanning-tree ensemble of a residue contact network. Residues are represented by Cα atoms, contacts define the edges of a weighted graph, and each spanning tree is a cycle-free topological microstate of the global communication network.

For a given weighted contact graph, the global spanning-tree partition function is evaluated exactly by Kirchhoff's Matrix-Tree Theorem,

$$Z = \det \widetilde{L}_w,$$

where $\widetilde{L}_w$ is any reduced weighted Kirchhoff (graph Laplacian) matrix. No stochastic sampling of the spanning-tree ensemble is required.

Allosteric channels are defined through simple paths between a source and a target residue. The probability that a specified path occurs in the spanning-tree ensemble is evaluated exactly by the Burton–Pemantle theorem and the Moore–Penrose pseudoinverse of the weighted Laplacian.

Because the number of simple paths grows combinatorially, channel thermodynamics uses a finite path-length expansion truncated at a maximum node length $L = 9$. Convergence is assessed independently against the exact full-network dynamic distance obtained from the Kirchhoff pseudoinverse.

The framework is applied to the oncogenic G12D mutation in KRAS, comparing the active-state structures

* **WT:** PDB [6GOD](https://www.rcsb.org/structure/6GOD)
* **G12D:** PDB [6GOF](https://www.rcsb.org/structure/6GOF)

using chain A of each structure.

---

## Quick verification

Any correct implementation of the contact graph should reproduce these counts at $r_c = 7.8$ Å, chain A:

| Quantity | 6GOD (WT) | 6GOF (G12D) |
| --- | ---: | ---: |
| Residues (Cα) | 172 | 172 |
| Contacts | 796 | 803 |
| Shared contacts | 793 | 793 |
| Cycle rank $\xi = E - N + 1$ | 625 | 632 |

The G12D substitution gains 10 contacts and loses 3, a net increase of seven. Cα RMSD between the two chains is 0.19 Å (TM-score 0.998), so the two graphs are near-identical, and any differences in the results come from the ensemble statistics rather than from a conformational change.

Path counts in the truncated channel ensemble ($L \le 9$) provide a second checkpoint:

| Channel | WT paths | G12D paths |
| --- | ---: | ---: |
| pre-Switch II, 55 → 60 | 638,606 | 642,205 |
| inter-lobe linker, 19 → 142 | 864,439 | 903,899 |
| canonical axis, 12 → 61 | 124,222 | 134,205 |

---

## Main findings

### Global thermodynamics is strongly conserved

At the reference parameters $r_c = 7.8$ Å and $kT = 1.0$ Å:

| Property | WT | G12D | Change |
| --- | ---: | ---: | ---: |
| Mean energy, $\langle E \rangle$ | 783 | 784 | +0.13% |
| Entropy, $S$ | 267 | 268 | +0.37% |
| Heat capacity, $C$ | 139.6 | 140.3 | +0.50% |

The G12D mutation produces very little change in the global thermodynamic state of the spanning-tree ensemble despite substantial changes in internal routing.

### Channel capacity is preserved while internal routing changes

Channel heat capacities are also stable under mutation, $|\Delta C| \le 0.08$ in every channel, but this stability hides large, opposing shifts in the energetic, topological, and cross-coupling contributions.

The canonical 12 → 61 channel (G12 → Q61 / Switch II) is the clearest example:

$$\Delta C_E = +2.3445, \qquad \Delta C_T = +2.4503, \qquad \Delta C_X = -4.8021,$$

giving a net change of only $\Delta C = -0.0072$.

### G12D redistributes residue-level allosteric importance

The mutation changes which intermediate residues carry the channel ensemble. Reported shifts include:

| Residue | Channel | $\Delta I_k$ |
| --- | --- | ---: |
| Q61 | pre-Switch II, 55 → 60 | +22.7% |
| A155 | inter-lobe linker, 19 → 142 | +29.1% |
| F156 | inter-lobe linker, 19 → 142 | +35.5% |
| G12 | P-loop, 6 → 11 | −17.8% |
| G12 | pre-Switch II, 55 → 60 | −25.9% |

The A155/F156 region of the inter-lobe linker and the Q61 Switch II axis take up traffic that G12 gives up. The mutation does not eliminate allosteric communication; it reorganizes path occupancy among intermediate residues while preserving global and channel-level thermodynamic capacity.

The long-range 12 → 170 channel is the least converged of the ten, reaching about 87% of the exact dynamic resistance at $L = 9$ and about 95% at $L = 10$. Its residue-level shifts should be read as provisional.

---

## Model

### Cα contact graph

The protein is represented by a weighted graph $G = (V, E)$, where vertices are Cα atoms and an edge is present when

$$d_{ij} \le r_c, \qquad r_c = 7.8\ \text{Å}.$$

Edge weights are

$$w_{ij} = \exp\left(-\frac{d_{ij}}{kT}\right), \qquad kT = 1.0\ \text{Å}.$$

Backbone contacts receive no separate parameter; they are governed by the same expression as all non-bonded edges.

### Global spanning-tree partition function

For a spanning tree $T$, the statistical weight is the product of its edge weights,

$$W(T) = \prod_{e \in T} w_e,$$

and the global partition function is

$$Z = \sum_T W(T) = \det \widetilde{L}_w.$$

Free energy, internal energy, entropy, and heat capacity all follow from $Z$.

### Path probabilities

For a specified simple path $p$, the Burton–Pemantle theorem gives the probability that all edges of $p$ occur simultaneously in a weighted spanning tree. Using the Moore–Penrose pseudoinverse $L_w^+$, an edge-response matrix $K_p$ is built for the edges of $p$, and

$$P(p) = \det K_p.$$

Each specified path is therefore evaluated exactly from the global Kirchhoff matrix rather than estimated by sampling random spanning trees.

### Channel heat-capacity decomposition

$$C_{AB} = C_E + C_T + C_X,$$

where $C_E$ is the energetic contribution from physical path distances, $C_T$ the topological contribution from spanning-tree degeneracy, and $C_X$ the energetic–topological cross-coupling. Mutation-induced changes satisfy $\Delta C = \Delta C_E + \Delta C_T + \Delta C_X$.

### Participation ratio

For normalized path weights $p_i$,

$$PR = \frac{1}{\sum_i p_i^2},$$

the effective number of paths participating in a channel. A small $PR$ indicates localization onto a few dominant routes.

### Residue-level allosteric importance

For residue $k$, the allosteric importance $I_k$ is the fraction of the channel path ensemble passing through that residue, and mutation-induced changes are reported as

$$\Delta I_k = 100 \times \frac{I_k^{\mathrm{G12D}} - I_k^{\mathrm{WT}}}{I_k^{\mathrm{WT}}}.$$

Positive values indicate recruitment into the mutant channel ensemble, negative values reduced participation. Only residues with a wild-type occupancy $I_k^{\mathrm{WT}} \ge 3\%$ are reported. $I_k$ is channel-specific, so a residue may shift in opposite directions in different channels.

---

## Allosteric channels analyzed

| # | Channel | Residues | Structural role |
| --- | --- | --- | --- |
| 1 | P-loop | 6 → 11 | Local phosphate-binding region |
| 2 | pre-Switch II | 55 → 60 | Pre-Switch II hinge |
| 3 | GBS | 110 → 117 | Guanine-base-binding region |
| 4 | SAK motif | 141 → 146 | Distal nucleotide-binding interface |
| 5 | Lobe linker | 19 → 142 | Inter-lobe bridge |
| 6 | G12 → Switch I | 12 → 35 | P-loop to Switch I |
| 7 | G12 → Q61 / Switch II | 12 → 61 | Canonical G12–Q61 signaling axis |
| 8 | G12 → α5 | 12 → 156 | P-loop to C-terminal α5 region |
| 9 | G12 → C-terminus | 12 → 170 | Long-range C-terminal corridor |
| 10 | Switch I → Switch II | 35 → 61 | Direct inter-switch coupling |

---

## Repository contents

```text
Spanning-Tree_Thermostatistics_of_Allostery/
│
├── 1_channel_convergence_analysis.ipynb
├── 2_compute_cv_table_2.ipynb
├── 3_cv_decomposition.ipynb
├── 4_participation_ratio.ipynb
├── 5_jsd_rewiring.ipynb
├── LICENSE
└── README.md
```

**`1_channel_convergence_analysis.ipynb`** — convergence of the truncated channel representation against the full-network Kirchhoff dynamic distance (Figure 2 of the paper).

**`2_compute_cv_table_2.ipynb`** — WT and G12D channel heat capacities (Table 2).

**`3_cv_decomposition.ipynb`** — the energetic, topological, and cross-coupling contributions $\Delta C_E$, $\Delta C_T$, $\Delta C_X$ (Table 2).

**`4_participation_ratio.ipynb`** — channel participation ratios.

**`5_jsd_rewiring.ipynb`** — Jensen–Shannon divergence and residue-level allosteric-importance changes (Table 3), computed from WT and G12D path-ensemble files.

---

## Reproducibility status

Two of the three layers reproduce from structure alone:

* the **global partition function** is evaluated exactly by Kirchhoff's Matrix-Tree Theorem, directly from the PDB files;
* the probability of each **specified path** is evaluated exactly by the Burton–Pemantle theorem;
* the complete **channel path ensemble** is a controlled finite path-length expansion, tested explicitly for convergence.

**Known gap.** `5_jsd_rewiring.ipynb` reads pre-generated WT and G12D path-ensemble files; the notebook that generates those ensembles from the PDB structures is not yet in the repository. Until it is added, Table 3 cannot be reproduced end-to-end from structure. This is the highest-priority addition.

---

## Running the notebooks

The notebooks are written for Google Colab. Open a notebook and use its **Open in Colab** button; each installs its dependencies and prompts for input files.

For calculations starting from structure, download the PDB files and upload them when prompted:

```bash
wget https://files.rcsb.org/download/6GOD.pdb
wget https://files.rcsb.org/download/6GOF.pdb
```

Dependencies:

```text
numpy
scipy
networkx
biopython
pandas
matplotlib
```

---

## Reference parameters

| Parameter | Reference value |
| --- | ---: |
| Protein chain | A |
| Cα contact cutoff, $r_c$ | 7.8 Å |
| Effective temperature, $kT$ | 1.0 Å |
| Maximum channel path length, $L_{\max}$ | 9 nodes |
| WT structure | 6GOD |
| G12D structure | 6GOF |

Parameter sensitivity is examined in the paper across $r_c = 7.5$–8.5 Å and $kT = 0.8$–1.2, where the principal channel and residue-level trends remain stable.

---

## Related work

* Ciftci F S, Erman B. *Generalizing the Gaussian Network Model: Spanning-Tree Thermodynamics Shows Entropy-Driven KRAS Activation.* **Proteins: Structure, Function, and Bioinformatics** (2026). [doi.org/10.1002/prot.70146](https://doi.org/10.1002/prot.70146)

---

## Citation

```bibtex
@article{SengulerCiftci2026SpanningTree,
  author  = {Senguler Ciftci, Fatma and Erman, Burak},
  title   = {Spanning-Tree Thermostatistics of Protein Allostery:
             An Exact Kirchhoff Framework with Application to Oncogenic KRAS},
  journal = {Physical Biology},
  year    = {2026},
  doi     = {10.1088/1478-3975/ae9e79}
}
```

---

## License

Released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), matching the licence of the published article.
