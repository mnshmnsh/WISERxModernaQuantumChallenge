# RNA Secondary Structure Prediction via Quantum–Classical Hybrid Optimization

**WISER × Moderna Quantum Challenge 2026 — Team QuantumCFO**

A quantum-classical hybrid pipeline for predicting minimum free energy (MFE) RNA secondary structure, formulated as a QUBO and solved with QAOA, VQE, and a Hadamard-transform interference decoder, benchmarked against ViennaRNA and scaled to a 44-nucleotide challenge sequence via a hierarchical divide-and-conquer solver.

---

## Table of Contents

- [1. The Challenge](#1-the-challenge)
- [2. Our Approach](#2-our-approach)
- [3. Repository Structure](#3-repository-structure)
- [4. Setup and Reproduction](#4-setup-and-reproduction)
- [5. Methods and Tools](#5-methods-and-tools)
- [6. Results](#6-results)
- [7. Scaling and Resource Analysis](#7-scaling-and-resource-analysis)
- [8. Limitations and Next Steps](#8-limitations-and-next-steps)
- [9. Optional Advanced Tasks Addressed](#9-optional-advanced-tasks-addressed)
- [10. Team](#10-team)
- [11. Data and Privacy](#11-data-and-privacy)

---

## 1. The Challenge

[WISER](http://www.thewiser.org/summer-program-2026) and Moderna posed the problem of predicting RNA secondary structure — specifically the minimum free energy (MFE) fold — using quantum or quantum-inspired optimization. mRNA secondary structure affects molecular stability, translation efficiency, and manufacturability of mRNA therapeutics, and the number of possible folds grows exponentially with sequence length, making exhaustive classical search intractable at scale.

The challenge asked participants to:
1. Formulate RNA folding as an optimization problem suitable for quantum methods.
2. Benchmark candidate structures against ViennaRNA's classical MFE reference.
3. Analyze how quantum resource requirements (qubits, circuit depth) scale with sequence length.
4. Reproduce known benchmark structures for small sequences, not necessarily outperform classical tools.

## 2. Our Approach

We formulate RNA secondary structure prediction as a Quadratic Unconstrained Binary Optimization (QUBO) problem, with coefficients grounded directly in ViennaRNA's thermodynamic model rather than hand-picked weights:

- **Decision variables:** one binary variable per chemically valid candidate base pair (Watson-Crick or wobble pairing, minimum hairpin loop size of 3).
- **Diagonal terms:** real ViennaRNA loop-closure energy per candidate pair.
- **Off-diagonal terms:** a stacking-cooperativity bonus between adjacent nested pairs — required for correctness, since omitting it collapses the QUBO optimum to the unfolded structure.
- **Constraints:** one-pair-per-base and no-pseudoknot penalties, applied strictly as cross-terms (only when conflicting variables are both selected).

Three solvers are benchmarked against exact brute force and ViennaRNA's classical MFE:

| Solver | Description |
|---|---|
| **QAOA** | Quantum Approximate Optimization Algorithm, p = 2 layers |
| **VQE** | Variational Quantum Eigensolver, hardware-efficient (RY + ring-CNOT) ansatz |
| **Interference+Greedy** | Hadamard transform → phase encoding of the exact Ising cost Hamiltonian → second Hadamard transform → classical greedy bit-flip refinement |

To bridge the gap between NISQ-scale qubit budgets (~16 qubits) and the 313-qubit direct requirement of the challenge's 44-nt example sequence, we implement a **hierarchical divide-and-conquer solver**: it recursively identifies candidate anchor helices, solves flanking and interior regions independently, and reassembles the structure by minimizing total ViennaRNA energy — recovering the exact 44-nt MFE structure using at most 16 qubits per subproblem.

## 3. Repository Structure

```
.
├── README.md                          # this file
├── rna_fold_qubo.ipynb                # implementation -- open and run in Colab
└── docs/
    ├── rna_fold_qubo_writeup.pdf       # full write-up (methodology, results, limitations)
    └── rna_fold_qubo_deck.pptx         # presentation deck
```

> If you're setting this repository up from scratch: place `rna_fold_qubo.ipynb` at the repo root, and the write-up PDF and deck under `docs/`, or adjust the paths above to match wherever you place them — just keep this README's links in sync with the actual layout.

## 4. Setup and Reproduction

1. Open `rna_fold_qubo.ipynb` in [Google Colab](https://colab.research.google.com/).
2. Run all cells top to bottom (Runtime → Run all). The first cell installs dependencies (`pennylane`, `ViennaRNA`) directly, so no separate environment setup is needed.
3. All 13 numbered sections — benchmarking, qubit scaling, solver comparison, the 44-nt hierarchical solve, encoding comparison, noise robustness, and the simulator backend comparison — execute in order and print/plot their results inline.

GitHub renders `.ipynb` files natively, so the code and (once committed) its output are readable directly in the repository without running anything. Re-running the notebook top to bottom is the fastest way to verify every result in the write-up is reproducible from scratch.

No GPU or real quantum hardware is required; all quantum circuits run on PennyLane's `lightning.qubit` simulator (noise-robustness studies use `default.mixed` for density-matrix simulation).

## 5. Methods and Tools

- **Quantum SDK:** PennyLane (`lightning.qubit` simulator)
- **Classical RNA folding / reference energies:** ViennaRNA Python package
- **Classical optimization:** SciPy (COBYLA, for QAOA/VQE parameter optimization)
- **QUBO encodings:**
  - *Pair encoding* — one qubit per candidate base pair (used by the QAOA/VQE/Interference+Greedy benchmarks).
  - *Partner-index encoding* — one binary register per base, log-encoding its partner index (included as a comparison point, not used by any solver).
  - *Helix-level encoding* — groups contiguous stacked pairs into helical runs, reducing qubit count 40–75% relative to pair encoding; this is what the hierarchical solver uses to reach the 44-nt result.

## 6. Results

| Metric | Result |
|---|---|
| Solver-comparison exact matches | **8/8** (Interference+Greedy, across sequences spanning 7–14 qubits) |
| Base validation set exact matches | **12/12** |
| 44-nt challenge sequence | **F1 = 1.00**, energy gap = 0.00 kcal/mol, exact structural match |
| Helix-level qubit reduction | **40–75%** vs. pair-level encoding |

**10-nt prefix benchmark** (`GGAGCAAAAC`, unfolded MFE, E = 0.000 kcal/mol):

| Solver | Gap (kcal/mol) | F1 | Runtime (s) |
|---|---|---|---|
| Brute Force | 0.000 | 1.00 | 0.000 |
| QAOA | 0.000 | 1.00 | 1.40 |
| VQE | 3.000 | 0.00 | 0.53 |
| Interference+Greedy | 0.000 | 1.00 | 1.42 |

QAOA and Interference+Greedy both recover the exact MFE structure. VQE does not — isolated as a barren-plateau-style trainability barrier in the RY+CNOT ansatz family (confirmed by comparing against random parameter search and an all-to-all entangling topology, both of which perform comparably to COBYLA-optimized VQE).

**44-nt challenge sequence** (`GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG`):

- ViennaRNA MFE: `.(((((((..((((...(((....)))...))))..))))))).` (E = −7.90 kcal/mol)
- Hierarchical solver result: identical structure, identical energy — exact match.
- Direct pair-level encoding requires 313 qubits; the hierarchical solver reaches this result using at most 16 qubits per subproblem.

Full results, including the ablation-style solver comparison across 8 sequences, the qubit-scaling table, noise-robustness study, and simulator backend timing, are in `docs/rna_fold_qubo_writeup.pdf`.

## 7. Scaling and Resource Analysis

Qubit count under the pair-level encoding scales roughly quadratically with sequence length:

| Length (nt) | Qubits | Exceeds 16-qubit budget |
|---|---|---|
| 10 | 4 | No |
| 20 | 59 | Yes |
| 36 | 161 | Yes |
| 44 | 313 | Yes |

The helix-level encoding used by the hierarchical solver reduces this substantially (e.g. 14 → 6 qubits on a 14-nt sequence). Noise-robustness testing (DepolarizingChannel, density-matrix simulation) shows Interference+Greedy maintains exact recovery up to p = 0.50 on an 8-qubit instance; exact noise studies are restricted to ≤10 qubits by the O(4^m) memory cost of density-matrix simulation. A `lightning.qubit` vs. `default.qubit` backend comparison shows a 2.9–3.9× speedup on a 14-qubit instance.

## 8. Limitations and Next Steps

- **Hierarchical solver:** bounded top-K anchor search is a heuristic, not exhaustive backtracking. A suboptimal anchor choice at one recursion level can foreclose the correct continuation at a deeper level — on a constructed stress test of 5 sequences with deeply nested triple-helix topologies (harder than any validated sequence), the solver recovers only 1/5 exactly. Closing this gap requires exhaustive backtracking or a Nussinov/Zuker-style dynamic-programming recursion over all split points.
- **VQE:** the hardware-efficient ansatz shows a barren-plateau-style trainability barrier on this problem class; QAOA and Interference+Greedy avoid this because their circuits are built directly from the cost Hamiltonian via time evolution rather than a generic decoupled ansatz.
- **Noise studies:** limited to small instances (≤10 qubits) by the memory cost of exact density-matrix simulation.
- **Model scope:** pseudoknots are explicitly excluded (see [Section 9](#9-optional-advanced-tasks-addressed)); the Turner energy model's internal loops, bulges, and terminal mismatches beyond simple stacking are not yet modeled.

**Next steps:** validate the interference decoder on real quantum hardware; add zero-noise extrapolation / probabilistic error cancellation for noisy-device fidelity; extend the QUBO to the full Turner energy model; replace heuristic anchor search with an exact DP recursion; test the hierarchical solver on 100–200 nt sequences.

## 9. Optional Advanced Tasks Addressed

- ✅ **Pseudoknot exclusion, explicitly justified:** the no-pseudoknot penalty term scopes the model to nested (non-crossing) secondary structures — the same scope ViennaRNA's own MFE folding uses — rather than leaving the exclusion implicit.
- ✅ **Multiple encodings compared:** pair encoding vs. partner-index encoding (qubit count and constraint-enforcement tradeoffs), plus a third helix-level encoding developed specifically to close the qubit-budget gap for longer sequences.
- ✅ **Robustness under noise:** depolarizing-noise study across noise levels p = 0.00–0.50 on an 8-qubit instance.

## 10. Team

**Team QuantumCFO**

| Member | Role and contribution |
|---|---|
| Hisham Mansour | QUBO formulation and thermodynamic grounding (pair and helix-level encodings); QAOA, VQE, and Interference+Greedy implementation; hierarchical divide-and-conquer solver. |
| Yasir Mansour | Scaling, noise-robustness, and solver-comparison studies; write-up and submission package. |

## 11. Data and Privacy

All sequences used are synthetic or drawn from the publicly provided challenge example sequence. No proprietary Moderna data, patient data, clinical data, or personally identifiable information is used at any point, consistent with the challenge's data privacy and eligibility restrictions.
