# pRosettaC — PROTAC Ternary Complex Modelling Pipeline

A computational pipeline for modelling PROTAC-induced ternary complexes between a kinase (POI), the VHL E3 ligase, and the bifunctional degrader molecule **SJF8240**. The goal is to predict which kinases are geometrically compatible with productive ubiquitination and therefore likely to be degraded.

---

## Background

PROTACs (PROteolysis TArgeting Chimeras) are bifunctional small molecules that recruit an E3 ubiquitin ligase to a protein of interest (POI), inducing its ubiquitination and subsequent proteasomal degradation. Whether degradation actually occurs depends critically on the geometry of the resulting ternary complex — the PROTAC must bridge the two proteins in a way that places a lysine residue on the POI within reach of the E2 enzyme (UBE2D).

This pipeline uses **[pRosettaC](https://github.com/LiorZ/pRosettaC)** (Zaidman et al. 2020) to sample ternary complex poses via global docking and constraint-guided local optimisation with Rosetta, then analyses the resulting poses for biological plausibility.

---

## Repository Structure

```
pRosettaC/
├── Structure_Prep/          # Step 1 — prepare input structures
│   ├── 01_download_structures.py   # Download and chain-extract PDB files
│   ├── 02_prepare_full_kinase.py   # Fix residue-number overlaps for RIPK2/MAPK14
│   └── 03_place_warhead.py         # Place SJF8240 warhead into each kinase ATP site
│
├── Rosetta_utils/           # Step 2 — run pRosettaC
│   ├── 01_install_rosetta          # Rosetta + pRosettaC installation guide
│   ├── 02_sds_to_params.py         # Generate Rosetta .params from SDF (RDKit-based)
│   ├── 03_run_local_prosettac.py   # Run pRosettaC locally (no cluster required)
│   ├── 04_run_all_pois.py          # Batch-run all POIs sequentially
│   ├── 05_gen_local_fasc.py        # Score docking PDBs → local.fasc with PyRosetta
│   └── 06_resume_pipeline.py       # Resume from constraint_generation + clustering
│
└── Analysis /               # Step 3 — analyse poses
    ├── 01_find_feasible_poses.py   # Scan Patchdock_Results for bridging poses
    ├── 02_recluster_feasible.py    # Re-cluster bridging poses by Cα RMSD
    ├── 03_batch_check_protac.py    # Quick PASS/MARGINAL/FAIL screen across all POIs
    └── 04_run_poi_full_analysis.py # Full feature extraction + pre-MD report
```

---

## Prerequisites

| Dependency | Purpose |
|---|---|
| [Rosetta](https://www.rosettacommons.org/software/license-and-download) | Global docking and local optimisation (free academic licence) |
| [pRosettaC](https://github.com/LiorZ/pRosettaC) | PROTAC-specific ternary complex sampling |
| [PyRosetta](https://www.pyrosetta.org) | Python interface to Rosetta (scoring, params) |
| RDKit | SDF→params conversion, warhead coordinate handling |
| BioPython | Sequence alignment, Cα superposition |
| NumPy / SciPy | Distance matrices, hierarchical clustering |

All Python dependencies can be installed into a conda environment:

```bash
conda create -n protac python=3.11
conda activate protac
conda install -c conda-forge rdkit biopython scipy numpy
pip install pyrosetta-installer   # or install PyRosetta wheel directly
```

Set environment variables before running:

```bash
export ROSETTA_HOME=/path/to/rosetta_src_<version>
export ROSETTA_BIN=$ROSETTA_HOME/main/source/bin
export PROSETTAC_HOME=~/pRosettaC       # upstream pRosettaC repo
```

---

## Pipeline Walkthrough

### Step 1 — Structure Preparation

```bash
# 1a. Download PDB files and extract the relevant chain for each POI and VHL
python Structure_Prep/01_download_structures.py

# 1b. Fix residue-number overlaps for RIPK2 and MAPK14
#     (their N-lobe numbering clashes with VHL residues 1-103)
python Structure_Prep/02_prepare_full_kinase.py

# 1c. Place the SJF8240 kinase warhead into each POI ATP site via MCS + superposition
python Structure_Prep/03_place_warhead.py
# Add --force <TARGET> to re-generate specific targets
```

After this step each POI directory contains:
- `<POI>_protein_clean.pdb` — cleaned kinase structure
- `<POI>_warhead_placed.sdf` — warhead in the ATP-site frame (no H)
- `<POI>_warhead_placed_H.sdf` — same with explicit hydrogens

### Step 2 — Rosetta / pRosettaC

```bash
# 2a. Generate Rosetta .params for the PROTAC (if not already done)
python Rosetta_utils/02_sds_to_params.py data/structures/PROTAC/SJF8240/SJF8240.sdf SJF data/ligand_params/SJF8240.params

# 2b. Run pRosettaC for a single POI (from that POI's run directory)
cd pRosettaC/runs/high_affinity_degraders/DDR2
python /path/to/Rosetta_utils/03_run_local_prosettac.py Protac_params.txt \
    --n_dock 200 --n_workers 8

# 2c. Or batch-run all configured POIs
python Rosetta_utils/04_run_all_pois.py

# 2d. If the pipeline was interrupted after docking but before clustering
python Rosetta_utils/06_resume_pipeline.py   # run from the POI run directory
```

Output per POI lands in `pRosettaC/runs/<category>/<POI>/`:
- `Patchdock_Results/` — raw docking PDBs (`pd.*_docking_????.pdb`) and constraint-generated ternary complexes (`combined_*_0001.pdb`)
- `Results/cluster*/` — pRosettaC clustered representatives

### Step 3 — Analysis

```bash
# 3a. Quick bridging screen: PASS / MARGINAL / FAIL for each POI
python "Analysis /03_batch_check_protac.py" DDR2 MET EPHB2
python "Analysis /03_batch_check_protac.py" --all

# 3b. Find all bridging poses (writes <poi>_feasible.txt to $HOME)
python "Analysis /01_find_feasible_poses.py" DDR2 MET

# 3c. Re-cluster the bridging poses by Cα RMSD and pick top-3 candidates for MD
python "Analysis /02_recluster_feasible.py" DDR2

# 3d. Full feature extraction + pre-MD structural report
python "Analysis /04_run_poi_full_analysis.py" DDR2
python "Analysis /04_run_poi_full_analysis.py" --all
```

Outputs of `04_run_poi_full_analysis.py`:
| File | Contents |
|---|---|
| `cluster_features.csv` | Per-cluster ML feature table (contacts, clashes, lysine geometry) |
| `cluster_features.json` | Same data in nested JSON |
| `premd_descriptors.csv` | Features for the top-3 MD candidates |
| `premd_report.txt` | Human-readable pre-MD structural check (chain completeness, clash counts, PROTAC geometry, verdict) |

---

## Configured POIs

| Target | Category | Notes |
|---|---|---|
| MET | high_affinity_degraders | |
| DDR2 | high_affinity_degraders | |
| RIPK2 | high_affinity_degraders | N-lobe absent; renumbered to res 200+ (`_full` run) |
| EPHB2 | high_affinity_degraders | |
| MAPK14 | high_affinity_degraders | N-lobe absent; renumbered to res 200+ (`_full` run) |
| AXL | high_affinity_no_degradation | |
| ABL1 | high_affinity_no_degradation | |
| EPHA2 | high_affinity_no_degradation | |
| MAP4K5 | high_affinity_no_degradation | |
| SLK | high_affinity_no_degradation | |
| TNIK | low_affinity_degraders | |
| PIP4K2C | low_affinity_degraders | Non-kinase phosphatase |

---

## Bridging Criteria

A ternary complex pose is scored:

| Verdict | Criterion |
|---|---|
| **PASS** | Closest PROTAC atom < 5 Å from both VHL (chain A res 1–103) and POI |
| **MARGINAL** | Closest PROTAC atom < 10 Å from both |
| **FAIL** | Outside both thresholds |

For N-lobe–truncated kinases (RIPK2, MAPK14, SLK, TNIK, MAP4K5) a stricter `poi_check_min` residue lower bound is used to avoid a systematic artefact where the VHL–kinase junction creates trivial false-positive bridging poses.

---

## Key Design Notes

- **Rosetta chain layout:** VHL = chain A res 1–103, POI = chain A (original or renumbered), ElonginB = chain B, ElonginC = chain C, PROTAC = chain X.
- **Score shifting:** pRosettaC's `clustering.py` filters out scores ≥ 0. If PyRosetta's `PackRotamersMover` produces positive scores, `06_resume_pipeline.py` shifts the entire `score.sc` by `max_score + 100` before clustering.
- **N-lobe renumbering:** RIPK2 and MAPK14 have kinase domain numbering that overlaps with VHL (res 1–103). `02_prepare_full_kinase.py` adds an offset (+192 / +196) so the full domain is retained without silently losing the N-lobe.

---

## Citation

If you use pRosettaC in your work, please cite:

> Zaidman, D. et al. *PROTACable integrates modeling and deep learning to enable data-driven PROTAC design.* **Nature Chemical Biology** (2023).
> Zaidman, D. et al. *PRosettaC: Rosetta based modeling of PROTAC mediated ternary complexes.* **J. Chem. Inf. Model.** 60, 4894–4903 (2020).
