# SpiPy

Computational pipeline for spiropyran photoswitches: **SMILES in → descriptors out**.

Manages HPC job submission (PBS via shell wrappers) for CREST conformer search, ORCA two-tier DFT, and Multiwfn wavefunction analysis on Metacentrum. All heavy computation runs on the cluster; this package handles input generation, job tracking, and output parsing.

See [`insilico-screening.md`](insilico-screening.md) for full project specification and [`PLAN.md`](PLAN.md) for the implementation plan.

---

## Working in an existing environment

```bash
module add mambaforge
conda activate spipy
pytest tests/ -v
```

---

## First-time setup

```bash
module add mambaforge
mamba env create -f environment.yml   # creates the spipy environment
conda activate spipy
pytest tests/ -v
```

After pulling changes that add new dependencies:

```bash
mamba env update -f environment.yml --prune
```
