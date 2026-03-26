# SpiPy — Implementation Plan

## Project Goal

Automated computational pipeline: **SMILES in → descriptors out** for spiropyran photoswitches. Given a spiropyran SMILES, the pipeline produces both closed (SP) and open (merocyanine, MC) forms, generates conformers, runs quantum chemistry calculations on an HPC cluster, and assembles a descriptor DataFrame for downstream ML.

Combinatorial library enumeration (Fischer indole, Sigma catalogue) is a later-stage layer on top of this core pipeline.

---

## Architecture Overview

```
SMILES string
    │
    ▼
Mol validation + SP 3D embedding
    │
    ├──► SP (closed form) geometries
    │
    └──► MC (open form) via automated ring-opening SMARTS
              │
              ▼
         Enumerate all E/Z isomers around formally double bonds
              │
              ▼
         Constrained 3D embedding per isomer
              │
              ▼
         Steric clash filter → valid MC geometries
    │
    ▼
Chemoinformatic descriptors (Mordred, morfeus, DScribe) ← local, no HPC
    │
    ▼
CREST conformer search (HPC via shell wrappers)
    │
    ▼
ORCA Tier 1 (small basis + implicit solvent) → energy window cutoff
    │
    ▼
ORCA Tier 2 (larger basis / better theory) → final wavefunctions
    │
    ▼
Multiwfn extraction (ESP, polarisability, charges, QTAIM)
    │
    ▼
Merged descriptor DataFrame (one row per molecule, all levels)
```

---

## Package Layout

```
spipy/
├── __init__.py
├── config.py            # TOML config loading + defaults
├── mol.py               # SMILES validation, SP↔MC transform, 3D embedding
├── isomers.py           # MC E/Z enumeration, steric clash filtering
├── runner.py            # Generic HPC job runner (shell wrapper interface)
├── crest.py             # CREST-specific: input prep, output parsing
├── orca.py              # ORCA .inp writing, two-tier logic, output parsing (cclib)
├── multiwfn.py          # Multiwfn batch scripts, text dump parsing
├── descriptors.py       # Mordred, morfeus, DScribe wrappers
├── pipeline.py          # Orchestrator: SMILES → full descriptor row
└── cli.py               # CLI entry point
tests/
├── conftest.py          # Shared fixtures (test SMILES, small molecules)
├── test_mol.py
├── test_isomers.py
├── test_runner.py
├── test_crest.py
├── test_orca.py
├── test_multiwfn.py
├── test_descriptors.py
└── test_pipeline.py
```

---

## Key Design Decisions

- **HPC submission:** The package does NOT manage PBS directly. User has existing wrapper scripts that accept xyz + basic parameters and handle all PBS decoration. We shell out to those wrappers.
- **Job tracking:** File-based JSON manifest per molecule. Each entry tracks: molecule ID, isomer label, tool (CREST/ORCA/Multiwfn), status (pending/running/done/failed), input/output file paths.
- **Persistence:** xyz files for geometries, Parquet for descriptor matrices, JSON for job state.
- **Config:** TOML file with per-stage sections (ORCA method/basis per tier, energy window thresholds, clash distance cutoff, etc.).
- **E/Z isomer handling:** Isolated module. Generates all 2^n combinations, filters by steric clash. Designed to be swappable if a smarter enumeration strategy is needed later.
- **Resumability:** Each pipeline stage is idempotent. If output files exist and are valid, skip. The orchestrator can resume from any stage.

---

## Work Packages

Each WP is scoped for a 1–3 hour co-coding session. Dependencies are listed — do not start a WP before its dependencies are complete.

---

### WP1 — Project Skeleton + Molecule Handling

**Depends on:** nothing

**Objective:** Establish the package, config system, and core molecule I/O. A worker finishing this WP should produce a package that can be installed (`pip install -e .`), with passing tests.

**Deliverables:**
1. `pyproject.toml` with dependencies: `rdkit`, `numpy`, `pandas`, `tomli` (or `tomllib` on 3.11+), `pytest`
2. `spipy/config.py` — load a TOML config file, provide sensible defaults. Config sections: `[mol]`, `[isomers]`, `[crest]`, `[orca]`, `[multiwfn]`, `[descriptors]`, `[pipeline]`.
3. `spipy/mol.py`:
   - `validate_smiles(smiles: str) -> Mol` — parse, sanitize, return RDKit Mol or raise
   - `embed_3d(mol: Mol) -> Mol` — ETKDG embedding + MMFF optimization
   - `mol_to_xyz(mol: Mol, path: Path) -> Path` — write xyz file
   - `sp_to_mc(mol: Mol) -> Mol` — SMARTS-based ring-opening transform (spiropyran → merocyanine). This is the most chemically demanding part of WP1. The SMARTS must break the C–O bond at the spiro centre and adjust bond orders in the resulting conjugated chain. Test with at least 2 known spiropyrans.
4. `tests/conftest.py` — fixtures with 2–3 known spiropyran SMILES (e.g., unsubstituted 1',3',3'-trimethylindolinospiropyran and one nitro-substituted variant)
5. `tests/test_mol.py` — validate round-trip SMILES→Mol→SMILES, 3D embedding produces reasonable geometry, SP→MC produces correct open form

**Notes for worker:**
- The SP→MC transform is the hardest part. The spiro C–O bond breaks, the pyran ring opens, and the result is a conjugated chain with formally alternating single/double bonds. Study the spiropyran ring-opening mechanism before writing the SMARTS. The worker should verify the product by visual inspection (draw the molecule) and by checking expected atom count / bond order changes.
- Do not over-engineer the config system. A simple dataclass with a `from_toml()` classmethod is sufficient.

---

### WP2 — Merocyanine E/Z Isomer Generation

**Depends on:** WP1

**Objective:** Given a merocyanine Mol object (from `sp_to_mc`), enumerate all E/Z isomers around the formally double bonds in the conjugated chain, embed each in 3D with constrained dihedrals, and filter out sterically impossible geometries.

**Deliverables:**
1. `spipy/isomers.py`:
   - `find_mc_double_bonds(mol: Mol) -> list[tuple[int, int]]` — identify atom index pairs of formally double bonds in the MC conjugated chain that are relevant for E/Z isomerism. Not all double bonds — only those in the open chain between the indoline and the phenolate/quinone moiety.
   - `enumerate_ez(mol: Mol, bonds: list[tuple[int, int]]) -> list[Mol]` — generate all 2^n combinations, setting each bond to E or Z. Return list of Mol objects with appropriate stereochemistry set.
   - `embed_constrained(mol: Mol) -> Mol` — 3D embedding respecting the assigned E/Z stereochemistry. Use RDKit constrained embedding or set dihedral angles post-embedding.
   - `filter_clashes(mols: list[Mol], min_distance: float) -> list[Mol]` — remove molecules where any non-bonded atom pair is closer than `min_distance` (configurable, default ~1.0 Å or based on vdW radii). Alternatively, use MMFF energy: discard if energy is above threshold relative to the lowest-energy isomer.
   - `generate_mc_isomers(mol: Mol) -> list[Mol]` — top-level function chaining the above steps
2. `tests/test_isomers.py` — test that known spiropyran produces expected number of valid isomers, that impossible geometries are filtered, that E/Z labels are correctly assigned

**Notes for worker:**
- The merocyanine conjugated chain typically has 3 bonds with E/Z freedom, giving 8 initial combinations. After steric filtering, expect 3–5 survivors (TTT, TTC, CCC are the commonly observed ones; the exact set depends on substitution).
- This module must be self-contained and swappable. If the steric clash approach doesn't work well, we may replace it with an explicit list of allowed configurations. Design the interface so that's easy.
- Label each isomer with its configuration string (e.g., "TTC") for downstream tracking.

---

### WP3 — Generic HPC Job Runner

**Depends on:** WP1

**Objective:** Provide a uniform interface for submitting jobs to the HPC cluster via the user's existing shell wrapper scripts, tracking their status, and collecting results.

**Deliverables:**
1. `spipy/runner.py`:
   - `Job` dataclass: molecule_id, isomer_label, tool (str), status (enum: pending/submitted/running/done/failed), input_path, output_path, submit_time, completion_time
   - `JobManifest` — read/write a JSON file tracking all jobs for a pipeline run. Support querying by status, molecule, tool.
   - `submit_job(command: list[str], job: Job) -> Job` — run `subprocess.run()` (or `Popen` for background) to call the user's wrapper script. Update job status. The command is constructed by the caller (CREST/ORCA/Multiwfn modules), this function just executes it and tracks state.
   - `poll_completion(manifest: JobManifest, check_file: Callable[[Job], bool]) -> JobManifest` — iterate over submitted jobs, check if expected output file exists (caller provides the check function), update status.
   - `wait_for_jobs(manifest: JobManifest, check_file: Callable, poll_interval: int) -> JobManifest` — blocking poll loop with configurable interval. Optional, for interactive use.
2. `tests/test_runner.py` — test manifest CRUD, test submit with a mock command (e.g., `touch output.xyz`), test polling logic

**Notes for worker:**
- Keep this generic. It knows nothing about CREST, ORCA, or Multiwfn. It runs shell commands and checks for output files.
- The user's wrapper scripts handle all PBS specifics (queue, walltime, modules). We just call them.
- Error handling: if the wrapper script returns non-zero exit code, mark job as failed and store stderr.
- Do NOT implement retry logic. If a job fails, the user investigates manually.

---

### WP4 — CREST Conformer Search

**Depends on:** WP2, WP3

**Objective:** For each valid MC isomer (and the SP closed form), run CREST via the HPC runner and parse the resulting conformer ensemble.

**Deliverables:**
1. `spipy/crest.py`:
   - `prepare_crest_input(mol: Mol, work_dir: Path, constraints: dict | None) -> Path` — write xyz file and any constraint file to work_dir. Return path to xyz.
   - `build_crest_command(xyz_path: Path, work_dir: Path, config: dict) -> list[str]` — construct the shell command to call the user's CREST wrapper. Config includes: solvent, method flags, constraint references.
   - `parse_crest_output(work_dir: Path) -> list[tuple[Mol, float]]` — read `crest_conformers.xyz` (multi-geometry xyz file), return list of (Mol, energy) pairs sorted by energy.
   - `filter_conformers(conformers: list[tuple[Mol, float]], energy_window: float) -> list[tuple[Mol, float]]` — keep conformers within energy_window (kJ/mol) of the lowest.
2. `tests/test_crest.py` — test xyz writing, command construction, output parsing (use a fixture `crest_conformers.xyz` file with known content)

**Notes for worker:**
- The multi-geometry xyz format: repeated blocks of `natoms\ncomment (energy)\nx y z coords...`. Parser must handle this robustly.
- Constraint handling for merocyanine: we may want to fix certain dihedrals during CREST to preserve the E/Z configuration. This is optional for v1 — flag it with a TODO if not implemented.
- Provide a small sample `crest_conformers.xyz` as a test fixture.

---

### WP5 — ORCA Two-Tier Pipeline

**Depends on:** WP3, WP4

**Objective:** Write ORCA input files, manage the two-tier optimization workflow, and parse results.

**Deliverables:**
1. `spipy/orca.py`:
   - `write_orca_input(mol: Mol, path: Path, config: dict) -> Path` — write `.inp` file. Config specifies: method, basis set, solvent model, job type (opt, sp, opt+freq), extra keywords. Geometry from Mol's 3D conformer.
   - `OrcaTier` dataclass: method, basis, solvent, job_type, extra_keywords
   - `build_orca_command(inp_path: Path, config: dict) -> list[str]` — shell command for user's ORCA wrapper
   - `parse_orca_output(out_path: Path) -> dict` — use cclib to extract: total energy, orbital energies (HOMO/LUMO/gap), dipole moment, partial charges, optimized geometry. Return as flat dict.
   - `tier1_filter(results: list[dict], energy_window: float) -> list[dict]` — keep structures within energy_window of the lowest Tier 1 energy
   - `run_two_tier(conformers: list[Mol], tier1: OrcaTier, tier2: OrcaTier, runner: JobManifest, energy_window: float)` — orchestrate: submit Tier 1 for all → wait → filter → submit Tier 2 for survivors → wait → parse
2. `tests/test_orca.py` — test .inp file generation (string comparison against expected), test output parsing (fixture .out file), test tier filtering logic

**Notes for worker:**
- ORCA .inp format: `! method basis solvent_keyword\n* xyz charge mult\ncoords\n*`. See ORCA manual for exact syntax.
- cclib is the primary parser. If cclib doesn't extract something we need, fall back to regex on the .out file.
- The two-tier orchestration is the most complex control flow. Keep it linear and explicit — no async, no clever abstractions. Submit Tier 1, wait, filter, submit Tier 2, wait, return.
- Charge and multiplicity: default to 0/1 (neutral singlet) but make configurable.

---

### WP6 — Multiwfn Extraction

**Depends on:** WP3, WP5

**Objective:** Run Multiwfn analyses on ORCA output files and parse the text dumps into Python data structures.

**Deliverables:**
1. `spipy/multiwfn.py`:
   - `write_multiwfn_script(analyses: list[str], input_file: Path, output_dir: Path) -> Path` — write a Multiwfn batch input script. Analyses are identified by name (e.g., "esp_surface", "polarisability", "hirshfeld_charges", "fukui", "bond_orders"). Each maps to a sequence of Multiwfn menu commands.
   - `build_multiwfn_command(script_path: Path, wfn_path: Path, config: dict) -> list[str]` — shell command for user's Multiwfn wrapper
   - `parse_multiwfn_output(output_dir: Path, analysis: str) -> dict` — parse a specific analysis text dump. Each analysis type has its own parser.
   - `extract_all(output_dir: Path, analyses: list[str]) -> dict` — run all parsers, merge into single dict
2. `tests/test_multiwfn.py` — test script generation, test each parser against fixture text dumps

**Notes for worker:**
- Multiwfn is menu-driven (stdin commands). Batch mode: pipe a sequence of numbers/commands via stdin. The script file is what gets piped in.
- Each analysis type produces different text output. Parsers will be regex-based. Start with: ESP surface statistics, polarisability tensor, Hirshfeld charges. Others can be added incrementally.
- Provide fixture text dumps for each analysis type in `tests/fixtures/`.
- This is inherently fragile (text parsing). Design parsers to fail loudly with clear error messages if the format doesn't match.

---

### WP7 — Chemoinformatic Descriptors (local, no HPC)

**Depends on:** WP1, WP2

**Objective:** Compute chemoinformatic, steric, and local-environment descriptors from RDKit 3D geometries. No HPC needed.

**Deliverables:**
1. `spipy/descriptors.py`:
   - `compute_mordred(mol: Mol) -> dict` — compute Mordred descriptors, return as dict. Handle NaN/missing values (Mordred returns `mordred.error.Missing` for some descriptors). Drop those.
   - `compute_morfeus(mol: Mol, spiro_idx: int, neighbor_idxs: list[int]) -> dict` — buried volume, Sterimol L/B1/B5 at the spiro centre. Requires atom indices.
   - `compute_soap(mol: Mol, center_idx: int) -> np.ndarray` — DScribe SOAP descriptor centred on the spiro carbon.
   - `find_spiro_center(mol: Mol) -> int` — SMARTS-based identification of the spiro carbon atom index.
   - `compute_all_descriptors(mol: Mol) -> dict` — combine all local descriptors into one dict
2. `tests/test_descriptors.py` — test each descriptor function returns expected types/shapes, test spiro centre finder on known molecules

**Notes for worker:**
- Mordred is a heavy dependency (~1800 descriptors). Some will be NaN or constant across the dataset. Don't filter here — just compute and return. Filtering happens downstream in ML feature selection.
- morfeus needs coordinates as numpy array + atom indices. Extract from RDKit Mol conformer.
- DScribe SOAP needs species list + positions as numpy arrays. Convert from RDKit.
- The spiro centre SMARTS: look for a carbon bonded to atoms in two distinct ring systems. This may need refinement — test on fixtures.

---

### WP8 — Pipeline Orchestrator + CLI

**Depends on:** WP1–WP7

**Objective:** Wire everything together. Single entry point that takes SMILES and produces a full descriptor DataFrame.

**Deliverables:**
1. `spipy/pipeline.py`:
   - `Pipeline` class:
     - `__init__(config_path: Path)` — load config
     - `run(smiles: str, work_dir: Path) -> pd.DataFrame` — full pipeline execution. Steps: validate → SP embed → MC transform → E/Z enumerate → local descriptors → CREST → ORCA Tier 1 → filter → ORCA Tier 2 → Multiwfn → merge all → return DataFrame
     - Stage-aware: check work_dir for existing outputs, skip completed stages
     - Each stage writes its outputs to a subdirectory of work_dir
   - `PipelineState` — tracks which stages are complete per molecule/isomer, persisted as JSON in work_dir
2. `spipy/cli.py`:
   - `spipy run <SMILES> --config config.toml --work-dir ./runs/` — run pipeline for one molecule
   - `spipy run --file molecules.csv --config config.toml --work-dir ./runs/` — batch mode
   - `spipy status --work-dir ./runs/` — show pipeline state (what's done, what's pending on HPC)
3. `tests/test_pipeline.py` — integration test with mock runner (no real HPC), verify full pipeline produces expected DataFrame shape

**Notes for worker:**
- The orchestrator is glue code. Keep it simple and linear. No DAG engine, no task queues.
- Batch mode: just loop over molecules. No parallelism within the Python process — HPC handles that.
- The status command is important for usability: user submits jobs, goes away, comes back later and checks status.

---

## Config File Structure

```toml
[mol]
# No options needed yet

[isomers]
clash_min_distance = 1.0  # Angstrom, for steric clash filter
# Or: clash_energy_threshold = 500.0  # kJ/mol relative to lowest

[crest]
wrapper_command = "/path/to/crest_wrapper.sh"
solvent = "chcl3"
energy_window = 30.0  # kJ/mol

[orca]
wrapper_command = "/path/to/orca_wrapper.sh"
charge = 0
multiplicity = 1

[orca.tier1]
method = "B3LYP"
basis = "def2-SVP"
solvent = "CPCM(CHCl3)"
job_type = "opt"
energy_window = 20.0  # kJ/mol, for tier1→tier2 filter

[orca.tier2]
method = "B3LYP"
basis = "def2-TZVP"
solvent = "CPCM(CHCl3)"
job_type = "opt freq"

[multiwfn]
wrapper_command = "/path/to/multiwfn_wrapper.sh"
analyses = ["esp_surface", "polarisability", "hirshfeld_charges"]

[descriptors]
mordred = true
morfeus = true
soap = true

[pipeline]
poll_interval = 60  # seconds
```

---

## Dependency Summary

```
WP1 (skeleton + mol)
 │
 ├── WP2 (E/Z isomers) ──────────────┐
 │    │                                │
 │    ├── WP7 (local descriptors)      │
 │    │                                │
 │    └── WP4 (CREST)                  │
 │         │                           │
 │         └── WP5 (ORCA two-tier)     │
 │              │                      │
 │              └── WP6 (Multiwfn)     │
 │                                     │
 ├── WP3 (job runner) ────────────────(used by WP4, WP5, WP6)
 │
 └── WP8 (orchestrator + CLI) ────────(after all above)
```

WP3 and WP2 can be developed in parallel after WP1.
WP7 can be developed in parallel with WP4–WP6 (no HPC dependency).

---

## Testing Strategy

- **Unit tests:** Each module has its own test file. Use small fixture molecules (2–3 SMILES) and fixture output files (sample CREST xyz, ORCA .out, Multiwfn text dumps).
- **No real HPC in tests:** Mock the shell wrapper calls in `runner.py`. Tests verify input generation and output parsing, not actual computation.
- **Integration test:** One test in `test_pipeline.py` runs the full pipeline with a mocked runner, verifying the DataFrame has expected columns and shape.
- **Fixture files:** `tests/fixtures/` contains sample output files from CREST, ORCA, and Multiwfn for parser testing.
