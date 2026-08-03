# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Code for the experiments in "Mapping the latent CRBN-molecular glue degrader interactome". It searches
a proteome-scale structural database for surface patches that geometrically/chemically mimic a known
protein-protein (or protein-ligand) interface, using MaSIF (surface fingerprint) descriptors, and then
scores and filters the resulting binder candidates ("bonsai" postprocessing/truncation workflow).

This repo (`masif_mimicry`) is **not standalone**: it is meant to be cloned *inside* a checkout of
[`masif_seed`](https://github.com/LPDI-EPFL/masif_seed) (see README.md Installation). Code under
`source/masif_mimicry/` imports `default_config.masif_opts` and `geometry.open3d_import`, which live in
the sibling `masif`/`masif_seed_search` source trees, not in this repo — always PYTHONPATH-injected at
runtime (see "Environments" below). Do not try to resolve these imports locally without that PYTHONPATH.

## Environments

There are two separate, incompatible runtime environments used by this pipeline. Getting them mixed up
is the most common source of import errors:

1. **MaSIF/Apptainer image** (`masif_mimicry.sif`, built from `Dockerfile`, `pymesh/pymesh` base): TensorFlow
   1.12, Open3D 0.8, PyMesh, MSMS, APBS, PDB2PQR, USalign, `reduce`. Used for steps 1–3 (PDB preprocessing,
   target site definition, mimicry search). Always invoked via `apptainer exec ... "${IMAGE}" python -u -m ...`.
2. **`MaSIF` conda env** (separate from the Apptainer image; see `scripts/bonsai_scripts/`): EvoEF2, STRIDE,
   DeepTMHMM, BioPython, RDKit. Used for step 4 (postprocessing metrics) and the bonsai truncation workflow.
   Created via `scripts/bonsai_scripts/slurm/conda.sh` / `install_dependencies.sh` and
   `scripts/bonsai_scripts/environment.yml`.

Search-step code (`search/`) additionally shells out to `/install/USalign/USalign` inside the Apptainer image
(hardcoded path in `search/workflow.py`), so it will not run outside that image.

## Pipeline architecture

The pipeline is a 4-stage SLURM workflow, all steps sharing config from `mimicry_pipeline.config.sh`
(sourced by every `*.slurm` file and by `submit_mimicry_pipeline.sh`). Edit that one file to point at a
new target/database/output location; do not hardcode paths in the slurm scripts.

```
1_preprocess_pdb.slurm        -> preprocess_pdb.sh (in Apptainer)
2_prepare_target_sites.slurm  -> python -m masif_mimicry.define_target_sites
3_run_masif_mimicry.slurm     -> python -m masif_mimicry.search           (SLURM array job)
4_postprocess_mimicry.slurm   -> python -m masif_mimicry.postprocess      (SLURM array job, conda env MaSIF)
```

`submit_mimicry_pipeline.sh` reads `RUN_PREPROCESS` / `RUN_TARGET_SITES` / `RUN_MIMICRY_SEARCH` /
`RUN_POSTPROCESS` booleans from the config to decide which steps to `sbatch`, chaining them with
`--dependency=afterok:<jobid>`. Steps 3 and 4 are SLURM array jobs, one array task per seed subset file
under `DATABASE_SUBSET_DIR/<array_id>` (a plain text file of seed PDB ids, one per line).

### Naming convention: P1 (seed) vs P2 (target)

Throughout `source/masif_mimicry/`, **P1** is a candidate/seed protein from the searched database and
**P2** is the fixed target protein whose interface you're trying to mimic. A "MaSIF PDB identifier" like
`021structure_C_AB` encodes `<pdb>_<p1chain(s)>_<p2chain(s)>` (see
`utils/parse_structures.target_chain_from_pdb_id` / `parse_partner_chain_ids`); single-domain seeds use
just `<pdb>_<chain>`. `p1`/`p2` "ppi_id" arguments select which side of that encoding to treat as source
vs. target for a given protein.

### Stage details

- **Step 1 (`preprocess_pdb.sh`)**: converts a raw PDB (+ optional ligand SDF/mol2) into MaSIF's expected
  `data_preparation/` tree (`00-raw_pdbs`, `01-benchmark_pdbs`, `01-benchmark_surfaces`,
  `04a-precomputation_9A`, `04b-precomputation_12A`, plus MaSIF-site interface predictions and
  ppi_search descriptors) — this is the same layout used for both the target protein and the seed
  database. `config/paths.py::set_params`/`resolve_database_paths` is what knows how to navigate this
  tree given either a database root or a `data_preparation` dir directly.
- **Step 2 (`define_target_sites.py`)**: picks patch centers ("sites") on the target surface to search
  against — either around a specific residue/atom (`--target_residue`/`--target_atom`), exhaustively
  (with `--downsample`/`--iface_cutoff`/`--top_iface_percent`/`--interface_only`), or from manually
  specified patch indices (`--site_idx`). Writes a shared `target_sites.json` manifest plus
  `target_vert/` patch files into `--target_run_dir`, so step 3 doesn't repeat target-side feature
  extraction per array task.
- **Step 3 (`search/workflow.py::run_search`)**: for each seed PDB, loads MaSIF descriptors
  (`search/features.get_features`), filters candidate patches by descriptor distance
  (`--desc_dist_cutoff`) against the target's patch descriptors, aligns candidates with Open3D
  ICP/RANSAC registration (`search/docking.py`), scores alignments
  (`search/scoring.compute_descriptor_score`, `--desc_dist_score_cutoff`), optionally filters clashing
  poses against the target's partner chains (`--count_clashes`), and clusters structurally similar hits
  (`utils/clustering.structural_clusters`). Writes one CSV per seed under
  `<target_run_dir>/<P1_id>/`.
- **Step 4 (`postprocess/`)**: `discovery.discover_deduplicated_rows` collects/dedupes the per-seed CSVs
  from step 3, then `pipeline.process_results_mimicry` computes per-hit metrics (clash counts, SASA,
  interface residues, iface aggregates in `postprocess/metrics/`) and merges in domain metadata from
  `--database_info_csv`, streaming results to one combined output CSV.

### Downstream: bonsai postprocessing/truncation workflow

`scripts/bonsai_scripts/` (see its own README) is a separate, later-stage SLURM workflow that takes the
best hits from step 4's output, computes additional metrics (secondary structure via STRIDE, membrane
topology via DeepTMHMM, geodesic length), and truncates binders to a target length using EvoEF2 folding
free energy calculations to find optimal truncation windows. Configured via
`scripts/bonsai_scripts/python/config.yaml` and launched with
`bash scripts/bonsai_scripts/postprocessing_truncation_workflow_apptainer.sh`. This uses the `MaSIF`
conda env, not the Apptainer image directly.

Other one-off scripts in `scripts/`:
- `domain_split.py` / `process_af_model.sh`: split an AlphaFold model into structured domains by PAE
  cutoff (used to build a searchable domain-level database from whole proteins).
- `results_explorer_app.py`: a Dash/py3Dmol app for interactively browsing a step-4 results CSV
  (3D view of target + matched binder poses). Has hardcoded paths at the top of the file that must be
  edited per-user/per-run before launching.

## Common commands

There is no build step, test suite, linter, or package manifest in this repo — it's a collection of
Python modules run via `python -m`, invoked from SLURM scripts inside the Apptainer image or conda env.
When making changes, verify by running the relevant stage's `python -m masif_mimicry.<stage>` command
directly (with `--help` to check argument wiring) rather than by unit tests.

```bash
# Edit shared config first, then submit whichever pipeline steps are enabled there:
vim mimicry_pipeline.config.sh   # toggle RUN_PREPROCESS / RUN_TARGET_SITES / RUN_MIMICRY_SEARCH / RUN_POSTPROCESS
./submit_mimicry_pipeline.sh                # submit
./submit_mimicry_pipeline.sh --overwrite    # also pass --overwrite to step 2 (define_target_sites)

# Monitor
squeue -u "$USER"
tail -f logs/masif_*.out
tail -f logs/mimicry_prepare_*.out
tail -f logs/mimicry_<jobid>_*.out
tail -f logs/postprocess_<jobid>_*.out

# Run a single stage manually inside the Apptainer image, e.g. to debug step 3 on one seed:
apptainer exec --bind /scratch:/scratch --bind /work:/work masif_mimicry.sif \
  python -u -m masif_mimicry.search --database_dir <db> --target_preprocess_dir <dir> \
  --target_run_dir <dir> --seed_pdb <PDB_chain> --desc_dist_cutoff 2.0 --desc_dist_score_cutoff 0.45
```

Rebuilding the Apptainer image (rarely needed):
```bash
apptainer pull masif_mimicry.sif docker://shxiao/masif_mimicry:latest
```

## Key tunable parameters (in `mimicry_pipeline.config.sh`)

- `DESC_DIST_CUTOFF` (recommended 1.5–2.0, lower = stricter/faster): descriptor-distance prefilter before
  expensive alignment.
- `DESC_DIST_SCORE_CUTOFF` (recommended 0.35–0.5, higher = stricter): post-alignment MaSIF mimicry score
  cutoff; hits below this are discarded.
- `iface_cutoff`/`--iface_cutoff` (0.0–0.95, higher = stricter): interface-propensity cutoff for patch
  selection; for degron interfaces MaSIF wasn't trained on ternary complexes, so `0.0` is recommended.
- `NUM_POINTS`/`num_sites`: number of target sites to search from (usually fine at low single digits).
- `TARGET_DOWNSAMPLE`/`downsample`: stride over candidate patches; increase to speed up large-database
  searches at the cost of recall.
