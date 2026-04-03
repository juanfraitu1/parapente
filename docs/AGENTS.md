# Repository Guidelines

## Project Structure and Module Organization
- Root scripts are standalone Python workflows (for example `shared_reads_detector.py`, `run_id14_analysis.py`, `plot_id14_summary.py`).
- Data and outputs live in `matrices/`, `results/`, `figures/`, `logs/`, and `legacy/` (archived experiments). Large `.bam`, `.fasta`, and `.bed` inputs are stored in the root and should be treated as read-only.

## Build, Test, and Development Commands
- `python3 shared_reads_detector.py --bam HSA_15.bam --seed_chrom NC_060925.1 --seed_start 121194173 --seed_end 121402237 --output results.bed` runs the production detector on a seed locus.
- `python3 shared_reads_detector.py ... --min_std_coverage 50 --max_iterations 1` enables the validated coverage filter and keeps discovery to component 1.
- `python3 run_id14_analysis.py` or `python3 run_overlap_analysis.py` run evaluation pipelines; `python3 plot_id14_summary.py` and `python3 create_presentation_figures.py` generate figures in `figures/`.

## Coding Style and Naming Conventions
- Python code uses 4-space indentation and `snake_case` for functions and variables. Script filenames are descriptive and use `snake_case`.
- No formatter or linter is configured in this repository. Follow PEP 8 where practical.
- Do not use emojis or em-dashes in any code, documentation, or outputs.

## Testing Guidelines
- There is no automated test framework in this workspace. Validation is performed against ground truth with command-line checks.
- Use `bedtools intersect` to validate new family detections against `ground_truth_80families.bed` as described in `CLAUDE.md`.

## Commit and Pull Request Guidelines
- Git history is not available in this workspace, so no existing commit convention can be summarized.
- Use concise, imperative commit subjects. PRs should link issues and summarize scripts, datasets, and outputs (figures or BED/TSV files).

## Data and Execution Notes
- The pipeline assumes BAM files are indexed (`.bai` present) and that large input files may be stored outside the repository as documented in `CLAUDE.md`.
- Prefer `--chromosomes` and `--max_iterations 1` for faster runs when appropriate.

## Documentation Summary
- `DOCUMENTATION_FINDINGS_AND_COMMANDS.md` consolidates findings on overextension and locus trimming/splitting, with copy-paste commands to reproduce ID_8 and ID_131 runs.
- Key conclusions: Jaccard flank trim helps but does not fix overextension alone; `--split-loci` yields gene-sized cores, while `jaccard_greedy` loses recall. Supplementary alignments do not materially affect results.
