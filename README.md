# RNA-seq Pipeline Reproduction

[![Upstream baseline](https://img.shields.io/badge/baseline-pinned-blue)](https://github.com/sarvezh1/rna-seq-pipeline/commit/7f56d14454948e295bae1db810a0f85e995a9d28)
[![Workflow](https://github.com/Aabida1/rnaseq-pipeline/actions/workflows/upstream-reproduction.yml/badge.svg)](https://github.com/Aabida1/rnaseq-pipeline/actions/workflows/upstream-reproduction.yml)

This repository documents and tests a reproducible run of the upstream [sarvezh1/rna-seq-pipeline](https://github.com/sarvezh1/rna-seq-pipeline) using the upstream project's **bundled technical FASTQ fixtures**.

## Scope

The first milestone is technical reproduction: check out an exact upstream commit, run its included sample sheet and *Bacillus subtilis* reference using Nextflow and Docker, and retain execution reports. The bundled FASTQs are for pipeline plumbing and **must not be used for biological conclusions**.

This repository is intended to keep the reproduced upstream baseline distinguishable from any later independent modifications. It does not copy or alter the upstream pipeline's source code. The CI workflow checks out the pinned upstream commit at runtime.

## Reproduced baseline

- **Upstream:** [sarvezh1/rna-seq-pipeline](https://github.com/sarvezh1/rna-seq-pipeline)
- **Pinned commit:** [`7f56d14454948e295bae1db810a0f85e995a9d28`](https://github.com/sarvezh1/rna-seq-pipeline/commit/7f56d14454948e295bae1db810a0f85e995a9d28)
- **Workflow:** FastQC, fastp, STAR, featureCounts, MultiQC and native Nextflow execution reports, as configured by the upstream project.
- **Test data:** Bundled upstream technical fixtures; no biological interpretation.
- **Automation:** [GitHub Actions](https://github.com/Aabida1/rnaseq-pipeline/actions).

## Verified technical test run

The first GitHub Actions run completed successfully on 11 October 2026.

| Check | Observed result |
|---|---|
| Nextflow | 26.04.6 |
| Workflow profile | Docker |
| Tasks | 7 succeeded, 0 failed |
| Processes completed | INPUT_CHECK, FASTQC, FASTP, STAR_GENOMEGENERATE, STAR_ALIGN, FEATURECOUNTS, MULTIQC |
| Reported workflow duration | 46.9 seconds |
| Pinned upstream revision | `7f56d14454948e295bae1db810a0f85e995a9d28` |

Evidence: [successful GitHub Actions run #1](https://github.com/Aabida1/rnaseq-pipeline/actions/runs/38166500002). The uploaded execution trace shows exit code 0 for all seven tasks. This verifies successful execution of the small technical fixture in that runner environment only. It does not establish human RNA-seq performance, biological validity, alignment quality thresholds, or accuracy of downstream scientific conclusions.

## Run locally

Prerequisites: Git, Java compatible with Nextflow, Nextflow 26.04.6, and Docker.

```bash
git clone https://github.com/sarvezh1/rna-seq-pipeline.git upstream-rna-seq-pipeline
cd upstream-rna-seq-pipeline
git checkout 7f56d14454948e295bae1db810a0f85e995a9d28
nextflow run . \
  -profile docker \
  --input samplesheet.csv \
  --genome_fasta test_data/reference/genome.fasta \
  --gtf test_data/reference/genes.gtf \
  --outdir results \
  -work-dir /tmp/rnaseq_reproduction_work
```

Refer to the upstream README for current parameter behavior and optional branches.

## Evidence and limitations

Each GitHub Actions run is the execution record for the pinned technical fixture. Download the `upstream-baseline-reports` artifact from the workflow run to inspect retained reports and outputs. Passing the fixture demonstrates only that the tested workflow executes on that environment; it does **not** validate performance on other datasets, establish biological findings, or independently verify every optional pipeline branch.

The reproduced code remains upstream-authored. Any later changes should be added as clearly labelled, independently authored improvements with their own tests, commits, and comparison against this baseline.

## Provenance and citation

See [REPRODUCTION_NOTES.md](REPRODUCTION_NOTES.md) for the baseline commit, run strategy, and attribution.

The upstream repository is distributed under the MIT License. Review and retain upstream copyright and license notices if you later copy or redistribute any upstream code. This reproduction repository's documentation and workflow are separate from the upstream software; the upstream license should not be interpreted as a claim that the original code was authored here.
