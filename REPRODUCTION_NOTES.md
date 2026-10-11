# Reproduction Notes

## Objective

Re-run the upstream RNA-seq workflow on its own bundled technical test data before using a biological dataset. The intent is to establish a traceable baseline that later changes can be compared against.

## Source pin

- Project: [sarvezh1/rna-seq-pipeline](https://github.com/sarvezh1/rna-seq-pipeline)
- Exact commit: [`7f56d14454948e295bae1db810a0f85e995a9d28`](https://github.com/sarvezh1/rna-seq-pipeline/commit/7f56d14454948e295bae1db810a0f85e995a9d28)
- Commit message at the time pinned: `docs: add citation metadata and Zenodo DOI`
- License declared upstream: MIT; see the upstream [LICENSE](https://github.com/sarvezh1/rna-seq-pipeline/blob/7f56d14454948e295bae1db810a0f85e995a9d28/LICENSE).

The CI job checks out this exact commit rather than tracking the moving `main` branch. Update the pin only in a dedicated, documented change.

## Technical fixture

The upstream README describes the included FASTQs as small technical fixtures, not data for biological interpretation. The supplied reference is *Bacillus subtilis* subsp. *subtilis* strain 168 (assembly GCF_000009045.1 / ASM904v1). The expected result of this milestone is evidence that pipeline stages execute and reports are produced, not a biological conclusion.

## Protocol

1. Set up Java, Nextflow 26.04.6, and Docker on the GitHub-hosted runner.
2. Check out the pinned upstream commit into a separate directory.
3. Run the upstream workflow using the bundled sample sheet and matching FASTA/GTF.
4. Confirm the expected native Nextflow report files are created.
5. Upload the report directory as a GitHub Actions artifact.

## Separation of baseline and improvements

- This repository does not vendor or rewrite upstream pipeline files.
- The original upstream repository is not modified by this project.
- A passing CI run is recorded by GitHub Actions and the report artifact.
- Future improvements should have explicit commits and tests; retain the upstream pin and document any deviations.

## Interpretation limits

A test fixture can verify workflow integration, file hand-offs, and report generation in a controlled small run. It does not prove that the pipeline produces valid biological results for human samples, that all optional branches work, or that downstream scientific assumptions are appropriate for a new study.

## Upstream citation

The upstream project provides its own citation metadata in [CITATION.cff](https://github.com/sarvezh1/rna-seq-pipeline/blob/7f56d14454948e295bae1db810a0f85e995a9d28/CITATION.cff). If you use the upstream software in a report, paper, or portfolio, consult that file and cite the upstream author/software. Do not claim the upstream implementation as original authorship.
