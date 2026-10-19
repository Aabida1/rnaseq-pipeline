# OncoGene Explorer

A Streamlit research application for exploratory breast-cancer differential-expression review and cross-cohort evidence comparison.

> Note: this repository is primarily the RNA-seq pipeline reproduction project. This subdirectory currently documents the OncoGene Explorer app; the runnable application source remains a separate project package and has not yet been copied into this repository.

## Planned features

- Cohort-specific volcano plots and DEG tables
- GSE183947 RNA-seq and user-uploaded cohort DEG results, including GSE42568 microarray
- Exploratory gene matching and cross-cohort evidence summaries
- Optional GMT pathway enrichment and source-traceable target evidence

## Scientific caveats

The app summarizes precomputed differential-expression results; it does not independently refit each cohort model. RNA-seq and microarray cohorts must be analyzed separately, and gene identifiers and contrast directions must be harmonized before interpreting cross-cohort evidence. Outputs are exploratory and are not clinical diagnoses or validated biomarker claims.

## Status

This README is a project pointer only. The GitHub connector available in this workflow does not expose repository creation or binary archive upload. To publish OncoGene Explorer as a standalone repository, create `Aabida1/OncoGene-Explorer` on GitHub, then push the project source and required small metadata files from the downloaded ZIP. Avoid committing the full processed-expression CSV if repository size or data redistribution terms are a concern; document how to obtain it instead.
