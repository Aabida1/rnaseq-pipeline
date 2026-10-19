# OncoGene Explorer

A separate Streamlit research app for exploratory breast-cancer differential-expression review and cross-cohort evidence comparison.

The app is maintained as a standalone project. Keep its independent analysis, runtime dependencies, and release lifecycle separate from this RNA-seq pipeline reproduction repository.

## Source project package

The current OncoGene Explorer source package is maintained in the project Library and can be added as a standalone repository once a repository is created for it. Do not place large processed datasets or GEO platform annotation files in Git; download or upload those separately as documented by the app.

## Core features

- Cohort-specific volcano plots and DEG tables
- GSE183947 RNA-seq and uploaded cohort DEG results such as GSE42568 microarray
- Exploratory gene matching and cross-cohort evidence summaries
- Optional GMT pathway enrichment and source-traceable target evidence

## Scientific caveats

The application summarizes precomputed differential-expression results; it does not independently refit each cohort model. RNA-seq and microarray cohorts must be analyzed separately, and gene identifiers and contrast directions must be harmonized before interpreting cross-cohort evidence. Outputs are exploratory and are not clinical diagnoses or validated biomarker claims.
