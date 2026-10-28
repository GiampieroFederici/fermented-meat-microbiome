# Fermented-meat microbiome analysis

Code accompanying *The microbiome of fermented meat: new metagenomic data from Italian Ciauscolo salami and a comparison with worldwide products*.

## Inputs

Obtain Tables S1–S10 supplied with the manuscript and place their XLSX files in a separate directory. The tables are not duplicated here. The preparation script reads them without modifying them and generates disposable inputs locally.

This repository adds only five pre-filter abundance matrices not recoverable from the published filtered profiles, the original 229-sequence marker alignment and IQ-TREE tree, and the basemap geometry needed for Figure 1. Other metadata, label maps, palette values and derived matrices are reconstructed from the supplementary tables.

## Run

Requirements: R with readxl, vegan, ape, ade4, ggplot2, patchwork and gtable; Arial and Cairo-capable graphics. Tested with R 4.6.1, readxl 1.5.0, vegan 2.7.5, ape 5.8.1, ade4 1.7.24, ggplot2 4.0.3, patchwork 1.3.2 and gtable 0.3.6.

From the repository root, prepare the inputs using the directory containing the ten XLSX files:

```sh
Rscript --vanilla scripts/00_prepare_inputs.R /path/to/supplementary_tables .
Rscript --vanilla scripts/03_numerical_audit.R .
Rscript --vanilla scripts/05_tree_audit.R .
Rscript --vanilla scripts/09_RV_ade4_verification.R .
Rscript --vanilla scripts/11_grouped_and_host_checks.R .
Rscript --vanilla scripts/01_community_figures.R .
Rscript --vanilla scripts/08_reorder_figure4.R .
Rscript --vanilla scripts/02_circular_tree.R .
Rscript --vanilla scripts/10_full_figure_S1.R .
Rscript --vanilla scripts/13_correct_figure1.R .
```

Run 01 before 08. Figures are written to `outputs`; Figure 1 is written to `03_RETTIFICHE_PRONTE/03_MANOSCRITTO`. Generated inputs and outputs are excluded by `.gitignore`.

## Scope

The scripts regenerate Figures 1, 4, 5, 6 and S1 and selected downstream numerical results. Source values for Figures 2 and 3 are provided in the supplementary tables. Raw-read profiling, marker reconstruction and IQ-TREE inference are not rerun. Accession status is provided in Table S3.

In the inspected Windows environment, some plotting processes can terminate abnormally after producing their outputs; reopen exported files to verify them independently.
