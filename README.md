# Analysis of single-cell neuron organoids data

This repository contains some code for the analysis of single-cell data of neuron organoids

## Files description

- `annotate_seu_objects.R`: R script for annotating clusters of all integrated [Seurat](https://satijalab.org/seurat/) objects

- `cluster-celltypes-neurorg-48`: TSV file containing manually annotated clusters for integrated samples at time 48

- `cluster-celltypes-neurorg-62`: TSV file containing manually annotated clusters for integrated samples at time 62

- `cluster-celltypes-neurorg-mutant`: TSV file containing manually annotated clusters for integrated samples with mutant genotype

- `cluster-celltypes-neurorg-WT`: TSV file containing manually annotated clusters for integrated samples with WT genotype

- `integration_genotype`: R script for integrating samples by genotype

- `integration_time`: R script for integrating samples by time

- `markers-general`: TSV file containing general marker genes for neurons and corresponding cell types

- `markers-specific`: TSV file containing specific marker genes for neurons and corresponding cell types

- `report_meeting_annotation.Rmd`: R Markdown template for generating a HTML report with all the markers expression (useful for cell type manual annotation)

- `report_meeting_genotype.Rmd`: R Markdown template for generating a HTML report with all the analysis of the samples integrated by genotype

- `report_meeting_time.Rmd`: R Markdown template for generating a HTML report with all the analysis of the samples integrated by time