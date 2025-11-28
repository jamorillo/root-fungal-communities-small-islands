# Root-associated fungal communities on five distant small islands

This repository contains the QIIME 2 code used in the analyses of:

**Morillo, J.A., Traveset, A., Correia, M., Hervías-Parejo, S., Heleno, R., Nogales, M., Hernández, M.M., Kaiser-Bunbury, C., Rodríguez-Echeverría, S.**  
*Root-associated fungal communities on five distant small islands.* (submitted, Fungal Ecology)

## Overview

The project explores root-associated fungal communities on five small, geographically distant islands, using ITS1 metabarcoding data.

This repository includes:

- QIIME 2 pipeline for ITS1 (single-end R1, ITS1F–ITS2)
- Export to R

## Repository structure

- `scripts/`
  - `ITS1_qiime2_pipeline.sh` – main QIIME 2 pipeline
- `metadata/`
- `docs/`

## QIIME 2 ITS1 pipeline

The main pipeline is implemented in:

- `scripts/ITS1_qiime2_pipeline.sh`

Key features:

- **Region**: ITS1 of fungi (EMP ITS1F / ITS2 primers)
- **Reads**: single-end R1 (Illumina)
- **Primer trimming**: `q2-cutadapt` (forward primer + reverse-complement of reverse primer)
- **Denoising / ASV inference**: `q2-dada2 denoise-single´
- **Taxonomic assignment**:
  - UNITE v9.0 (dynamic, all taxa)
- **Filtering**:
  - Keep only Fungi
  - Keep ASVs classified at phylum level
- **Export**:
  - ASV table and representative sequences

## Requirements

- QIIME 2 (tested with:
  - `qiime2-amplicon-2024.5`
- `seqkit` (optional, for primer-checking)
- UNITE classifiers in QIIME 2 format



