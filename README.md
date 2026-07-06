# Rhodopsin search in assemblies

This repository contains a Snakemake pipeline for detecting and annotating rhodopsin-like ORFs from genome assemblies.

## Required inputs

1. Assemblies in FASTA format: `input/{assembly}.fasta`
2. Databases:
   - BLAST database sequence file: `databases/{database}.faa`
   - Marker file used for workflow discovery: `databases/{database}.faa.pdb`

## Main outputs

- Final tables: `output/{assembly}.tsv`
- Intermediate results:
  - `analysis/getorf/{assembly}.fasta`
  - `analysis/opsinmaphmm/{assembly}/`
  - `analysis/blast/{assembly}/{database}.txt`

## Deploy

```bash
snakedeploy deploy-workflow https://github.com/BejaLab/opsinsearch-assemblies project_name --branch main
```

## Run

Using the provided profile:

```bash
snakemake --cores 8
```

## Config parameters (`config/config.yaml`)

- `pos`: reference residue positions used by mapping logic
- `hmm_score`: minimum HMM score threshold
- `blast_evalue`: maximum BLAST e-value threshold
- `mapping_score`: minimum residue mapping score threshold
