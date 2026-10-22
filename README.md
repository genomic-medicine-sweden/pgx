# pgx

Pharmacogenomics pipeline implemented in hydra

![Lint](https://github.com/genomic-medicine-sweden/pgx/actions/workflows/lint.yaml/badge.svg?branch=develop)
![Snakefmt](https://github.com/genomic-medicine-sweden/pgx/actions/workflows/snakefmt.yaml/badge.svg?branch=develop)

![pycodestyle](https://github.com/genomic-medicine-sweden/pgx/actions/workflows/pycodestyl.yaml/badge.svg?branch=develop)
![pytest](https://github.com/genomic-medicine-sweden/pgx/actions/workflows/pytest.yaml/badge.svg?branch=develop)

![integration test](https://github.com/genomic-medicine-sweden/pgx/actions/workflows/integration1.yaml/badge.svg?branch=develop)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## :speech_balloon: Introduction

## :heavy_exclamation_mark: Dependencies

To run this workflow, the following tools need to be available:

![python](https://img.shields.io/badge/python-3.8-blue)
[![snakemake](https://img.shields.io/badge/snakemake-6.8.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![singularity](https://img.shields.io/badge/singularity-3.7-blue)](https://sylabs.io/docs/)

## :school_satchel: Preparations

### Sample data

1. Add all sample ids to `samples.tsv` in the column `sample`.
2. Add all sample data information to `units.tsv`. Each row represents a `fastq` file pair with
corresponding forward and reverse reads. Also indicate the sample id, run id and lane number, adapter.

### Reference data

1. You need a BAM file marked for duplicates
2. Reference genome
3. dbSNP database in VCF file format

`data/dbsnp/dbsnp_138.pgx_targets.vcf.gz` is dbSNP 138 (hg19) limited to the PGX target regions in
`data/genomic_regions/exons_variants_pharmacogenomics_18_06_2019_ex_cyp2d6.bed`, padded by 200 bp:

```bash
awk 'BEGIN{OFS="\t"} !/^#/ {s=$2-200; if (s<0) s=0; print $1, s, $3+200}' $BED > pgx_targets_pad200.bed
bcftools view -T pgx_targets_pad200.bed -Oz -o dbsnp_138.pgx_targets.vcf.gz dbsnp_138.vcf.gz
tabix -p vcf dbsnp_138.pgx_targets.vcf.gz
```

The padding covers the `get_padded_bed` padding (100 bp) used for variant calling.

## :white_check_mark: Testing

The workflow repository contains a small test dataset `.tests/integration` which can be run like so:

```bash
cd .tests/integration
snakemake -s ../../Snakefile -j1 --use-singularity

# alternative command:
snakemake -s ../../workflow/Snakefile -j1 --use-singularity --configfile config/config.yaml

# generate DAG:
snakemake --cores 1 -s workflow/Snakefile --configfile config/config.yaml --rulegraph | dot -Tsvg > ./images/dag.svg
```

## :rocket: Usage

The workflow is designed for WGS data meaning huge datasets which require a lot of compute power. For
HPC clusters, it is recommended to use a cluster profile and run something like:

```bash
snakemake -s /path/to/Snakefile --profile my-awesome-profile
```

### Without internet access

The hydra-genetics module versions (git tags) are set under `modules` in `config/config.yaml`. The modules are fetched from GitHub by default. To use local clones instead, set
`hydra_local_path` in the config to a folder with clones of `prealignment`, `alignment` and `snv_indels`.
The module versions are read from the local clones by tag, so the clones must contain the tags set under
`modules`.

```yaml
hydra_local_path: "/path/to/hydra-genetics"  # with /path/to/hydra-genetics/prealignment etc.
```

The reference paths in `config/config.yaml` are absolute. To use another location, add
`config/config.data.hg19.yaml` and set `PROJECT_REF_DATA` to a folder with `ref_data/hg19/hg19.with.mt.*`:

```bash
snakemake --configfiles config/config.yaml config/config.data.hg19.yaml --config PROJECT_REF_DATA=/path/to/references ...
```

## :judge: Rule Graph

![rule_graph](https://github.com/genomic-medicine-sweden/pgx/blob/develop/images/dag.svg)
