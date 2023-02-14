# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Joel Ås, Massimiliano Volpe, Chelsea Ramsin & Lauri Mesilaakso"
__copyright__ = "Copyright 2022, Joel Ås, Massimiliano Volpe, Chelsea Ramsin & Lauri Mesilaakso"
__email__ = "massimiliano.volpe@liu.se, chelsea.ramsin@regionostergotland.se & lauri.mesilaakso@regionostergotland.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics import min_version as hydra_min_version

hydra_min_version("0.11.0")

min_version("6.8.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

design = pd.read_table(config.get("reference", {}).get("design_bed", ""), dtype=str, index_col=0)

### Read and validate units file

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

### Reading choromosomes from design bed file


def get_choromosomes(genomic_regions: pd.DataFrame) -> typing.List[str]:
    """
    function used to extract all unique chromosomes found in design bed file
    Args:
        genomic_regions (pd.DataFrame):

    Returns:
        typing.List[str]: List of strings with all unique chromosomes found in regions
    """
    return list(set([f"{str(region.Index)}" for region in genomic_regions.itertuples()]))


### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    chr="[^_]+",
    type="N|T|R",


def compile_output_list(wildcards):
    output_files = ["pgx/reform_genomic_region/get_padded_bed/padded_bait_interval.bed"]
    output_files += ["pgx/reform_genomic_region/get_padded_baits/padded_bait_interval.list"]
    output_files += [
        "alignment/samtools_extract_reads/%s_%s_%s.bam" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += ["snv_indels/bed_split/design_bedfile_%s.bed" % (c) for c in get_choromosomes(design)]
    output_files += [
        "alignment/picard_mark_duplicates/%s_%s_%s.bam" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += [
        "snv_indels/haplotypecaller/%s_%s_%s.vcf" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += [
        "pgx/variant_filtration/%s_%s_%s.filtered.vcf" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += [
        "pgx/variant_annotator/%s_%s_%s.annotated.vcf" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += [
        "pgx/detected_variants/%s_%s_%s.annotated.csv" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += [
        "pgx/reform_genomic_region/sample_target_list/%s_%s_%s.target_interval.list" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += [
        "pgx/depth_of_coverage/depth_of_baits/%s_%s_%s.output.gdf" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += [
        "pgx/depth_of_coverage/depth_of_targets/%s_%s_%s.depth_at_missing.gdf" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += [
        "pgx/get_clinical_guidelines/%s_%s_%s.output.csv" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += [
        "pgx/get_interaction_guidelines/%s_%s_%s.output.csv" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    output_files += [
        "pgx/append_id_to_gdf/%s_%s_%s.depth_at_missing_annotated.gdf" % (sample, t, c)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for c in get_choromosomes(design)
    ]
    return output_files
