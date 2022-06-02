# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Massimiliano Volpe & Lauri Mesilaakso"
__copyright__ = "Copyright 2022, Massimiliano Volpe & Lauri Mesilaakso"
__email__ = "massimiliano.volpe@liu.se & lauri.mesilaakso@regionostergotland.se"
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

if os.path.isfile("config/config.yaml"):

    configfile: "config/config.yaml"


elif os.path.isfile("config.yaml"):

    configfile: "config.yaml"


elif not workflow.overwrite_configfiles:
    raise FileExistsError("No config file found in working directory or passed as argument!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    type="N|T|R",


def compile_output_list(wildcards):
    output_files = [
        "pgx/padd_target_regions/%s_%s_padded_bait_interval.bed" % (sample, t)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
    ]
    return output_files
