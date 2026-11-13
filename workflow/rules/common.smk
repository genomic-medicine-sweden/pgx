__author__ = "Joel Ås, Massimiliano Volpe, Chelsea Ramsin & Lauri Mesilaakso"
__copyright__ = "Copyright 2022, Joel Ås, Massimiliano Volpe, Chelsea Ramsin & Lauri Mesilaakso"
__email__ = "massimiliano.volpe@liu.se, chelsea.ramsin@regionostergotland.se & lauri.mesilaakso@regionostergotland.se"
__license__ = "GPL-3"

import json
import linecache
import pathlib
import re
import sys
import typing

import pandas as pd
import yaml
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics import min_version as hydra_min_version
from hydra_genetics.utils.misc import extract_chr
from hydra_genetics.utils.misc import get_module_snakefile
from hydra_genetics.utils.misc import replace_dict_variables

# TODO: restore hydra_min_version("4.1.1") once a hydra-genetics release with Snakemake 9 /
# Python 3.12 support exists. The migration branch pinned in requirements.txt reports
# 3.4.1.devN, so the 4.1.1 floor cannot be asserted during the migration.
hydra_min_version("3.4.0")

min_version("9.0.0")

# Set and validate config file
if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


config = replace_dict_variables(config)

validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


# Read and validate samples file
samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

design = pd.read_table(config.get("reference", {}).get("design_bed", ""), dtype=str, index_col=0)

# Read and validate units file
units = (
    pd.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

with open(config["output"]) as output:
    if config["output"].endswith("json"):
        output_spec = json.load(output)
    elif config["output"].endswith("yaml") or config["output"].endswith("yml"):
        output_spec = yaml.safe_load(output.read())

validate(output_spec, schema="../schemas/output_files.schema.yaml")


# Reading choromosomes from design bed file
def get_choromosomes(genomic_regions: pd.DataFrame) -> typing.List[str]:
    """
    function used to extract all unique chromosomes found in design bed file
    Args:
        genomic_regions (pd.DataFrame):

    Returns:
        typing.List[str]: List of strings with all unique chromosomes found in regions
    """
    return list(set([f"{str(region.Index)}" for region in genomic_regions.itertuples()]))


# Set wildcard constraints
wildcard_constraints:
    sample="|".join(re.escape(s) for s in samples.index),
    chr="[^_]+",
    type="N|T|R",


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec.get("directory", "./"))
    output_files = []
    types = set([unit.type for unit in units.itertuples()])

    for f in output_spec["files"]:
        # Please remember to add any additional values down below
        # that the output strings should be formatted with.
        outputpaths = set(
            [
                f["output"].format(sample=sample, type=unit_type)
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
                if unit_type in set(f["types"]).intersection(types)
            ]
        )
        for op in outputpaths:
            output_files.append(outdir / pathlib.Path(op))
    return output_files


class _IdentityLineMap(dict):
    """Line mapping for generated code, whose compiled and source line numbers are the same."""

    def __contains__(self, lineno):
        return True

    def __missing__(self, lineno):
        return lineno


def generate_copy_rules(output_spec):
    output_directory = pathlib.Path(output_spec.get("directory", "./"))
    rulestrings = []

    for f in output_spec["files"]:
        if f["input"] is None:
            continue

        rule_name = "_copy_{}".format("_".join(re.split(r"\W{1,}", f["name"].strip().lower())))
        input_file = pathlib.Path(f["input"])
        output_file = output_directory / pathlib.Path(f["output"])

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        rule_code = "\n".join(
            [
                f'@workflow.rule(name="{rule_name}")',
                f'@workflow.input("{input_file}")',
                f'@workflow.output("{output_file}")',
                f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
                f'@workflow.container("{copy_container}")',
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb={mem_mb}, '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
                f'@workflow.shellcmd("{copy_container}")',
                "@workflow.run\n",
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, "
                "log, rule, conda_env, container_img, singularity_args, use_singularity, "
                "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
                "shadow_dir, edit_notebook, conda_base_path, basedir, sourcecache_path, "
                "runtime_sourcecache_path, runtime_paths, __is_snakemake_rule_func=True):",
                '\tshell("(cp {input[0]} {output[0]}) &> {log}", bench_record=bench_record, '
                "bench_iteration=bench_iteration)\n\n",
            ]
        )

        rulestrings.append(rule_code)

    code = "\n".join(rulestrings)
    source_name = "copy_result_files"

    # Snakemake 9 derives rule.run_func_src for the "code" rerun trigger by looking the compiled
    # filename up in workflow.linemaps and reading the lines back through linecache. Neither knows
    # about source we compile ourselves, so register it in both. The generated source is its own
    # original, so compiled and source line numbers coincide.
    linecache.cache[source_name] = (len(code), None, code.splitlines(True), source_name)
    workflow.linemaps[source_name] = _IdentityLineMap()

    exec(compile(code, source_name, "exec"), workflow.globals)


generate_copy_rules(output_spec)
