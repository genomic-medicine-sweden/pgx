# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Lauri Mesilaakso"
__copyright__ = "Copyright 2022, Lauri Mesilaakso"
__email__ = "lauri.mesilaakso@regionostergotland.se"
__license__ = "GPL-3"


rule padd_target_regions:
    input:
        target_regions=config.get("padd_target_regions", {}).get("target_regions", ""),
    output:
        padded_target_regions=temp("pgx/padd_target_regions/padded_bait_interval.bed"),
    params:
        padding=config.get("padd_target_regions", {}).get("padding", "100"),
    log:
        "pgx/padd_target_regions/padded_bait_interval.output.log",
    benchmark:
        repeat(
            "pgx/padd_target_regions/padded_bait_interval.output.benchmark.tsv",
            config.get("padd_target_regions", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("padd_target_regions", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("padd_target_regions", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("padd_target_regions", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("padd_target_regions", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("padd_target_regions", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("padd_target_regions", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("padd_target_regions", {}).get("container", config["default_container"])
    message:
        "{rule}: padd bed file on {input}"
    script:
        "../scripts/padd_target_regions.py"
