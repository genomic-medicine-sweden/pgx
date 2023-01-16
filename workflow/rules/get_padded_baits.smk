__author__ = "Massimiliano Volpe"
__copyright__ = "Copyright 2022, Massimiliano Volpe"
__email__ = "massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule get_padded_baits:
    input:
        target_bed=config.get("get_padded_baits", {}).get("target_regions", ""),
    output:
        interval="pgx/get_padded_baits/padded_bait_interval.list",
    params:
        padding=config.get("get_padded_baits", {}).get("padding", "0"),
        file_format=config.get("get_padded_baits", {}).get("file_format", "list"),
    log:
        "pgx/get_padded_baits/padded_bait_interval.output.log",
    benchmark:
        repeat(
            "pgx/get_padded_baits/padded_bait_interval.output.benchmark.tsv",
            config.get("get_padded_baits", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("get_padded_baits", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("get_padded_baits", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("get_padded_baits", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("get_padded_baits", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("get_padded_baits", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("get_padded_baits", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("get_padded_baits", {}).get("container", config["default_container"])
    conda:
        "../envs/get_padded_baits.yaml"
    message:
        "{rule}: padd bed file on {input}"
    script:
        "../scripts/reform_genomic_region.py"
