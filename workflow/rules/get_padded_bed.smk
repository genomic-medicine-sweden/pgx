__author__ = "Massimiliano Volpe"
__copyright__ = "Copyright 2023, Massimiliano Volpe"
__email__ = "massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule get_padded_bed:
    input:
        target_bed=config.get("get_padded_bed", {}).get("target_regions", ""),
    output:
        interval="pgx/get_padded_bed/padded_bait_interval.bed",
    params:
        padding=config.get("get_padded_bed", {}).get("padding", "0"),
        file_format=config.get("get_padded_bed", {}).get("file_format", "bed"),
        extra=config.get("get_padded_bed", {}).get("extra", ""),
    log:
        "pgx/get_padded_bed/padded_bait_interval.output.bed.log",
    benchmark:
        repeat(
            "pgx/get_padded_bed/padded_bait_interval.bed.output.benchmark.tsv",
            config.get("get_padded_bed", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("get_padded_bed", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("get_padded_bed", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bed_get_padded_bed", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("get_padded_bed", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("get_padded_bed", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("get_padded_bed", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("get_padded_bed", {}).get("container", config["default_container"])
    conda:
        "../envs/get_padded_bed.yaml"
    message:
        "{rule}: padd bed file on {input}"
    script:
        "../scripts/reform_genomic_region.py"
