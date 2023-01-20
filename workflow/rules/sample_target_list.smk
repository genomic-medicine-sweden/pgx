__author__ = "Massimiliano Volpe"
__copyright__ = "Copyright 2022, Massimiliano Volpe"
__email__ = "massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule sample_target_list:
    input:
        target_bed=config.get("sample_target_list", {}).get("target_rsid", ""),
        detected_variants="pgx/detected_variants/{sample}_{type}_{chr}.output.csv",
    output:
        interval="pgx/sample_target_list/{sample}_{type}_{chr}.output.tsv",
    params:
        padding=config.get("sample_target_list", {}).get("padding", "0"),
        file_format=config.get("sample_target_list", {}).get("file_format", "list"),
    log:
        "pgx/sample_target_list/{sample}_{type}_{chr}.output.log",
    benchmark:
        repeat(
            "pgx/sample_target_list/{sample}_{type}_{chr}.output.benchmark.tsv",
            config.get("sample_target_list", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sample_target_list", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sample_target_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sample_target_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sample_target_list", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sample_target_list", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sample_target_list", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sample_target_list", {}).get("container", config["default_container"])
    conda:
        "../envs/sample_target_list.yaml"
    message:
        "{rule}: reform variants on {input}"
    script:
        "../scripts/reform_genomic_region.py"
