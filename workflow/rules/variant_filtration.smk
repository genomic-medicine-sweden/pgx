__author__ = "Lauri Mesilaakso & Massimiliano Volpe"
__copyright__ = "Copyright 2022, Lauri Mesilaakso & Massimiliano Volpe"
__email__ = "Lauri.Mesilaakso@regionostergotland.se & massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule variant_filtration:
    input:
        vcf="snv_indels/haplotypecaller/{sample}_{type}.merged.vcf",
    output:
        filtered_vcf="pgx/variant_filtration/{sample}_{type}.filtered.vcf",
    params:
        read_ratio=config.get("variant_filtration", {}).get("read_ratio", ""),
        read_depth=config.get("variant_filtration", {}).get("read_depth", ""),
    log:
        "pgx/variant_filtration/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "pgx/variant_filtration/{sample}_{type}.output.benchmark.tsv",
            config.get("variant_filtration", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("variant_filtration", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("variant_filtration", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("variant_filtration", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("variant_filtration", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("variant_filtration", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("variant_filtration", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("variant_filtration", {}).get("container", config["default_container"])
    message:
        "{rule}: filter vcf by read ratio and read depth on {input}"
    script:
        "../scripts/variant_filtration.py"
