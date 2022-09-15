__author__ = "Lauri Mesilaakso & Massimiliano Volpe"
__copyright__ = "Copyright 2022, Lauri Mesilaakso & Massimiliano Volpe"
__email__ = "Lauri.Mesilaakso@regionostergotland.se & massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule filtering_variant_filtration:
    input:
        vcf="snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf",
    output:
        filtered_vcf="filtering/variant_filtration/{sample}_{type}_{chr}.filtered.vcf",
    params:
        read_ratio=config.get("filtering_variant_filtration", {}).get("read_ratio", ""),
        read_depth=config.get("filtering_variant_filtration", {}).get("read_depth", ""),
    log:
        "pgx/filtering_variant_filtration/{sample}_{type}_{chr}.output.log",
    benchmark:
        repeat(
            "pgx/filtering_variant_filtration/{sample}_{type}_{chr}.output.benchmark.tsv",
            config.get("filtering_variant_filtration", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("filtering_variant_filtration", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("filtering_variant_filtration", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("filtering_variant_filtration", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("filtering_variant_filtration", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("filtering_variant_filtration", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("filtering_variant_filtration", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("filtering_variant_filtration", {}).get("container", config["default_container"])
    conda:
        "../envs/filtering.yaml"
    message:
        "{rule}: Filtering vcf by read ratio and read depth on pgx/{rule}/{wildcards.sample}_{wildcards.type}_{wildcards.chr}.input"
    script:
        "../scripts/variant_filtration.py"
