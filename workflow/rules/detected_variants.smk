__author__ = "Massimiliano Volpe"
__copyright__ = "Copyright 2022, Massimiliano Volpe"
__email__ = "massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule detected_variants:
    input:
        vcf="pgx/variant_annotator/{sample}_{type}_{chr}.output.vcf",
    output:
        csv="pgx/detected_variants/{sample}_{type}_{chr}.output.csv",
    params:
        target_bed=config.get("detected_variants", {}).get("target_rsid", ""),
    log:
        "pgx/detected_variants/{sample}_{type}_{chr}.output.log",
    benchmark:
        repeat(
            "pgx/detected_variants/{sample}_{type}_{chr}.output.benchmark.tsv",
            config.get("detected_variants", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("detected_variants", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("detected_variants", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("detected_variants", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("detected_variants", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("detected_variants", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("detected_variants", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("detected_variants", {}).get("container", config["default_container"])
    conda:
        "../envs/detected_variants.yaml"
    message:
        "{rule}: get variants with target rsIDs on {input}"
    script:
        "../scripts/get_target_variants.py"
