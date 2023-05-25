__author__ = "Chelsea Ramsin"
__copyright__ = "Copyright 2023, Chelsea Ramsin"
__email__ = "chelsea.ramsin@regionostergotland.se"
__license__ = "GPL-3"


rule get_clinical_guidelines:
    input:
        found_variants="pgx/detected_variants/{sample}_{type}.annotated.csv",
    output:
        csv="pgx/get_clinical_guidelines/{sample}_{type}.output.csv",
    params:
        haplotype_definitions=config.get("get_clinical_guidelines", {}).get("haplotype_definitions", ""),
        clinical_guidelines=config.get("get_clinical_guidelines", {}).get("clinical_guidelines", ""),
        haplotype_activity=config.get("get_clinical_guidelines", {}).get("haplotype_activity", ""),
        hidden_haplotypes=config.get("get_clinical_guidelines", {}).get("hidden_haplotypes", ""),
        extra=config.get("get_clinical_guidelines", {}).get("extra", ""),
    log:
        "pgx/get_clinical_guidelines/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "pgx/get_clinical_guidelines/{sample}_{type}.output.benchmark.tsv",
            config.get("get_clinical_guidelines", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("get_clinical_guidelines", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("get_clinical_guidelines", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("get_clinical_guidelines", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("get_clinical_guidelines", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("get_clinical_guidelines", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("get_clinical_guidelines", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("get_clinical_guidelines", {}).get("container", config["default_container"])
    message:
        "{rule}: given detected variants, get possible Haplotype combinations on {input}"
    script:
        "../scripts/get_clinical_guidelines.py"
