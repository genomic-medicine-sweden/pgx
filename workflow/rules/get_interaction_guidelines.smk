__author__ = "Chelsea Ramsin"
__copyright__ = "Copyright 2023, Chelsea Ramsin"
__email__ = "Chelsea.Ramsin@regionostergotland.se"
__license__ = "GPL-3"


rule get_interaction_guidelines:
    input:
        diploids="pgx/get_clinical_guidelines/{sample}_{type}.output.csv",
    output:
        csv=temp("pgx/get_interaction_guidelines/{sample}_{type}.output.csv"),
    params:
        interacting_targets=config.get("get_interaction_guidelines", {}).get("interaction_guidelines"),
        extra=config.get("get_interaction_guidelines", {}).get("extra", ""),
    log:
        "pgx/get_interaction_guidelines/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "pgx/get_interaction_guidelines/{sample}_{type}.output.benchmark.tsv",
            config.get("get_interaction_guidelines", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("get_interaction_guidelines", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("get_interaction_guidelines", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("get_interaction_guidelines", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("get_interaction_guidelines", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("get_interaction_guidelines", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("get_interaction_guidelines", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("get_interaction_guidelines", {}).get("container", config["default_container"])
    message:
        "{rule}: given Haplotype Combinations, get possible interactions betweens these on {input}"
    script:
        "../scripts/get_interaction_guidelines.py"
