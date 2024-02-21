__author__ = "Chelsea Ramsin"
__copyright__ = "Copyright 2023, Chelsea Ramsin"
__email__ = "Chelsea.Ramsin@regionostergotland.se"
__license__ = "GPL-3"


rule get_clinical_recommendations:
    input:
        found_variants="pgx/detected_variants/{sample}_{type}.annotated.csv",
        missed_variants="pgx/append_id_to_gdf/{sample}_{type}.depth_at_missing_annotated.gdf",
        clinical_guidelines="pgx/get_clinical_guidelines/{sample}_{type}.output.csv",
        interactions="pgx/get_interaction_guidelines/{sample}_{type}.output.csv",
    output:
        recommendations=temp("pgx/get_clinical_recommendations/{sample}_{type}_pgx_clinical_recommendations.txt"),
    params:
        haplotype_definitions=config.get("get_clinical_guidelines", {}).get("haplotype_definitions", ""),
        analysed_variants=config.get("get_clinical_recommendations", {}).get("analysed_variants", ""),
    log:
        "pgx/get_clinical_recommendations/{sample}_{type}_pgx_clinical_recommendations.txt.output.log",
    benchmark:
        repeat(
            "pgx/get_clinical_recommendations/{sample}_{type}_pgx_clinical_recommendations.txt.benchmark.tsv",
            config.get("get_clinical_recommendations", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("get_clinical_recommendations", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pgx_clinical_recommendations", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pgx_clinical_recommendations", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pgx_clinical_recommendations", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pgx_clinical_recommendations", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pgx_clinical_recommendations", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("pgx_clinical_recommendations", {}).get("container", config["default_container"])
    message:
        "{rule}: Generates clinical recommendations per sample pgx/{rule}/{wildcards.sample}_{wildcards.type}.input"
    script:
        "../scripts/get_clinical_recommendations.py"
