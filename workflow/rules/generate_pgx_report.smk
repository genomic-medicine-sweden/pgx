__author__ = "Chelsea Ramsin"
__copyright__ = "Copyright 2023, Chelsea Ramsin"
__email__ = "Chelsea.Ramsin@regionostergotland.se"
__license__ = "GPL-3"


rule generate_pgx_report:
    input:
        found_variants="pgx/detected_variants/{sample}_{type}.annotated.csv",
        missed_variants="pgx/append_id_to_gdf/{sample}_{type}.depth_at_missing_annotated.gdf",
        clinical_guidelines="pgx/get_clinical_guidelines/{sample}_{type}.output.csv",
        interactions="pgx/get_interaction_guidelines/{sample}_{type}.output.csv",
    output:
        report="pgx/generate_pgx_report/{sample}_{type}_pgx_report.txt",
    params:
        haplotype_definitions=config.get("get_clinical_guidelines", {}).get("haplotype_definitions", ""),
        report_template=config.get("generate_pgx_report", {}).get("report_template", ""),
    log:
        "pgx/generate_pgx_report/{sample}_{type}_pgx_report.txt.output.log",
    benchmark:
        repeat(
            "pgx/generate_pgx_report/{sample}_{type}_pgx_report.txt.benchmark.tsv",
            config.get("generate_pgx_report", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("generate_pgx_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("generate_pgx_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("generate_pgx_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("generate_pgx_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("generate_pgx_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("generate_pgx_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("generate_pgx_report", {}).get("container", config["default_container"])
    conda:
        "../envs/generate_pgx_report.yaml"
    message:
        "{rule}: Generates markdown report per sample pgx/{rule}/{wildcards.sample}_{wildcards.type}.input"
    script:
        "../scripts/generate_pgx_report.py"
