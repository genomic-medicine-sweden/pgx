__author__ = "Chelsea Ramsin"
__copyright__ = "Copyright 2023, Chelsea Ramsin"
__email__ = "chelsea.ramsin@regionostergotland.se"
__license__ = "GPL-3"


rule generate_pgx_report:
    input:
        missed_variants="pgx/append_id_to_gdf/{sample}_{type}.depth_at_missing_annotated.gdf",
        recommendations="pgx/get_clinical_recommendations/{sample}_{type}_pgx_clinical_recommendations.txt",
        detected_variants="pgx/detected_variants/{sample}_{type}.annotated.csv",
    output:
        html_report=temp("pgx/generate_pgx_report/{sample}_{type}_pgx_report.html"),
        depth_table=temp("pgx/generate_pgx_report/{sample}_{type}_pgx_depth_table.xlsx"),
    params:
        haplotype_definitions=config.get("get_clinical_guidelines", {}).get("haplotype_definitions", ""),
        html_template=config.get("generate_pgx_report", {}).get("html_template", ""),
        dbsnp=config.get("reference", {}).get("dbsnp", ""),
        ref=config.get("reference", {}).get("fasta", ""),
        read_depth=config.get("variant_filtration", {}).get("read_depth", ""),
    log:
        "pgx/generate_pgx_report/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "pgx/generate_pgx_report/{sample}_{type}.output.benchmark.tsv",
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
    message:
        "{rule}: Generate pgx report on pgx/{rule}/{wildcards.sample}_{wildcards.type}.input"
    script:
        "../scripts/generate_pgx_report.py"
