__author__ = "Massimiliano Volpe"
__copyright__ = "Copyright 2023, Massimiliano Volpe"
__email__ = "massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule append_id_to_gdf:
    input:
        gdf="pgx/depth_of_coverage/depth_of_targets/{sample}_{type}_{chr}.depth_at_missing.gdf",
    output:
        out="pgx/append_id_to_gdf/{sample}_{type}_{chr}.depth_at_missing_annotated.gdf",
    params:
        extra=config.get("append_id_to_gdf", {}).get("extra", ""),
        target_bed=config.get("append_id_to_gdf", {}).get("target_rsid", ""),
        #module=workflow.source_path("../scripts/get_target_variants.py"),
    log:
        "pgx/append_id_to_gdf/{sample}_{type}_{chr}.output.log",
    benchmark:
        repeat(
            "pgx/append_id_to_gdf/{sample}_{type}_{chr}.output.benchmark.tsv",
            config.get("append_id_to_gdf", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("append_id_to_gdf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("append_id_to_gdf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("append_id_to_gdf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("append_id_to_gdf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("append_id_to_gdf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("append_id_to_gdf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("append_id_to_gdf", {}).get("container", config["default_container"])
    conda:
        "../envs/append_id_to_gdf.yaml"
    message:
        "{rule}: add variant id to appropriate location in {input}"
    script:
        "../scripts/append_rsid_to_gdf.py"
