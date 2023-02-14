__author__ = "Massimiliano Volpe"
__copyright__ = "Copyright 2022, Massimiliano Volpe"
__email__ = "massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule detected_variants:
    input:
        input_f="pgx/variant_annotator/{sample}_{type}_{chr}.annotated.vcf",
    output:
        output_f="pgx/detected_variants/{sample}_{type}_{chr}.annotated.csv",
    params:
        extra=config.get("detected_variants", {}).get("extra", ""),
        target_bed=config.get("reference", {}).get("design_rsid", ""),
        file_type=config.get("detected_variants", {}).get("file_type", ""),
    log:
        "pgx/detected_variants/{sample}_{type}_{chr}.annotated.csv.log",
    benchmark:
        repeat(
            "pgx/detected_variants/{sample}_{type}_{chr}.annotated.csv.benchmark.tsv",
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
        "../envs/get_target_variants.yaml"
    message:
        "{rule}: get variants with target rsIDs on {input}"
    script:
        "../scripts/get_target_variants.py"


rule append_id_to_gdf:
    input:
        input_f="pgx/depth_of_coverage/depth_of_targets/{sample}_{type}_{chr}.depth_at_missing.gdf",
    output:
        output_f="pgx/append_id_to_gdf/{sample}_{type}_{chr}.depth_at_missing_annotated.gdf",
    params:
        extra=config.get("append_id_to_gdf", {}).get("extra", ""),
        target_bed=config.get("reference", {}).get("design_rsid", ""),
        file_type=config.get("append_id_to_gdf", {}).get("file_type", ""),
    log:
        "pgx/append_id_to_gdf/{sample}_{type}_{chr}.output.log",
    benchmark:
        repeat(
            "pgx/append_id_to_gdf/{sample}_{type}_{chr}.depth_at_missing_annotated.gdf.benchmark.tsv",
            config.get("append_id_to_gdf", {}).get("benchmark_repeats", 1),
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
        "../envs/get_target_variants.yaml"
    message:
        "{rule}: add variant id to appropriate location in {input}"
    script:
        "../scripts/get_target_variants.py"
