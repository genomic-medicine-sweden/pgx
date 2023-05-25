__author__ = "Massimiliano Volpe"
__copyright__ = "Copyright 2023, Massimiliano Volpe & Chelsea Ramsin & Lauri Mesilaakso"
__email__ = "massimiliano.volpe@liu.se & chelsea.ramsin@regionostergotland.se & lauri.mesilaakso@regionostergotland.se"
__license__ = "GPL-3"


rule depth_of_baits:
    input:
        intervals="pgx/reform_genomic_region/get_padded_baits/padded_bait_interval.list",
        fasta=config.get("reference", {}).get("fasta", ""),
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
    output:
        gdf="pgx/depth_of_coverage/depth_of_baits/{sample}_{type}.output.gdf",
    params:
        extra=config.get("depth_of_baits", {}).get("extra", ""),
    log:
        "pgx/depth_of_coverage/depth_of_baits/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "pgx/depth_of_coverage/depth_of_baits/{sample}_{type}.output.benchmark.tsv",
            config.get("depth_of_baits", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("depth_of_baits", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("depth_of_baits", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("depth_of_baits", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("depth_of_baits", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("depth_of_baits", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("depth_of_baits", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("depth_of_baits", {}).get("container", config["default_container"])
    message:
        "{rule}: get read depth of baits on {input.bam}"
    wrapper:
        "v1.14.1/bio/gatk/depthofcoverage"


rule depth_of_targets:
    input:
        intervals="pgx/reform_genomic_region/sample_target_list/{sample}_{type}.target_interval.list",
        fasta=config.get("reference", {}).get("fasta", ""),
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
    output:
        "pgx/depth_of_coverage/depth_of_targets/{sample}_{type}.depth_at_missing.gdf",
    params:
        extra=config.get("depth_of_targets", {}).get("extra", ""),
    log:
        "pgx/depth_of_coverage/depth_of_targets/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "pgx/depth_of_coverage/depth_of_targets/{sample}_{type}.output.benchmark.tsv",
            config.get("depth_of_targets", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("depth_of_targets", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("depth_of_targets", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("depth_of_targets", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("depth_of_targets", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("depth_of_targets", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("depth_of_targets", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("depth_of_targets", {}).get("container", config["default_container"])
    message:
        "{rule}: get read depth of variant locations at wildtrype-called positions on {input.bam}"
    wrapper:
        "v1.14.1/bio/gatk/depthofcoverage"
