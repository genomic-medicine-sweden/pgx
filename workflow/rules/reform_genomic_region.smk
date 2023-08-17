__author__ = "Massimiliano Volpe"
__copyright__ = "Copyright 2023, Massimiliano Volpe"
__email__ = "massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule get_padded_bed:
    input:
        target_bed=config.get("reference", {}).get("design_bed", ""),
    output:
        interval=temp("pgx/reform_genomic_region/get_padded_bed/padded_bait_interval.bed"),
    params:
        extra=config.get("get_padded_bed", {}).get("extra", ""),
        padding=config.get("get_padded_bed", {}).get("padding", "0"),
        file_format=config.get("get_padded_bed", {}).get("file_format", "bed"),
    log:
        "pgx/reform_genomic_region/get_padded_bed/padded_bait_interval.bed.log",
    benchmark:
        repeat(
            "pgx/reform_genomic_region/get_padded_bed/padded_bait_interval.bed.benchmark.tsv",
            config.get("get_padded_bed", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("get_padded_bed", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("get_padded_bed", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bed_get_padded_bed", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("get_padded_bed", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("get_padded_bed", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("get_padded_bed", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("get_padded_bed", {}).get("container", config["default_container"])
    message:
        "{rule}: padd bed file on {input}"
    script:
        "../scripts/reform_genomic_region.py"


rule get_padded_baits:
    input:
        target_bed=config.get("reference", {}).get("design_bed", ""),
    output:
        interval=temp("pgx/reform_genomic_region/get_padded_baits/padded_bait_interval.list"),
    params:
        extra=config.get("get_padded_baits", {}).get("extra", ""),
        padding=config.get("get_padded_baits", {}).get("padding", "0"),
        file_format=config.get("get_padded_baits", {}).get("file_format", "list"),
    log:
        "pgx/reform_genomic_region/get_padded_baits/padded_bait_interval.list.log",
    benchmark:
        repeat(
            "pgx/reform_genomic_region/get_padded_baits/padded_bait_interval.benchmark.tsv",
            config.get("get_padded_baits", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("get_padded_baits", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("get_padded_baits", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("get_padded_baits", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("get_padded_baits", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("get_padded_baits", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("get_padded_baits", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("get_padded_baits", {}).get("container", config["default_container"])
    message:
        "{rule}: padd bed file on {input}"
    script:
        "../scripts/reform_genomic_region.py"


rule sample_target_list:
    input:
        target_bed=config.get("reference", {}).get("design_rsid", ""),
        detected_variants="pgx/detected_variants/{sample}_{type}.annotated.csv",
    output:
        interval=temp("pgx/reform_genomic_region/sample_target_list/{sample}_{type}.target_interval.list"),
    params:
        extra=config.get("sample_target_list", {}).get("extra", ""),
        padding=config.get("sample_target_list", {}).get("padding", "0"),
        file_format=config.get("sample_target_list", {}).get("file_format", "list"),
    log:
        "pgx/reform_genomic_region/sample_target_list/{sample}_{type}.target_interval.list.log",
    benchmark:
        repeat(
            "pgx/reform_genomic_region/sample_target_list/{sample}_{type}.target_interval.list.benchmark.tsv",
            config.get("sample_target_list", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sample_target_list", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sample_target_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sample_target_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sample_target_list", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sample_target_list", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sample_target_list", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sample_target_list", {}).get("container", config["default_container"])
    message:
        "{rule}: reform variants on {input}"
    script:
        "../scripts/reform_genomic_region.py"
