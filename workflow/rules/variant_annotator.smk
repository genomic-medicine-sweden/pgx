__author__ = "Massimiliano Volpe"
__copyright__ = "Copyright 2022, Massimiliano Volpe"
__email__ = "massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule variant_annotator:
    input:
        vcf="filtering/variant_filtration/{sample}_{type}_{chr}.filtered.vcf",
        aln="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        ref=config.get("reference", {}).get("fasta", ""),
        db=config.get("reference", {}).get("dbsnp", ""),
    output:
        vcf="pgx/variant_annotator/{sample}_{type}_{chr}.output.vcf",
    params:
        extra="",  # optional "--resource-allele-concordance -A Coverage --expression db.END",
        java_opts="",  # optional
    log:
        "pgx/variant_annotator/{sample}_{type}_{chr}.output.log",
    benchmark:
        repeat(
            "pgx/variant_annotator/{sample}_{type}_{chr}.output.benchmark.tsv",
            config.get("variant_annotator", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("variant_annotator", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("variant_annotator", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("variant_annotator", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("variant_annotator", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("variant_annotator", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("variant_annotator", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("variant_annotator", {}).get("container", config["default_container"])
    conda:
        "../envs/variant_annotator.yaml"
    message:
        "{rule}: Annotating vcf on pgx/{rule}/{wildcards.sample}_{wildcards.type}_{wildcards.chr}.input"
    wrapper:
        "v1.14.1/bio/gatk/variantannotator"