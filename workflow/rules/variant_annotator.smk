__author__ = "Massimiliano Volpe"
__copyright__ = "Copyright 2022, Massimiliano Volpe"
__email__ = "massimiliano.volpe@liu.se"
__license__ = "GPL-3"


rule variant_annotator:
    input:
        vcf="pgx/variant_filtration/{sample}_{type}.filtered.vcf",
        aln="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("pgx/variant_annotator/{sample}_{type}.annotated.vcf"),
    params:
        dbsnp=config.get("reference", {}).get("dbsnp", ""),
        extra=config.get("variant_annotator", {}).get("extra", ""),  # optional "--resource-allele-concordance -A Coverage --expression db.END",
        java_opts=config.get("variant_annotator", {}).get("java_opts", ""),  # optional extra java options
    log:
        "pgx/variant_annotator/{sample}_{type}.annotated.vcf.log",
    benchmark:
        repeat(
            "pgx/variant_annotator/{sample}_{type}.annotated.vcf.benchmark.tsv",
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
    message:
        "{rule}: annotate vcf on {input}"
    shell:
        "gatk VariantAnnotator"
        " --variant {input.vcf}"
        " --input {input.aln}"
        " --reference {input.ref}"
        " --dbsnp {params.dbsnp}"
        " {params.extra}"
        " --output {output.vcf}"
        " >& {log}"
