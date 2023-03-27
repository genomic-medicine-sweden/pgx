__author__ = "Chelsea Ramsin"
__copyright__ = "Copyright 2023, Chelsea Ramsin"
__email__ = "Chelsea.Ramsin@regionostergotland.se"
__license__ = "GPL-3"


rule generate_pgx_report:
    input:
        found_variants="pgx/detected_variants/{sample}_{type}.annotated.csv",
        missed_variants = "pgx/append_id_to_gdf/{sample}_{type}.depth_at_missing_annotated.gdf",
        clinical_guidelines="pgx/get_clinical_guidelines/{sample}_{type}.output.csv",
        interactions="pgx/get_interaction_guidelines/{sample}_{type}.output.csv",
        depth_at_baits="pgx/depth_of_coverage/depth_of_baits/{sample}_{type}.output.gdf",
        report_template="workflow/scripts/generate_pgx_report.Rmd",
    output:
        html="pgx/generate_pgx_report/{sample}_{type}_pgx_report.html",
    params:
    	haplotype_definitions=config.get("get_clinical_guidelines", {}).get("haplotype_definitions", ""),
    	dbsnp=config.get("reference", {}).get("dbsnp", ""),
        ref=config.get("reference", {}).get("fasta", ""),
        name=config.get("generate_pgx_report", {}).get("name", ""),
        address=config.get("generate_pgx_report", {}).get("address", ""),
        email=config.get("generate_pgx_report", {}).get("email", ""),
        phone_number=config.get("generate_pgx_report", {}).get("phone_number", ""),
        design_bed=config.get("reference", {}).get("design_bed",""),
    log:
        "pgx/generate_pgx_report/{sample}_{type}_pgx_report.html.output.log",
    benchmark:
        repeat(
            "pgx/generate_pgx_report/{sample}_{type}_pgx_report.html.benchmark.tsv",
            config.get("generate_pgx_report", {}).get("benchmark_repeats", 1)
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
    shell:
        """
        WKDIR=$(pwd)  # Needed since Rscript will set wd to location of file not session
        echo "the working dir is: $WKDIR"
        INTDIR=$(echo {output.html} | head -c -6)
        DBSNP=$(basename {params.dbsnp})

        Rscript -e "\
        library(rmdformats); \
        rmarkdown::render(\
        input = '{input.report_template}', \
        params = \
        list(\
        title  = 'Farmakogenomisk analys av {wildcards.sample}', \
        found_variants = '$WKDIR/{input.found_variants}', \
        missed_variants = '$WKDIR/{input.missed_variants}', \
        haplotype_definitions = '$WKDIR/{params.haplotype_definitions}', \
        clinical_guidelines = '$WKDIR/{input.clinical_guidelines}', \
        interaction_guidelines = '$WKDIR/{input.interactions}', \
        data_location = '$WKDIR/data', \
        depth_file = '$WKDIR/{input.depth_at_baits}', \
        sample = '{wildcards.sample}', \
        dbsnp = '$DBSNP', \
        ref = '{params.ref}', \
        name = '{params.name}', \
        adress = '{params.address}', \
        mail = '{params.email}', \
        phone = '{params.phone_number}'\
        ), \
        output_file = '$WKDIR/{output.html}', \
        output_format = c('readthedown'))" &> {log}
        """
