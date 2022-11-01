from pysam import VariantFile


def filter_variants(vcf, read_ratio, depth, output):
    """
    Soft filter all variants with suspicious read ratio and insufficient read-depth
    """
    vcf_in = VariantFile(vcf)
    new_header = vcf_in.header
    new_header.filters.add(f"AR{read_ratio}", None, None,
                           f"Ratio of ref/alt reads lower than {read_ratio}")
    new_header.filters.add(f"DP{depth}", None, None,
                           f"DP is lower than {depth}x")
    vcf_out = VariantFile(output, "w", header=new_header)

    for record in vcf_in.fetch():
        ad = record.samples[0]["AD"]
        #  No multiallelic split

        if record.info["DP"] < depth:
            record.filter.add("DP100")

        elif len(ad) == 2:
            n_ref, n_alt = ad
            if n_alt / (n_ref + n_alt) < read_ratio:
                record.filter.add(f"AR{read_ratio}")
            else:
                record.filter.add("PASS")
        vcf_out.write(record)


def main():
    filter_variants(snakemake.input["vcf"], snakemake.params.read_ratio,
                    snakemake.params.read_depth,
                    snakemake.output["filtered_vcf"])


if __name__ == '__main__':
    main()
