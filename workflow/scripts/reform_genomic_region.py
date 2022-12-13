import pandas as pd
import argparse
import sys


def reform(target_bed, output_f, detected_variants, padding, file_format):
    targets = pd.read_csv(target_bed, sep="\t",
                          names=["CHROM", "START", "END", "ID", "GENE"],
                          dtype={"START": int, "END": int})
    if detected_variants != "":
        detected_rsid = pd.read_csv(detected_variants, sep="\t").ID
        targets = targets[~targets.ID.isin(detected_rsid)]

    bed = file_format == "bed"
    with open(output_f, "w+") as f:
        for i, row in targets.iterrows():
            chrom, start, end, id = row[0:4]
            if start > end:
                end, start = start, end
            start -= padding
            end += padding
            if bed:
                f.write(f"chr{chrom}\t{start}\t{end}\t{id}\n")
            else:
                f.write(f"chr{chrom}:{start}-{end}\n")


def main():
    target_bed = snakemake.params["target_bed"]
    output_file = snakemake.output["interval"]
    detected_variants = snakemake.input["detected_variants"]
    padding = snakemake.params["padding"]
    file_format = snakemake.params["file_format"]
    reform(target_bed, output_file, detected_variants, int(padding), file_format)


if __name__ == '__main__':
    main()
