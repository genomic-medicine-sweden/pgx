import pandas as pd


def reform(target_bed, output_f, detected_variants, padding, file_format):
    targets = pd.read_csv(target_bed,
                          sep="\t",
                          names=["CHROM", "START", "END", "ID", "GENE"],
                          dtype={
                              "START": int,
                              "END": int
                          })
    if detected_variants is not None:
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
                f.write(f"{chrom}\t{start}\t{end}\t{id}\n")

            else:
                f.write(f"{chrom}:{start}-{end}\n") # noqa


def main():
    target_bed = snakemake.input["target_bed"]
    output_file = snakemake.output["interval"]
    detected_variants = snakemake.input.get("detected_variants", None)
    padding = snakemake.params["padding"]
    file_format = snakemake.params["file_format"]
    reform(target_bed, output_file, detected_variants, int(padding),
           file_format)


if __name__ == '__main__':
    main()
