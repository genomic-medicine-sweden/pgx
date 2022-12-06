import pandas as pd
import gzip
import re
import os


class GDF:
    """
    We are assuming single sample GDF
    """

    def __init__(self, filename):
        self.data = pd.read_csv(filename, sep="\t")
        self.data_cols = self.data.columns.values

        pos_df = pd.DataFrame(
            self.data.Locus.apply(lambda x: x.split(":")).to_list(),
            columns=["CHROM", "POS"])
        pos_df.POS = pos_df.POS.astype('int64')
        self.data = pd.concat([self.data, pos_df], axis=1)

    def rsid_per_position(self, target_bed):

        def _annotate(x, targets):
            try:
                ids = targets[(targets.CHROM == x.CHROM)
                              & (targets.START <= x.POS) &
                              (x.POS <= targets.END)].ID.to_list()
                return ", ".join(ids)
            except IndexError:
                return "-"

        targets = pd.read_csv(target_bed,
                              sep="\t",
                              names=["CHROM", "START", "END", "ID", "GENE"])
        targets["save"] = targets.START
        idx_swap = targets.START > targets.END
        targets.loc[idx_swap, "START"] = targets.loc[idx_swap, "END"]
        targets.loc[idx_swap, "END"] = targets.loc[idx_swap, "save"]
        targets["CHROM"] = targets.CHROM.apply(lambda x: f"chr{x}")
        self.data["ID"] = self.data.apply(lambda x: _annotate(x, targets),
                                          axis=1)

    def write_proccessed_gdf(self, filename, annotate=True):
        if annotate:
            self.data.to_csv(filename, sep="\t", index=False)
        else:
            self.data.to_csv(filename,
                             sep="\t",
                             columns=self.data_cols,
                             index=False)


class VCF:
    """
    We are assuming single sample VCF
    """

    def __init__(self, filename):
        self.meta = []
        self.data = pd.DataFrame()
        self.original_header = []
        self.read_vcf(filename)

    def read_vcf(self, filename):
        if ".gz" in filename:
            f = gzip.open(filename, "rt")
        else:
            f = open(filename, "r")

        lines = f.readlines()
        lines = [lr.strip() for lr in lines]
        i = None
        for i, line in enumerate(lines):
            if re.search("^#CHROM", line):
                break

        if i is None:
            raise ImportError("No lines in: " + filename)

        self.meta = lines[:i - 1]
        data = [lr.split("\t") for lr in lines[i:]]
        self.data = pd.DataFrame(data[1:], columns=data[0])
        self.original_header = self.data.columns.values
        if not self.data.empty:
            sample_column = self.data.columns.values[-1]
            # max_len_idx = self.data.FORMAT.str.len().idxmax
            # format_columns = self.data.FORMAT[max_len_idx].split(":")
            format_columns = self.data.FORMAT[0].split(":")
            format_split = pd.DataFrame(self.data[sample_column].apply(
                lambda x: x.split(":")).to_list(),
                                        columns=format_columns)
            self.data = pd.concat([self.data, format_split], axis=1)

    def filter_snp(self, filter_file, exclude=True, column="SNP"):
        filter_snps = pd.read_csv(filter_file, sep="\t")[column]
        if exclude:
            self.data = self.data[~self.data.ID.isin(filter_snps)]
        else:
            self.data = self.data[self.data.ID.isin(filter_snps)]

    def write_vcf(self, filename_out):
        with open(filename_out, "w+") as f:
            for line in self.meta:
                f.write(line + "\n")

            self.data.to_csv(f,
                             mode="a",
                             sep="\t",
                             columns=self.original_header,
                             index=False)


class VariantQCCollection:

    def __init__(self, target_bed, vcf_filename):
        self.bed_targets = pd.read_csv(target_bed, names=None, sep="\t")
        self.bed_targets.columns = ["#CHROM", "START", "END", "ID", "GENE"]
        self.sample = os.path.basename(vcf_filename)
        self.vcf = VCF(vcf_filename)
        self.detected_variants = []

    def detect_variants(self):
        """
        Select all variants within VCF with ID same as in target_bed
        """
        self.detected_variants = self.vcf.data.ID[
            (self.vcf.data.ID.isin(self.bed_targets.ID))
            & (self.vcf.data.FILTER == "PASS")].tolist()

    def write_detected_variant_qc(self, output_file):
        """
        Write detected variants to csv
        """
        if not self.detected_variants:
            self.detect_variants()

        current_variants = self.vcf.data[self.vcf.data.ID.isin(
            self.detected_variants)]
        current_variants = current_variants.merge(
            self.bed_targets[["ID", "GENE"]], on="ID")

        current_variants.to_csv(output_file, index=False, sep="\t")


def main():
    bed_f = snakemake.params["target_bed"]
    vcf_f = snakemake.input["vcf"]
    output_f = snakemake.output["csv"]
    var_collect = VariantQCCollection(bed_f, vcf_f)
    var_collect.write_detected_variant_qc(output_f)


if __name__ == '__main__':
    main()
