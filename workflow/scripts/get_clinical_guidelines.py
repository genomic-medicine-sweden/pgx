import pandas as pd
import re
import numpy as np


class ArrangeHaplotype:
    """
    Get possible haplotypes from variant combinations detected in file.
    Add clinical guidelines based on the haplotypes detected
    """

    def __init__(self, detected_variants, haplotype_definitions,
                 activity_scores, clinical_guidelines):

        self.detected_variants = pd.read_csv(detected_variants, sep="\t")
        self.haplotype_definitions = pd.read_csv(haplotype_definitions,
                                                 sep="\t")
        self.activity_scores = pd.read_csv(activity_scores, sep="\t")
        self.clinical_guidelines = pd.read_csv(clinical_guidelines, sep="\t")
        if not self.detected_variants.empty:
            self.merge_data()

    def merge_data(self):
        """
        Join possible haplotypes containing ids
        :return:
        """

        self.detected_variants["multival_haplotype"] = \
            self.detected_variants.ID.apply(
                lambda x: list(self.haplotype_definitions["HAPLOTYPE"][
                                   self.haplotype_definitions.ID == x
                                   ])
            )
        self.detected_variants["CN"] = self.detected_variants["GT"].apply(
            lambda x: sum(map(int, re.split('[/|]+', x))))

    def get_haplotypes(self):
        """
        Return tree-structure of possible combinations of haplotypes explaning seen variants.
        Assumptions: Haplotype explaning most variants goes first, if any of variants in haplotype
         is zero the all haplotypes containing that variant is removed for futher chocies.
        :return: Gene - haplotype-tree dict
        """

        def remove_hap(x, y):
            return x.remove(y) if y in x else x

        def _get_haplotypes(variant_subdf, current_haplotype, depth=2):
            idx = variant_subdf["multival_haplotype"].apply(
                lambda x: current_haplotype in x)
            variant_subdf.loc[idx, "CN"] -= 1

            if any(variant_subdf["CN"] == 0):
                variant_subdf["multival_haplotype"] = variant_subdf[
                    "multival_haplotype"].apply(
                        lambda x: remove_hap(x, current_haplotype))

            variant_subdf = variant_subdf[variant_subdf["CN"] != 0]

            if depth == 1:
                if len(variant_subdf) == 0:
                    return [current_haplotype, True]
                else:
                    return [current_haplotype, False]

            if len(variant_subdf) == 0 or not any(
                    variant_subdf["multival_haplotype"].apply(
                        lambda x: bool(x))):
                wt_haplotype = "WT"
                return [
                    current_haplotype,
                    [
                        _get_haplotypes(variant_subdf.copy(), wt_haplotype,
                                        depth - 1)
                    ]
                ]

            remaining_haplo = set([
                element for haplolist in variant_subdf["multival_haplotype"]
                for element in haplolist
            ])
            return [
                current_haplotype,
                [
                    _get_haplotypes(variant_subdf.copy(), hap, depth - 1)
                    for hap in remaining_haplo
                ]
            ]

        genes = set(self.detected_variants["GENE"])
        full_mat = {}
        for gene in genes:
            gene_subset = self.detected_variants[self.detected_variants.GENE ==
                                                 gene]
            candidate_haplotypes = np.array(
                list(
                    set([
                        element
                        for haplolist in gene_subset["multival_haplotype"]
                        for element in haplolist
                    ])))

            order = list(
                reversed(
                    np.argsort([
                        sum(gene_subset["multival_haplotype"].apply(
                            lambda y: x in y)) for x in candidate_haplotypes
                    ])))
            candidate_haplotypes = candidate_haplotypes[order]

            gene_haps = []
            for current_haplotype in candidate_haplotypes:
                gene_haps.append(
                    _get_haplotypes(gene_subset.copy(), current_haplotype))
                idx = gene_subset["multival_haplotype"].apply(
                    lambda x: current_haplotype in x)
                cn = gene_subset.loc[idx, "CN"]
                if any([(c - 1) == 0 for c in cn]):
                    gene_subset["multival_haplotype"] = gene_subset[
                        "multival_haplotype"].apply(
                            lambda x: remove_hap(x, current_haplotype))

            full_mat.update({gene: gene_haps})

        return full_mat

    def get_haplotype_dataframe(self):  # Wow what a mess
        hap_mat = self.get_haplotypes()

        def _haplot_to_row(hap, gene):

            def prim_haplot_to_row(hap, gene):
                current_hap = (f"{gene}-1" if hap[0] == "WT" else hap[0])
                if not type(hap[1]) is list:
                    return [current_hap, gene, hap[1]]
                else:
                    next_hap = [
                        prim_haplot_to_row(c_hap, gene) for c_hap in hap[1]
                    ]
                    return [[current_hap] + c_hap for c_hap in next_hap][0]

            return [prim_haplot_to_row(c_hap, gene) for c_hap in hap]

        out = []
        for gene, hap in hap_mat.items():
            if len(hap) != 0:
                out += _haplot_to_row(hap, gene)

        hap_df = pd.DataFrame(
            out, columns=["Haplotype1", "Haplotype2", "gene", "pass"])
        hap_df = hap_df[hap_df["pass"]]

        return hap_df[["gene", "Haplotype1", "Haplotype2"]]

    def get_clinical_guidelines_table(self):
        if self.detected_variants.empty:
            columns = [
                "gene", "Haplotype1", "Haplotype2", "HAPLOTYPE1",
                "ACTIVITY_SCORE1", "HAPLOTYPE2", "ACTIVITY_SCORE2",
                "Genotype_activity", "Gene", "Activity", "Guideline"
            ]
            return pd.DataFrame(columns=columns)
        hap_df = self.get_haplotype_dataframe()
        hap_df = hap_df.merge(self.activity_scores,
                              how="left",
                              left_on="Haplotype1",
                              right_on="HAPLOTYPE")
        hap_df = hap_df.merge(self.activity_scores,
                              how="left",
                              left_on="Haplotype2",
                              right_on="HAPLOTYPE",
                              suffixes=("1", "2"))
        hap_df["Genotype_activity"] = hap_df["ACTIVITY_SCORE1"] + hap_df[
            "ACTIVITY_SCORE2"]

        hap_df = hap_df.merge(self.clinical_guidelines,
                              how="left",
                              left_on=["gene", "Genotype_activity"],
                              right_on=["Gene", "Activity"])

        return hap_df

    def get_wildtypes(self, hap_df):
        hap_genes = list(hap_df.gene.values)
        for gene in set(self.clinical_guidelines.Gene):
            if hap_df.empty or gene not in hap_genes:
                gene_df = pd.DataFrame({
                    "gene": [gene],
                    "Haplotype1": [gene + "-1"],
                    "Haplotype2": [gene + "-1"],
                    "HAPLOTYPE1": [gene + "-1"],
                    "ACTIVITY_SCORE1": [1],
                    "HAPLOTYPE2": [gene + "-1"],
                    "ACTIVITY_SCORE2": [1],
                    "Genotype_activity": [2.0]
                })
                gene_df = gene_df.merge(self.clinical_guidelines,
                                        how="left",
                                        left_on=["gene", "Genotype_activity"],
                                        right_on=["Gene", "Activity"])
                hap_df = pd.concat([hap_df, gene_df], ignore_index=True)
        return hap_df


def main():
    haplotype_definitions = snakemake.params["haplotype_definitions"]
    haplotype_activity = snakemake.params["haplotype_activity"]
    hidden_haplotypes = snakemake.params["hidden_haplotypes"]
    clinical_guidelines = snakemake.params["clinical_guidelines"]
    variant_csv = snakemake.input["found_variants"]
    output = snakemake.output["csv"]

    ah = ArrangeHaplotype(variant_csv, haplotype_definitions,
                          haplotype_activity, clinical_guidelines)

    df = ah.get_clinical_guidelines_table()
    df = ah.get_wildtypes(df)
    columns = [
        "gene", "Haplotype1", "Haplotype2", "HAPLOTYPE1", "ACTIVITY_SCORE1",
        "HAPLOTYPE2", "ACTIVITY_SCORE2", "Genotype_activity", "Gene",
        "Activity", "Guideline"
    ]
    if not df.empty:
        hidden_haplotypes = pd.read_csv(hidden_haplotypes, sep="\t")
        hidden_haplotypes["comb"] = hidden_haplotypes[[
            "Haplotype1", "Haplotype2"
        ]].apply(lambda x: "".join(sorted(x.tolist())), axis=1)
        df["comb"] = df[["Haplotype1", "Haplotype2"
                         ]].apply(lambda x: "".join(sorted(x.tolist())),
                                  axis=1)
        print(df["comb"])
        print(hidden_haplotypes["comb"])
        df = df[~df["comb"].isin(hidden_haplotypes["comb"])]
    df.to_csv(output, sep="\t", index=False, columns=columns)


if __name__ == '__main__':
    main()
