import pandas as pd
import numpy as np
from itertools import product


class GetInteractions:

    def __init__(self, variants_file, guideline_file):
        self.variants_df = pd.read_csv(variants_file, sep="\t")
        self.variants_df = self.variants_df.sort_values(by=["gene"])
        self.interaction_guidelines_df = pd.read_csv(guideline_file, sep="\t")

    def get_possible_interactions(self):
        # TODO: refactor inefficient and ugly
        interacting_genes = set(
            self.interaction_guidelines_df[["gene1", "gene2"]].values.flat)
        if self.variants_df.empty:
            return self.interaction_guidelines_df[np.zeros(len(
                self.interaction_guidelines_df),
                                                           dtype=bool)]
        self.variants_df = self.variants_df[self.variants_df.gene.isin(
            interacting_genes)]

        if self.variants_df.empty:
            return self.interaction_guidelines_df[np.zeros(len(
                self.interaction_guidelines_df),
                                                           dtype=bool)]

        index_combinations = product(range(len(self.variants_df)),
                                     range(len(self.variants_df)))

        prev_idx = np.zeros(len(self.interaction_guidelines_df), dtype=bool)
        haplotypes = [""] * len(self.interaction_guidelines_df)
        for i, j in index_combinations:
            if not i == j:
                gene1 = self.variants_df.iloc[[i]].gene.values[0]
                gene2 = self.variants_df.iloc[[j]].gene.values[0]

                activity1 = int(self.variants_df.iloc[[i]].Genotype_activity)
                activity2 = int(self.variants_df.iloc[[j]].Genotype_activity)

                all_haplotypes = ",".join([
                    self.variants_df.iloc[[i]]["Haplotype1"].values.flat[0],
                    self.variants_df.iloc[[i]]["Haplotype2"].values.flat[0],
                    self.variants_df.iloc[[j]]["Haplotype1"].values.flat[0],
                    self.variants_df.iloc[[j]]["Haplotype2"].values.flat[0]
                ])

                idx = (self.interaction_guidelines_df.gene1 == gene1) & \
                      (self.interaction_guidelines_df.gene2 == gene2) & \
                      (self.interaction_guidelines_df.activity_1 == activity1) & \
                      (self.interaction_guidelines_df.activity_2 == activity2)

                if any(idx):
                    prev_idx += idx
                    for idx_i in idx:
                        if idx_i:
                            haplotypes[idx_i] = all_haplotypes

        haplotypes = [h for h in haplotypes if h != ""]
        interactions = self.interaction_guidelines_df[prev_idx]
        interactions["haplotypes"] = haplotypes
        return interactions

    def run_and_write(self, output):
        interactions = self.get_possible_interactions()
        interactions.to_csv(output, index=None, sep="\t")


def main():
    interaction_guidelines = snakemake.params["interacting_targets"]
    diploids = snakemake.input["diploids"]
    output = snakemake.output["csv"]
    get_interactions = GetInteractions(diploids, interaction_guidelines)
    get_interactions.run_and_write(output)


if __name__ == '__main__':
    main()
