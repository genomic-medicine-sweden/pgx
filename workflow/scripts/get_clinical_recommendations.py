#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:05:11 2023

"""

import pandas as pd
import numpy as np


def risk_duplication(x):
    min_dup_var = min(abs(np.array([1 / 4, 1 / 3, 2 / 3, 3 / 4]) - x))
    min_diploid = min(abs(np.array([0, 1 / 2, 1]) - x))
    if (min_dup_var > min_diploid):
        return False
    return True


def get_haplotypes(id_, haplotype_definitions):
    haplo = haplotype_definitions[haplotype_definitions['ID'] ==
                                  id_]["HAPLOTYPE"]
    return ('/'.join(haplo))


def get_recommendations(found_variants, haplotype_definitions,
                        clinical_guidelines_p, interaction_guidelines_p,
                        report, analysed_variants):
    detected_variants = pd.read_csv(found_variants, sep='\t')
    haplotype_definitions = pd.read_csv(haplotype_definitions, sep='\t')
    clinical_guidelines = pd.read_csv(clinical_guidelines_p, sep='\t')
    interaction_guidelines = pd.read_csv(interactions, sep='\t')
    zygosities = ['Hetero-', 'Homo-']

    if detected_variants.empty:
        detected_variants_present = []
        faulty_haplotypes = []
    else:
        detected_variants['Zygosity'] = detected_variants['GT'].map(
            lambda row: zygosities[int(row.split('/')[0])])
        detected_variants['Position'] = detected_variants[
            '#CHROM'] + ':' + detected_variants['POS'].astype(str)
        detected_variants[['Ref.reads', 'Alt.reads'
                           ]] = detected_variants['AD'].str.split(',',
                                                                  expand=True)
        detected_variants['Alt.reads'] = detected_variants['Alt.reads'].astype(
            int)
        detected_variants['Ref.reads'] = detected_variants['Ref.reads'].astype(
            int)
        detected_variants['Haplotype'] = detected_variants['ID'].map(
            lambda x: get_haplotypes(x, haplotype_definitions))
        columns = detected_variants.columns[2:].drop('PL')
        detected_variants_present = detected_variants[columns]
        detected_variants_present[
            'Variantfrekvens'] = detected_variants_present['AF']
        detected_variants_present[
            "Möjlig Duplikation"] = detected_variants_present[
                'Variantfrekvens'].apply(lambda x: risk_duplication(x))
        faulty_haplotypes = pd.Series(
            np.where(~detected_variants_present['Möjlig Duplikation'], '',
                     detected_variants_present['Haplotype']))
        faulty_haplotypes = faulty_haplotypes.map(
            lambda row: row.split("/")).explode().unique()
        order_columns = [
            "GENE", "ID", "Haplotype", "Position", "Zygosity",
            "Variantfrekvens", "Möjlig Duplikation"
        ]
        verbose_columns = [
            "Gen", "rsID", "Möjliga Haplotyper", "Position", "Zygositet",
            "Variantfrekvens", "Möjlig Duplikation"
        ]
        detected_variants_present = detected_variants_present[order_columns]
        detected_variants_present.columns = verbose_columns

    clin_columns = [
        "gene", "Haplotype1", "Haplotype2", "Guideline", "Activity"
    ]
    verbose_columns = [
        "Gen", "Haplotyp 1", "Haplotyp 2", "Klinisk Rekommendation",
        "Aktivitet"
    ]
    clinical_guidelines_present = clinical_guidelines[clin_columns]
    clinical_guidelines_present.columns = verbose_columns

    faulty_haplotype_recommendation = 'Ingen rekommendation ges pga obalans i heterozygositet'

    clinical_guidelines_present.loc[
        clinical_guidelines_present['Haplotyp 1'].isin(faulty_haplotypes),
        'Klinisk Rekommendation'] = faulty_haplotype_recommendation
    clinical_guidelines_present.loc[
        clinical_guidelines_present['Haplotyp 2'].isin(faulty_haplotypes),
        'Klinisk Rekommendation'] = faulty_haplotype_recommendation

    for haplotype in faulty_haplotypes:
        if haplotype == "":
            continue
        interaction_guidelines.loc[
            interaction_guidelines['haplotypes'].str.contains(haplotype),
            'Guideline'] = faulty_haplotype_recommendation

    clinical_guidelines_present[
        'Klinisk Rekommendation'] = clinical_guidelines_present[
            'Klinisk Rekommendation'].str.replace(r'<[^>]+>', '', regex=True)
    interaction_guidelines['Guideline'] = interaction_guidelines[
        'Guideline'].str.replace(r'<[^>]+>', '', regex=True)

    with open(analysed_variants, 'r') as reader:
        analysed = reader.read()

    gene1_interaction = interaction_guidelines['gene1'].values[0]
    gene2_interaction = interaction_guidelines['gene2'].values[0]
    genotype_interaction = interaction_guidelines['haplotypes'].values[0]
    recommendation_interaction = interaction_guidelines['Guideline'].values[0]

    genes = clinical_guidelines_present['Gen'].values.tolist()

    with open(report, 'w') as writer:
        for gene in genes:
            if gene not in interaction_guidelines.values:
                rec_per_gene = clinical_guidelines_present.loc[clinical_guidelines_present['Gen'] == gene]
                haplotype_1 = rec_per_gene['Haplotyp 1'].values[0]
                haplotype_2 = rec_per_gene['Haplotyp 2'].values[0]
                genotype = f"{haplotype_1}/{haplotype_2}"
                recommendation = clinical_guidelines_present.loc[
                    clinical_guidelines_present['Gen'] ==
                    gene]['Klinisk Rekommendation'].values[0]

                writer.write(f"{gene} \nGenotype: {genotype}")
                writer.write(f"\nKlinisk rekommendation: {recommendation}")

        writer.write(f"\n\n{gene1_interaction}-{gene2_interaction}")
        writer.write(f"\nGenotype: {genotype_interaction}")
        writer.write(f"\nKlinisk rekommendation: {recommendation_interaction}\n")
        writer.write(f"\n{analysed}")


if __name__ == "__main__":
    found_variants = snakemake.input.found_variants
    missed_variants = snakemake.input.missed_variants
    clinical_guidelines = snakemake.input.clinical_guidelines
    interactions = snakemake.input.interactions
    haplotype_definitions = snakemake.params.haplotype_definitions
    analysed_variants = snakemake.params.analysed_variants
    recommendations = snakemake.output.recommendations

    get_recommendations(found_variants, haplotype_definitions,
                        clinical_guidelines, interactions, recommendations,
                        analysed_variants)
