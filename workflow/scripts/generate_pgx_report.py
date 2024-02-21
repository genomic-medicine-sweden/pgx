#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 11:34:29 2023

"""

import pandas as pd
import jinja2


def highlight_greaterthan(s, threshold, column):
    is_max = pd.Series(data=False, index=s.index)
    is_max[column] = s.loc[column] < threshold
    return [
        'background-color: red ; font-weight: bold ; color: white'
        if is_max.any() else '' for v in is_max
    ]


def highlight_risk_genotype(s, column):
    is_max = pd.Series(data=False, index=s.index)
    is_max[column] = s[column].str.contains('-[^1]') is not True
    return [
        'background-color: red ; font-weight: bold ; color: white'
        if not is_max.any() else '' for v in is_max
    ]


def create_depth_table(missing, haplotype_definitions, threshold):
    missing = missing[['ID', 'Locus', 'Total_Depth']]
    depth_table = missing.merge(haplotype_definitions, on=['ID'])
    depth_table.rename(columns={
        'ID': 'rsID',
        'HAPLOTYPE': 'Haplotyp',
        'Locus': 'Position',
        'Total_Depth': 'Läsdjup'
    },
                       inplace=True)

    new_col_order = ['rsID', 'Haplotyp', 'Position', 'Läsdjup']
    depth_table = depth_table[new_col_order]

    depth_table = depth_table.groupby(['Position', 'Läsdjup',
                                       'rsID'])['Haplotyp'].apply(
                                           '/'.join).reset_index()
    depth_table['sort'] = depth_table['Position'].str.extract(
        r'(\d+)', expand=False).astype(int)
    depth_table = depth_table.sort_values('sort')
    depth_table = depth_table.drop('sort', axis=1)

    if depth_table['Läsdjup'].le(threshold).any():
        warning_depth = "!"
    else:
        warning_depth = ""

    styled_depth_table = depth_table.style.apply(highlight_greaterthan,
                                                 threshold=threshold,
                                                 column=['Läsdjup'],
                                                 axis=1).hide(axis="index")
    return styled_depth_table, warning_depth


def create_styled_clinical_recommendations(rekommendation):

    recommendation_restructured = []

    for element in rekommendation.split('<br><br>'):
        if len(element) > 0:
            split_recommendation = element.split('<br>')
            gene = split_recommendation[0]
            genotype = split_recommendation[1].split(' ')[1]
            rec = split_recommendation[2]
            recommendation_restructured.append([gene, genotype, rec])

    pandas_recc = pd.DataFrame(
        recommendation_restructured,
        columns=['Gen', 'Genotyp', 'Klinisk rekommendation'])

    if pandas_recc['Genotyp'].str.contains('-[^1]').any() is not True:
        warning_rec = ""
    else:
        warning_rec = "!"

    styled_recc = pandas_recc.style.apply(highlight_risk_genotype,
                                          column=['Genotyp'],
                                          axis=1).hide(axis="index")

    return styled_recc, warning_rec


def generate_html_report(styled_depth_table, missing, template, dbsnp,
                         ref_path, rekommendation, analysed_variants,
                         detected_variants, haplotype_definitions,
                         warning_depth, warning_rec):
    sample_column = missing.columns[3]
    sample = sample_column.split('_')[2]
    sample_id = sample_column.split('_')[3]

    with open(template, 'r') as f:
        template = jinja2.Template(source=f.read())
    styled_depth_report = template.render(
        sample=sample,
        sample_id=sample_id,
        dbsnp=dbsnp,
        ref_path=ref_path,
        rekommendation=rekommendation.to_html(index=False),
        analysed_variants=analysed_variants,
        detected_variants=detected_variants.to_html(index=False),
        haplotype_definitions=haplotype_definitions.to_html(index=False),
        table=styled_depth_table.set_uuid('läsdjup').to_html(index=False),
        warning_depth=warning_depth,
        warning_rec=warning_rec)

    return styled_depth_report


if __name__ == "__main__":
    missed_variants = snakemake.input.missed_variants
    detected_variants = snakemake.input.detected_variants
    dbsnp = snakemake.params.dbsnp
    ref = snakemake.params.ref
    haplotype_definitions = snakemake.params.haplotype_definitions
    read_depth = snakemake.params.read_depth
    html_template = snakemake.params.html_template
    depth_table = snakemake.output.depth_table
    html_report = snakemake.output.html_report
    recommendations = snakemake.input.recommendations

    haplotype_definitions = pd.read_csv(haplotype_definitions, sep='\t')
    missed_variants = pd.read_csv(missed_variants, sep='\t')
    detected_variants = pd.read_csv(detected_variants, sep='\t')

    rekommendation = []
    with open(recommendations) as f:
        rec = f.read().split('Analyserade varianter:')

    rekommendation = rec[0].replace('\n', '<br>')
    analysed_variants = rec[1].replace('\n', '<br>')

    clinical_recommendations = create_styled_clinical_recommendations(
        rekommendation)

    styled_rec = clinical_recommendations[0]
    warning_rec = clinical_recommendations[1]

    unstyled_depth_table = create_depth_table(missed_variants,
                                              haplotype_definitions,
                                              threshold=read_depth)

    styled_depth_table = unstyled_depth_table[0]
    warning_depth = unstyled_depth_table[1]
    styled_depth_report = generate_html_report(
        styled_depth_table, missed_variants, html_template, dbsnp, ref,
        styled_rec, analysed_variants, detected_variants,
        haplotype_definitions, warning_depth, warning_rec)

    with open(html_report, 'w') as w:
        w.write(styled_depth_report)

    styled_depth_table.to_excel(depth_table, index=False)
