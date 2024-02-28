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


def create_depth_table(missing, haplotype_definitions):
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

    styled_depth_table = depth_table.style.apply(highlight_greaterthan,
                                                 threshold=100,
                                                 column=['Läsdjup'],
                                                 axis=1).hide(axis="index")

    return styled_depth_table


def generate_html_report(styled_depth_table, missing, template, dbsnp,
                         ref_path):
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
        table=styled_depth_table.to_html(index=False))

    return styled_depth_report


if __name__ == "__main__":
    missed_variants = snakemake.input.missed_variants
    dbsnp = snakemake.params.dbsnp
    ref = snakemake.params.ref
    haplotype_definitions = snakemake.params.haplotype_definitions
    html_template = snakemake.params.html_template
    excel_report = snakemake.output.excel_report
    html_report = snakemake.output.html_report

    haplotype_definitions = pd.read_csv(haplotype_definitions, sep='\t')
    missed_variants = pd.read_csv(missed_variants, sep='\t')

    styled_depth_table = create_depth_table(missed_variants,
                                            haplotype_definitions)
    styled_depth_report = generate_html_report(styled_depth_table,
                                               missed_variants, html_template,
                                               dbsnp, ref)

    with open(html_report, 'w') as w:
        w.write(styled_depth_report)

    styled_depth_table.to_excel(excel_report, index=False)
