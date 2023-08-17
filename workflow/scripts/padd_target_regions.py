# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
"""
Add padding to regions defined in a bed file and switch start coordinate to
be first in each row
"""

__author__ = "Joel Ås & Lauri Mesilaakso"
__copyright__ = "Copyright 2022, Joel Ås & Lauri Mesilaakso"
__email__ = "lauri.mesilaakso@regionostergotland.se"
__license__ = "GPL-3"

from pathlib import Path
import logging
from typing import List
import pandas as pd


def switch_coordinates_if_reversed(coordinates: List) -> List:
    """
    Switch start and end coordinates when they are reversed, i.e. (start > end)
    Args:
        coordinates (List): A list of coordinates and their IDs

    Returns:
        List: A list of coordinates where start < end
    """
    chrom, start, end, region_id = coordinates
    if start > end:
        logging.info(f"Original start coordinate: {str(start)}"
                     f", end coordinate: {str(end)}")
        end, start = start, end
    return [chrom, start, end, region_id]


def read_df(input_bed_file: Path) -> pd.DataFrame:
    """Read an input 4 column bed file into a pandas dataframe

    Args:
        input_bed_file (Path): Path to input 4 column bed file

    Returns:
        pd.DataFrame: Parsed pandas dataframe containing the bed file data
    """
    return pd.read_csv(input_bed_file,
                       sep="\t",
                       names=["CHROM", "START", "END", "ID"],
                       dtype={
                           "START": int,
                           "END": int
                       })


def add_padding(target_bed_df: pd.DataFrame, padding: int) -> pd.DataFrame:
    """
    Add padding to target start and end coordinates

    Args:
        target_bed_df (pd.DataFrame): Bed file as a pandas dataframe
        padding (int): The number bases of padding to add to start
                        and end coordinates

    Returns:
        pd.DataFrame: Pandas dataframe with padded start and end coordinates
    """
    rows: list = []
    for _, row in target_bed_df.iterrows():
        chrom, start, end, region_id = switch_coordinates_if_reversed(row[0:4])
        start -= padding
        end += padding
        if (not str(chrom).startswith('chr')):
            chrom = f"chr{chrom}"
        rows.append([chrom, start, end, region_id])
    return pd.DataFrame(rows, columns=["CHROM", "START", "END", "ID"])


if __name__ == '__main__':
    target_bed = snakemake.input.target_regions
    output_file_name = snakemake.output.padded_target_regions
    padding = snakemake.params.padding

    logging.basicConfig(level=logging.INFO, filename=snakemake.log[0])

    logging.info(f"Read {target_bed}")
    bed_df = read_df(Path(target_bed))
    padded_coordinates = add_padding(bed_df, int(padding))

    padded_coordinates.to_csv(Path(output_file_name),
                              index=False,
                              sep='\t',
                              header=False)
    logging.info(f"Padded bed file written: {output_file_name}")
