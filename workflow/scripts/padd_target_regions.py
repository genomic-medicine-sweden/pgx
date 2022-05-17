# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
"""
Add padding to regions defined in a bed file and switch start coordinate to be first in each row
"""

__author__ = "Joel Ås & Lauri Mesilaakso"
__copyright__ = "Copyright 2022, Joel Ås & Lauri Mesilaakso"
__email__ = "lauri.mesilaakso@regionostergotland.se"
__license__ = "GPL-3"

from pathlib import Path
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
        end, start = start, end
    return [chrom, start, end, region_id]


def add_padding(target_bed: Path, padding: int) -> pd.DataFrame:
    """
    Add padding to target start and end coordinates

    Args:
        target_bed (Path): Path to input bed file
        padding (int): The number bases of padding to add to start and end coordinates

    Returns:
        pd.DataFrame: Pandas dataframe with padded start and end coordinates
    """
    targets : pd.DataFrame = pd.read_csv(target_bed,
                                        sep="\t",
                                        names=["CHROM", "START", "END", "ID"],
                                        dtype={"START": int, "END": int})
    rows : list = []
    for _, row in targets.iterrows():
        chrom, start, end, region_id = switch_coordinates_if_reversed(row[0:4])
        start -= padding
        end += padding
        rows.append([chrom,start,end,region_id])
    return pd.DataFrame(rows, columns=["CHROM", "START", "END", "ID"])


def write_output(input_table: pd.DataFrame, output_f: Path) -> None:
    """
    Write pandas df table into a bed file or in "chr{chrom}:{start}-{end}" format
    Args:
        input_table (pd.DataFrame): Pandas data frame with columns:
                                    ["CHROM", "START", "END", "ID"]
        output_f (Path): Path to output file
    """

    with open(output_f, "w", encoding="utf-8") as f_out:
        for _, row in input_table.iterrows():
            chrom, start, end, region_id = row[0:4]
            f_out.write(f"chr{chrom}\t{start}\t{end}\t{region_id}\n")


def main() -> None:
    """
    Run the main program
    """
    target_bed       = snakemake.input.target_regions
    output_file_name = snakemake.output.padded_target_regions
    padding          = snakemake.params.padding

    padded_coordinates = add_padding(Path(target_bed), int(padding))
    write_output(padded_coordinates, Path(output_file_name))

if __name__ == '__main__':
    main()
