# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
"""
Test functions used for padding bed files
"""

__author__ = "Lauri Mesilaakso"
__copyright__ = "Copyright 2022, Lauri Mesilaakso"
__email__ = "lauri.mesilaakso@regionostergotland.se"
__license__ = "GPL-3"

import pytest
import padd_target_regions as ptr


def test_switch_coordinates_no_switch():
    coordinates = ["chr1", 1, 2, "ID1"]
    assert ptr.switch_coordinates_if_reversed(coordinates) == coordinates


def test_switch_coordinates_switch_needed():
    coordinates = ["chr1", 2, 1, "ID1"]
    assert ptr.switch_coordinates_if_reversed(coordinates) == ["chr1", 1, 2, "ID1"]

@pytest.fixture(scope="session")
def bed_file(tmp_path_factory):
    fn = tmp_path_factory.mktemp("data") / "test.bed"
    content = "chr1\t200\t300\tID1"
    p.write_text(content)
    return fn


def test_add_padding(bed_file):
    bed_fn = bed_file()
    padded_bed = add_padding(bed_fn, 100)
