import unittest
from dataclasses import dataclass
from pathlib import Path
import pandas as pd
import tempfile
import padd_target_regions as ptr
from typing import List


class TestPaddTargetRegions(unittest.TestCase):

    def test_switch_coordinates_if_reversed(self):

        @dataclass
        class TestCase:
            name: str
            input: List
            expected: List

        testcases = [
            TestCase(
                name="Switch coordinates successfully",
                input=["chr1", 2, 1, "ID1"],
                expected=["chr1", 1, 2, "ID1"],
            ),
            TestCase(
                name="Don't switch coordinates when it's not necessary",
                input=["chr1", 1, 2, "ID1"],
                expected=["chr1", 1, 2, "ID1"],
            ),
        ]

        for case in testcases:
            actual = ptr.switch_coordinates_if_reversed(case.input)
            self.assertEqual(
                case.expected,
                actual,
                "failed test '{}': expected {}, got {}".format(
                    case.name, case.expected, actual),
            )


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
