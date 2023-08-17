import unittest
from dataclasses import dataclass
from pathlib import Path
import pandas as pd
from tempfile import NamedTemporaryFile
import csv
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


class TestReadDf(unittest.TestCase):

    def test_read_df(self):

        @dataclass
        class TestCase:
            name: str
            input: Path
            expected: pd.DataFrame

        with NamedTemporaryFile(mode="w") as tmp_f:
            # Generate a bed file with known contents
            tsv_writer = csv.writer(tmp_f, delimiter="\t")
            tsv_writer.writerow(["chr1", "200", "300", "ID1"])
            tmp_f.seek(0)

            # Generate a pandas df with expected contents
            padded_bed = pd.DataFrame([["chr1", "200", "300", "ID1"]],
                                      columns=["CHROM", "START", "END", "ID"])
            padded_bed = padded_bed.astype({"START": int, "END": int})

            testcases = [
                TestCase(
                    name="Read a bed file successfully",
                    input=Path(tmp_f.name),
                    expected=padded_bed,
                ),
            ]

            for case in testcases:
                actual = ptr.read_df(case.input)
                # breakpoint()
                self.assertTrue(
                    pd.testing.assert_frame_equal(actual, case.expected) is
                    None,
                    "failed test '{}': expected {}, got {}".format(
                        case.name, case.expected, actual),
                )


class TestAddPadding(unittest.TestCase):

    def test_add_padding(self):

        @dataclass
        class TestCase:
            name: str
            input: pd.DataFrame
            padding: int
            expected: pd.DataFrame

        unpadded_bed = pd.DataFrame([["chr1", "200", "300", "ID1"]],
                                    columns=["CHROM", "START", "END", "ID"])
        unpadded_bed = unpadded_bed.astype({"START": int, "END": int})

        unpadded_bed_chr_missing = pd.DataFrame(
            [["1", "200", "300", "ID1"]],
            columns=["CHROM", "START", "END", "ID"])
        unpadded_bed_chr_missing = unpadded_bed.astype({
            "START": int,
            "END": int
        })

        padded_bed = pd.DataFrame([["chr1", "100", "400", "ID1"]],
                                  columns=["CHROM", "START", "END", "ID"])
        padded_bed = padded_bed.astype({"START": int, "END": int})

        testcases = [
            TestCase(
                name="Pad coordinates successfully",
                input=unpadded_bed,
                padding=100,
                expected=padded_bed,
            ),
            TestCase(
                name="Pad coordinates successfully",
                input=unpadded_bed_chr_missing,
                padding=100,
                expected=padded_bed,
            ),
        ]

        for case in testcases:
            actual = ptr.add_padding(case.input, case.padding)
            self.assertTrue(
                pd.testing.assert_frame_equal(actual, case.expected) is None,
                "failed test '{}': expected {}, got {}".format(
                    case.name, case.expected, actual),
            )
