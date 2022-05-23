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
