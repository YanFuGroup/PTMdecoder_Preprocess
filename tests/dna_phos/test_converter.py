import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "external" / "PTMdecoder" / "utils" / "convert_MaxQuant_to_PTMdecoder.py"
FIXTURE = Path(__file__).with_name("fictional_msms_split.txt")


class TestMaxQuantConversion(unittest.TestCase):
    def convert(self, directory: Path, strict: bool = False):
        pep_spec = directory / "pepSpecFile.txt"
        filtered = directory / "filtered_result.txt"
        command = [
            sys.executable,
            str(SCRIPT),
            "-i",
            str(FIXTURE),
            "-o",
            str(pep_spec),
            "-r",
            str(filtered),
            "-s",
            "_HCDFT",
        ]
        if strict:
            command.append("--strict-con-filter")
        subprocess.run(command, check=True, capture_output=True, text=True)
        with filtered.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        return pep_spec.read_text(encoding="utf-8"), rows

    def test_default_outputs_and_filters(self):
        with tempfile.TemporaryDirectory() as tmp:
            pep_spec, rows = self.convert(Path(tmp))

        self.assertEqual(
            pep_spec,
            "AS\nTest_A_HCDFT.mgf\t101\nTY\nTest_B_HCDFT.mgf\t105\n",
        )
        self.assertEqual([row["Scan"] for row in rows], ["101", "102", "105"])
        self.assertEqual(rows[0]["protein"], "P001")
        self.assertEqual(rows[0]["modification"], "Phospho (STY)")
        self.assertEqual(rows[0]["modificationlocation"], "2")
        self.assertEqual(rows[0]["precursor_neutral_mass"], "200.000000")
        self.assertEqual(rows[1]["modification"], "-")
        self.assertEqual(rows[2]["massdiff"], "0")
        self.assertEqual(rows[2]["num_match_ions"], "0")

    def test_strict_contaminant_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            pep_spec, rows = self.convert(Path(tmp), strict=True)

        self.assertEqual(pep_spec, "TY\nTest_B_HCDFT.mgf\t105\n")
        self.assertEqual([row["Scan"] for row in rows], ["102", "105"])


if __name__ == "__main__":
    unittest.main()
