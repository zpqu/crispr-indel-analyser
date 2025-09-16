# Copyright (C) 2025 Zhipeng Qu
#
# Permission is hereby granted, free of charge, to any person obtaining a copy # of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Unit tests for IndelAnalyser."""

import pytest
from unittest.mock import patch
from pathlib import Path
import json
import pandas as pd
from src.analysis.indel_analyser import IndelAnalyser


class TestIndelAnalyser:
    """Test suite for IndelAnalyser class."""

    def test_init_paths(self, indel_analyser, temp_demux_dir, result_dir):
        """Test that paths are correctly set during initialisation."""
        assert indel_analyser.demux_dir == temp_demux_dir
        assert indel_analyser.result_dir == result_dir


    def test_match_flanks_exact_match(self, indel_analyser):
        """Test _match_flanks finds exact flank sequences."""
        seq = "TTAAACCCGGGAA"
        matched, up_end, down_start = indel_analyser._match_flanks(seq, "TTA", "GAA")
        assert matched
        assert up_end == 3
        assert down_start == 10

    def test_match_flanks_with_mismatch(self, indel_analyser):
        """Test flank matching with one mismatch."""
        indel_analyser.mismatch = 1
        seq = "TAAAACCCGGGAA"
        matched, up_end, down_start = indel_analyser._match_flanks(seq, "TTA", "GAT")
        assert matched
        assert up_end == 3
        assert down_start == 10

    def test_analyse_no_edit_txt_output(self, indel_analyser, result_dir):
        """Test no-edit detection with TXT output."""
        result = indel_analyser.analyse("sample1", output_format="txt")

        assert result is not None
        assert result["num_other"] == 1
        assert result["num_ins"] == 0
        assert result["num_del"] == 0
        assert result["per_ins"] == 0
        assert result["per_del"] == 0

        output_file = result_dir / "sample1_summary.txt"
        assert output_file.exists()
        content = output_file.read_text()
        assert "sample1\t0\t0\t1\t0" in content

    def test_analyse_insertion_detection(self, indel_analyser, result_dir):
        """Test insertion is detected and recorded."""
        result = indel_analyser.analyse("sample1_ins", output_format="json")

        assert result is not None
        assert result["num_other"] == 0
        assert result["num_ins"] > 0
        assert result["num_del"] == 0
        assert result["per_ins"] > 0
        assert result["per_del"] == 0
        assert result["pos_ins"]

    def test_analyse_deletion_detection(self, indel_analyser, result_dir):
        """Test deletion is detected and recorded."""
        result = indel_analyser.analyse("sample1_del", output_format="json")

        assert result is not None
        assert result["num_other"] == 0
        assert result["num_ins"] == 0
        assert result["num_del"] > 0
        assert result["per_ins"] == 0
        assert result["per_del"] > 0
        assert result["pos_del"]

    def test_analyse_json_output_structure(self, indel_analyser, result_dir):
        """Test JSON output contains structured data."""
        indel_analyser.analyse("sample1_ins", output_format="json")

        json_file = result_dir / "sample1_ins_summary.json"
        assert json_file.exists()

        with open(json_file) as f:
            data = json.load(f)
            assert "sample" in data
            assert "pos_ins" in data
            assert isinstance(data["pos_ins"], dict)

    def test_analyse_missing_sample_returns_none(
            self, 
            indel_analyser, 
            processed_meta_data, 
            tmp_path, 
            result_dir, 
        ):
        """Test missing FASTQ is handled gracefully."""
        expected_path = tmp_path / "missing" / "sample1.fq.gz"
        analyser = IndelAnalyser(
            meta_data=processed_meta_data, 
            demux_dir=tmp_path / "missing", 
            result_dir=result_dir, 
        )
        with patch("builtins.print") as mock_print:
            result = analyser.analyse("sample1")
            mock_print.assert_any_call(
                f"[WARNING] FASTQ not found: {expected_path}")
            assert result is None

    def test_summarise_all_generates_summary_csv(self, indel_analyser, result_dir):
        """Test summarise_all creates combined CSV."""
        indel_analyser.analyse("sample1")
        indel_analyser.analyse("sample1_ins")

        summary_csv = result_dir / "summary.csv"
        indel_analyser.summarise_all(summary_csv)

        assert summary_csv.exists()
        df = pd.read_csv(summary_csv)
        assert len(df) == 2
        assert "sample" in df.columns
        assert "num_ins" in df.columns

