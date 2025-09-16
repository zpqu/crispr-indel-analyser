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

"""Unit tests for FASTQDemultiplexer."""

from pathlib import Path
from unittest.mock import patch, Mock
import pytest
from src.demux.fastq_demux import FASTQDemultiplexer


class TestFASTQDemultiplexer:
    """Test suite for FASTQDemultiplexer."""

    def test_init_with_path_string(self, sample_meta_csv):
        """Test initialisation with string path."""
        demux = FASTQDemultiplexer(sample_meta_csv, output_dir="test_output")
        assert demux.output_dir == Path("test_output")

    def test_init_with_pathlib_path(self, sample_meta_csv):
        """Test initialisation with Path object."""
        demux = FASTQDemultiplexer(sample_meta_csv, output_dir=Path("out"))
        assert demux.output_dir == Path("out")

    def test_find_subseq_with_mismatch_exact_match(self, sample_meta_csv):
        """Test _find_subseq_with_mismatch finds exact match."""
        demux = FASTQDemultiplexer(sample_meta_csv, mismatch=0)
        fp = "AGGTCA"
        seq = "NNNAGGTCANNN"
        positions = demux._find_subseq_with_mismatch(seq, fp)
        assert 3 in positions
        assert len(positions) == 1

    def test_find_subseq_with_mismatch_one_mismatch(self, sample_meta_csv):
        """Test _find_subseq_with_mismatch with one mismatch."""
        demux = FASTQDemultiplexer(sample_meta_csv, mismatch=1)
        fp = "ATCCG"
        seq = "NNNATTCGNNN"
        positions = demux._find_subseq_with_mismatch(seq, fp)
        assert 3 in positions
        assert len(positions) == 1

    def test_find_subseq_no_match(self, sample_meta_csv):
        """Test returns empty list when no match."""
        demux = FASTQDemultiplexer(sample_meta_csv, mismatch=1)
        fp = "ATCG"
        seq = "NNNGGGGNNN"
        positions = demux._find_subseq_with_mismatch(seq, "ATCG") 
        assert positions == []

    def test_match_read_forward_orientation(self, sample_meta_csv):
        """Test forward read matching: fp ... rp."""
        demux = FASTQDemultiplexer(
            sample_meta_csv,
            mismatch=0,
            min_dist=5,
            max_dist=100
        )
        # fp + spacer + target + rp
        seq = "AGGTCA" + "NNN" + "AAACCCGGG" + "CTAGCT"
        assert demux._match_read(seq) == "sample1"

    def test_match_read_with_umi(self, sample_meta_csv):
        """Test matching read with UMI between fp and target."""
        demux = FASTQDemultiplexer(
            sample_meta_csv,
            mismatch=0,
            min_dist=5,
            max_dist=100
        )
        # fp + UMI + target + rp
        seq = "AGGTCA" + "AAACCCGGG" + "CTAGCT"
        assert demux._match_read(seq) == "sample1"

    def test_match_read_reverse_complement(self, sample_meta_csv):
        """Test reverse-complement matching: rp_rc ... fp_rc."""
        demux = FASTQDemultiplexer(
            sample_meta_csv,
            mismatch=0,
            min_dist=5,
            max_dist=100
        )
        # rp_rc = AGCTAG, fp_rc = TGACCT
        rev_seq = "AGCTAGNNNAGCTGCANNNTGACCT"
        assert demux._match_read(rev_seq) == "sample1"

    def test_match_read_too_close_fails(self, sample_meta_csv):
        """test fails if fp and rp are too close (below min_dist)."""
        demux = FASTQDemultiplexer(
            sample_meta_csv,
            mismatch=0,
            min_dist=10,
            max_dist=100
        )
        # dist between fp and rp is 9 < min_dist
        seq = "AGGTCA" + "AAACCCGGG" + "CTAGCT"
        assert demux._match_read(seq) == "unknown"

    def test_match_read_too_far_fails(self, sample_meta_csv):
        """test fails if fp and rp are too far (above min_dist)."""
        demux = FASTQDemultiplexer(
            sample_meta_csv,
            mismatch=0,
            min_dist=5,
            max_dist=50
        )
        # dist between fp and rp is 100 > max_dist
        long_spacer = "N" * 100
        seq = "AGGTCA" + long_spacer + "GTAGCT"
        assert demux._match_read(seq) == "unknown"

    def test_demultiplex_reads_uncompressed_fastq(
        self, 
        sample_meta_csv, 
        temp_fastq,
        tmp_path,
    ):
        """Test demultiplexing uncompressed FASTQ."""
        output_dir = tmp_path / "demultiplexed"
        demux = FASTQDemultiplexer(
            sample_meta_csv,
            mismatch=0,
            output_dir=output_dir
        )
        demux.demultiplex(temp_fastq)
        demux.close()
        # Check ouput file was created
        output_file = output_dir / "unknown.fq.gz"
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        
    def test_demultiplex_reads_gzipped_fastq(
        self, 
        sample_meta_csv, 
        temp_fastq_gz,
        tmp_path,
    ):
        """Test demultiplexing gzipped FASTQ."""
        output_dir = tmp_path / "demultiplexed"
        demux = FASTQDemultiplexer(
            sample_meta_csv,
            mismatch=0,
            output_dir=output_dir,
        )
        demux.demultiplex(temp_fastq_gz)
        demux.close()
        output_file = output_dir / "unknown.fq.gz"
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        
    def test_counts_dict_returns_copy(self, sample_meta_csv):
        """Test counts_dict returns a copy, not reference."""
        demux = FASTQDemultiplexer(sample_meta_csv)
        demux.counts["sample1"] = 42
        demux.counts["unknown"] = 8
        original = demux.counts_dict
        assert isinstance(original, dict)

        # Modify returned dict should not affect internal state
        original["sample1"] = 99999

        # Internal counts unchanged
        assert demux.counts_dict["sample1"] != 99999
        assert demux.counts["sample1"] == 42

    def test_property_counts_dict_type_hint(self, sample_meta_csv):
        """Test counts_dict has correct type annotations."""
        demux = FASTQDemultiplexer(sample_meta_csv)
        counts = demux.counts_dict
        assert isinstance(counts, dict)
        for k, v in counts.items():
            assert isinstance(k, str)
            assert isinstance(v, int)

    def test_get_writer_creates_dir(self, sample_meta_csv):
        """Test _get_writer creates output directory."""
        with patch("pathlib.Path.mkdir") as mock_mkdir, \
                patch("gzip.open", return_value=Mock()) as mock_gzip:
            
            demux = FASTQDemultiplexer(sample_meta_csv, output_dir="mock_dir")
            writer = demux._get_writer("sample1")

            mock_mkdir.assert_called_once_with(parents=True, exist_ok=True)
            mock_gzip.assert_called_once()
            assert writer is not None

