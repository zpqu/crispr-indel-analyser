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

"""Unit tests for ParasailAligner."""

from crispr_indel_analyser.analysis.aligner import ParasailAligner


class TestParasailAligner:
    """Test suite for ParasailAligner class."""

    def setup_method(self):
        """Initialize aligner with default parameters before each test."""
        self.aligner = ParasailAligner(
            match=2, 
            mismatch=-1, 
            gap_open=10, 
            gap_extend=1,
        )

    def test_align_global_identical_sequences(self):
        """Test alignment of identical sequences produces no gaps."""
        query = "AAACCCGGG"
        ref = "AAACCCGGG"
        aligned_query, aligned_ref = self.aligner.align_global(query, ref) 
        assert aligned_query == query
        assert aligned_ref == ref
        assert "-" not in aligned_query
        assert "-" not in aligned_ref

    def test_align_with_insertion_in_read(self):
        """Test alignment detects insertion (gap in reference)."""
        query = "AAATTTCCCCTGGG"  # extra 'TTTCCCC' inserted
        ref = "AAATGGG"
        aligned_query, aligned_ref = self.aligner.align_global(query, ref)
        assert len(aligned_query) == len(aligned_ref)
        assert "-"in aligned_ref  # gap in ref
        assert "TTTCCCC" in aligned_query.replace("-", "")

    def test_align_with_deletion_in_read(self):
        """Test alignment detects deletion (gap in query)."""
        query = "AAAGGG"  # missing 'TTTCCCC'
        ref = "AAATTTCCCCGGG"
        aligned_query, aligned_ref = self.aligner.align_global(query, ref)
        assert len(aligned_query) == len(aligned_ref)
        assert "-"in aligned_query  # gap in query
        assert "TTTCCCC" in aligned_ref.replace("-", "")

    def test_align_mismatch_only_no_gap(self):
        """Test alignment with mismatches but no indels."""
        query = "AAGCCCTGG"
        ref = "AAACCCGGG"
        aligned_query, aligned_ref = self.aligner.align_global(query, ref)
        assert "-" not in aligned_query
        assert "-" not in aligned_ref
        assert len(aligned_query) == len(query)
        mismatches = sum(1 for a, b in zip(aligned_query, aligned_ref) if a != b)
        assert mismatches > 0

    def test_align_short_query(self):
        """Test alignment when query is much shorter than reference."""
        query = "CCCGGG"
        ref = "AAATTTCCCGGGTTTAAA"
        aligned_query, aligned_ref = self.aligner.align_global(query, ref)
        assert len(aligned_query) == len(aligned_ref)
        assert len(aligned_ref) == len(ref)
        assert query in aligned_query.replace("-", "")

    def test_align_empty_sequence(self):
        """Test behavior with empty input."""
        query = ""
        ref = "AAACCC"
        aligned_query, aligned_ref = self.aligner.align_global(query, ref)
        assert aligned_query == "-" * len(ref)
        assert aligned_ref == ref
        
    def test_align_reverse_complement_pattern(self):
        """Test alignment can handle reverse-complement like patterns."""
        query = "CCCGGGTTT"
        ref = "AAACCCGGG"
        aligned_query, aligned_ref = self.aligner.align_global(query, ref)
        assert isinstance(aligned_query, str)
        assert isinstance(aligned_ref, str)
        assert len(aligned_query) == len(aligned_ref)

    def test_custom_scoring_parameters(self):
        """Test that custom match/mismatch/gap scores are used."""
        # Strict aligner: penalize gaps more
        strict_aligner = ParasailAligner(
            match=3, 
            mismatch=-3, 
            gap_open=15, 
            gap_extend=2, 
        )
        query = "AAAATCCGGG"  # one mismatch, one insertion
        ref = "AAACCCGGG"
        aligned_query, aligned_ref = strict_aligner.align_global(query, ref)
        assert len(aligned_query) == len(aligned_ref)
        matches = sum(1 for a, b in zip(aligned_query, aligned_ref) if a == b)
        assert matches == 8

