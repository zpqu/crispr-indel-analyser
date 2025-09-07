# tests/test_helpers.py
"""Tests for helpers module."""

# Copyright (C) 2025 Zhipeng Qu
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
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

import pytest
from src.utils.helpers import reverse_complement, clean_sequence


def test_reverse_complement():
    """Test reverse complement."""
    assert reverse_complement("A") == "T"
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("GCAT") == "ATGC"
    assert reverse_complement("") == ""  # test empty string
    assert reverse_complement("atcg") == "CGAT"  # test lowercase


def test_clean_sequence():
    """Test sequence cleaning."""
    assert clean_sequence("AAT!@#CC") == "AATCC"
    assert clean_sequence("NNGGTT") == "GGTT"
    assert clean_sequence("123") == ""
    assert clean_sequence("") == ""
    assert clean_sequence(None) == ""
    assert clean_sequence(123) == ""
    assert clean_sequence("atcg") == "ATCG"
    assert clean_sequence("  ") == ""
    assert clean_sequence("NaN") == ""

