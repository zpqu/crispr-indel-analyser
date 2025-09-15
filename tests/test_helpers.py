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

"""Tests for helpers module."""

import pytest
import math
from src.utils.helpers import reverse_complement, clean_sequence, hamming_distance


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


# unit test for hamming_distance
@pytest.mark.parametrize("s1, s2, expected", [
    # same strings
    ("", "", 0),
    ("a", "a", 0),
    ("abc", "abc", 0),

    # different string with same length
    ("a", "b", 1),
    ("abc", "abd", 1),
    ("kittens", "sitting", 3),
    ("1011101", "1001001", 2),

    # spectial chars
    ("a b", "a-b", 1),
    ("!@#", "!@$", 1),

    # Unicode
    ("café", "cafe", 1),
    ("naïve", "naive", 1), 
])
def test_hamming_distance_same_length(s1, s2, expected):
    """Test hamming distance for strings of the same length."""
    assert hamming_distance(s1, s2) == expected


@pytest.mark.parametrize("s1, s2", [
    ("a", "ab"),
    ("", "a"),
    ("long", "short"),
    ("abc", "abcd"),
    ("longer", "short"),
])
def test_hamming_distance_different_lengths(s1, s2):
    """Test that strings of different lengths return inf."""
    result = hamming_distance(s1, s2)
    assert result == float('inf')
    assert math.isinf(result)
    assert result > 0


def test_hamming_distance_infinity_properties():
    """Test specific properties of the returned infinity value."""
    result = hamming_distance("a", "ab")

    assert math.isinf(result)
    assert result == float('inf')
    assert result > 0
    assert result > 1000000

