# tests/test_helpers.py
#
# Copyright (C) 2025 Zhipeng Qu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Tests for helpers module."""

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
