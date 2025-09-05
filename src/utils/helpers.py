# src/utils/helpers.py
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
"""Utility functions for DNA sequence manipulation."""

from typing import Union


def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of a DNA sequence.

    Args:
        seq (str): Input DNA sequence.

    Returns:
        str: Reverse complement.
    """
    if not isinstance(seq, str):
        return ""
    trans = str.maketrans("ACGT", "TGCA")
    return seq.upper().translate(trans)[::-1]


def clean_sequence(seq: Union[str, None]) -> str:
    """Cleans a DNA sequence by removing invalid characters.

    Args:
        seq (str or None): Input sequence.

    Returns:
        str: Cleaned ATCG-only sequence.
    """
    if (
        not isinstance(seq, str)
        or not seq.strip()
        or seq.lower() in {'nan', 'none', 'null'}
    ):
        return ""
    return ''.join(c.upper() for c in seq if c.upper() in 'ATCG')
