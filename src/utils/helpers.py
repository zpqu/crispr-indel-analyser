# src/utils/helpers.py
"""Utility functions for DNA sequence manipulation."""

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

