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

"""Utility functions for DNA sequence manipulation."""

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


def clean_sequence(seq: str | None) -> str:
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


def hamming_distance(s1: str, s2: str) -> int | float:
    """Calculates Hamming distance between two strings.

    Args:
        s1: First string.
        s2: Second string.

    Returns:
        Integer distance or float('inf') if length mismatch.

    Raises:
        TypeError if s1 or s2 is not string type.
    """
    if not isinstance(s1, str) or not isinstance(s2, str):
        raise TypeError("Both arguments must be strings")

    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

