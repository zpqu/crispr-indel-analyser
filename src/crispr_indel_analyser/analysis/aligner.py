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

"""High-performance sequence aligner using parasail."""

import parasail


class ParasailAligner:
    """Wrapper for parasail global (Needleman-Wunsch) alignment."""

    def __init__(
        self, 
        match: int = 2, 
        mismatch: int = -1,
        gap_open: int = 10, 
        gap_extend: int = 1,
    ):
        """
        Args:
            match: Match score (e.g., 2)
            mismatch: Mismatch penalty (e.g., -1)
            gap_open: Gap opening penalty (positive int, e.g., 10)
            gap_extend: Gap extension penalty (positive int, e.g., 1)
        """
        self.matrix = parasail.matrix_create("ACGT", match, mismatch)
        self.gap_open = gap_open
        self.gap_extend = gap_extend

    def align_global(self, query: str, ref: str) -> tuple:
        """Performs global alignment.

        Args:
            query: Query sequence.
            ref: Reference sequence.

        Returns:
            Tuple of query sequences and aligned reference.
        """
        if not query or not ref:
            return query or "-" * len(ref), ref or "-" * len(query)

        # Needleman-Wunsch (guaranteed trackback)
        try:
            result = parasail.nw_trace_scan_sat(
                query, 
                ref,
                self.gap_open, 
                self.gap_extend, 
                self.matrix,
            )
            if result.traceback:
                return result.traceback.query, result.traceback.ref
        except Exception as e:
            print(f"[DEBUG] nw_trace failed: {e}")

        # Ultimate fallback
        return query, ref 

