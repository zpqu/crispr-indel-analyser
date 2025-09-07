# tests/test_preprocess.py
#
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
"""Unit tests for metadata preprocessing module.

This module tests the functionality of `load_and_process_meta_csv` from
`src.preprocess.meta_preprocessor`. It verifies that:
- CSV is correctly parsed
- Sequences are cleaned (non-ATCG characters removed)
- Reverse complements are accurately generated
- Output structure matches expected nested dictionary format

Uses pytest fixtures defined in `conftest.py` to ensure consistency
and reduce code duplication.
"""

from src.preprocess.meta_preprocessor import load_and_process_meta_csv
from tests.conftest import sample_meta_csv


def test_load_and_process_meta_csv(sample_meta_csv):
    """Tests that metadata CSV is correctly loaded and processed.

    Validates:
        - All samples are present in output
        - All expected sequence fields (original and _rc) are generated
        - Reverse complements are correct
        - Sequence cleaning is applied

    Uses `sample_meta_csv` fixture to simulate file input without
    relying on disk I/O.

    The test checks both structural integrity and content accuracy.

    Fixture:
        sample_meta_csv: StringIO object containing test CSV data.

    Asserts:
        - Correct number of samples (2)
        - Each sample has 10 sequence fields
        - Original sequences match input
        - Reverse complements are correctly computed

    Example:
        Input row:
            sample1,AAACCCGGG,TT,AA,AGGTCA,CTAGCT

        Expected output:
            'target': 'AAACCCGGG'
            'target_rc': 'CCCGGGTTT'
            'fp_rc': 'TGACCT'  # revcomp of 'AGGTCA'
            'rp_rc': 'AGCTAG'  # revcomp of 'CTAGCT'

    See `processed_meta_data` fixture in conftest.py for full expected output.
    """
    meta_data = load_and_process_meta_csv(sample_meta_csv)

    # Structural checks
    assert 'sample1' in meta_data
    assert 'sample2' in meta_data
    assert len(meta_data) == 2

    # Field count check
    assert len(meta_data['sample1']) == 10  # 5 original + 5 rc
    assert len(meta_data['sample2']) == 10

    # Content validation for sample1
    assert meta_data['sample1']['target'] == 'AAACCCGGG'
    assert meta_data['sample1']['target_rc'] == 'CCCGGGTTT'
    assert meta_data['sample1']['up'] == 'TT'
    assert meta_data['sample1']['up_rc'] == 'AA'
    assert meta_data['sample1']['down'] == 'AA'
    assert meta_data['sample1']['down_rc'] == 'TT'
    assert meta_data['sample1']['fp'] == 'AGGTCA'
    assert meta_data['sample1']['fp_rc'] == 'TGACCT'
    assert meta_data['sample1']['rp'] == 'CTAGCT'
    assert meta_data['sample1']['rp_rc'] == 'AGCTAG'

    # Content validation for sample2
    assert meta_data['sample2']['target'] == 'TGCAGCTA'
    assert meta_data['sample2']['target_rc'] == 'TAGCTGCA'
    assert meta_data['sample2']['up'] == 'CA'
    assert meta_data['sample2']['up_rc'] == 'TG'
    assert meta_data['sample2']['down'] == 'TG'
    assert meta_data['sample2']['down_rc'] == 'CA'
    assert meta_data['sample2']['fp'] == 'GATCCA'
    assert meta_data['sample2']['fp_rc'] == 'TGGATC'
    assert meta_data['sample2']['rp'] == 'CGAAGT'
    assert meta_data['sample2']['rp_rc'] == 'ACTTCG'
