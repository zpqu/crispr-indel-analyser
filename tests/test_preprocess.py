# tests/test_preprocess.py
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
