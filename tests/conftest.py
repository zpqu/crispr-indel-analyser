# tests/conftest.py
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

"""Pytest fixtures for CRISPR pipeline unit tests.

This module defines reusable test fixtures for:
- Sample metadata CSV (simulated file input)
- Expected processed reference and barcode data
These fixtures are used across multiple test modules to ensure consistency
and reduce code duplication.

All sequence data is synthetic and for testing only.
"""
import pytest
from io import StringIO


@pytest.fixture
def sample_meta_csv():
    """Provides a StringIO object simulating a metadata CSV file."""
    csv = """sample,target,up,down,fp,rp
sample1,AAACCCGGG,TT,AA,AGGTCA,CTAGCT
sample2,TGCAGCTA,CA,TG,GATCCA,CGAAGT"""
    return StringIO(csv)


@pytest.fixture
def processed_meta_data():
    """Provides expected output dictionary for metadata preprocessing."""
    return {
        'sample1': {
            'target': 'AAACCCGGG',
            'target_rc': 'CCCGGGTTT',
            'up': 'TT',
            'up_rc': 'AA',
            'down': 'AA',
            'down_rc': 'TT',
            'fp': 'AGGTCA',
            'fp_rc': 'TGACCT',  # reverse_complement("AGGTCA")
            'rp': 'CTAGCT',
            'rp_rc': 'AGCTAG'   # reverse_complement("CTAGCT")
        },
        'sample2': {
            'target': 'TGCAGCTA',
            'target_rc': 'TAGCTGCA',  # reverse_complement("TGCAGCTA")
            'up': 'CA',
            'up_rc': 'TG',           # reverse_complement("CA")
            'down': 'TG',            # corrected from 'TC' in your original
            'down_rc': 'CA',         # reverse_complement("TG")
            'fp': 'GATCCA',
            'fp_rc': 'TGGATC',       # reverse_complement("GATCCA")
            'rp': 'CGAAGT',
            'rp_rc': 'ACTTCG'        # reverse_complement("CGAAGT")
        }
    }
