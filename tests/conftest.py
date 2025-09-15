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

"""Pytest fixtures for CRISPR pipeline unit tests.

This module defines reusable test fixtures for:
- Sample metadata CSV (simulated file input)
- Expected processed reference and barcode data
These fixtures are used across multiple test modules to ensure consistency
and reduce code duplication.

All sequence data is synthetic and for testing only.
"""

import gzip
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


@pytest.fixture
def temp_fastq(tmp_path):
    """Creates a temporary FASTQ file for testing."""
    fastq = tmp_path / "test.fq"
    fastq.write_text("@read1\nATCG\n+\nIIII\n")
    return fastq


@pytest.fixture
def temp_fastq_gz(tmp_path):
    """Creates a temporary gzipped FASTQ file."""
    fastq_gz = tmp_path / "test.fq.gz"
    with gzip.open(fastq_gz, "wt") as f:
        f.write("@read1\nATCG\n+\nIIII\n")
    return fastq_gz

