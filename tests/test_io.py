# tests/test_io.py
"""Tests for I/O modules."""

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

import pandas as pd
import pytest
import gzip
from pathlib import Path 
from io import StringIO
from src.io.writer import write_table
from src.io.fastq import read_fastq


def test_write_table(tmp_path):
    """Test write_table creates correct file."""
    df = pd.DataFrame([{'A': 1, 'B': 2}])
    filepath = tmp_path / "test.txt"
    write_table(df, filepath)
    assert filepath.exists()

    with filepath.open() as f:
        content = f.read()
    assert "A\tB" in content
    assert "1\t2" in content


def test_read_fastq(tmp_path):
    """Test read_fastq reads valid FASTQ."""
    fastq_content = "@read1\nATCG\n+\nIIII\n"
    fastq_path = tmp_path / "test.fastq"
    fastq_path.write_text(fastq_content)

    records = list(read_fastq(fastq_path))
    assert len(records) == 1
    assert str(records[0].seq) == "ATCG"
    assert records[0].id == "read1"
    


def test_read_fastq_gz(tmp_path):
    """Test read_fastq reads gzipped FASTQ."""
    fastq_content = "@read1\nATCG\n+\nIIII\n"
    fastq_path = tmp_path / "test.fastq.gz"

    with gzip.open(fastq_path, 'wt') as f:
        f.write(fastq_content)

    records = list(read_fastq(fastq_path))
    assert len(records) == 1
    assert str(records[0].seq) == "ATCG"
    assert records[0].id == "read1"

