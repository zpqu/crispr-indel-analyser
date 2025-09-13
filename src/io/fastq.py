# src/io/fastq.py
"""FASTQ reading utilities."""

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

from pathlib import Path
from Bio import SeqIO
import gzip


def read_fastq(fastq_path: str | Path):
    """Reads FASTQ file (supports .gz)

    Args:
        fastq_path: Path to FASTQ file.

    Yields:
        SeqRecord: Biopython sequence record.

    Raises:
        FileNotFoundError: If the specified file does not exist.
    """
    path = Path(fastq_path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {fastq_path}")
    
    opener = gzip.open if path.suffix == ".gz" else open
    mode = "rt"
    with opener(fastq_path, mode) as f:
        for record in SeqIO.parse(f, "fastq"):
            yield record

