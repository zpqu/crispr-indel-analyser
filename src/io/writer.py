# src/io/writer.py
"""File writing utilities."""

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
import gzip
from pathlib import Path


def open_gzip_writer(filepath: str | Path) -> gzip.GzipFile:
    """ Opens a gzipped writer.

    Args:
        filepath: Output path.

    Returns:
        Gzip file handle.
    """
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)
    return gzip.open(path, 'wt')


def write_table(df: pd.DataFrame, filepath: str | Path, sep: str = '\t') -> None:
    """Writes DataFrame to file.

    Args:
        df: Data to write.
        filepath: Output path.
        sep: Column separator.
    """
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(path, sep=sep, index=False)

