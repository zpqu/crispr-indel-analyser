# src/preprocess/meta_preprocessor.py
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

"""Preprocesses sample metadata CSV with target, flanking, and primer sequences.

This module provides a robust function to load and process a metadata CSV
containing sample-specific sequences (target, up/down flanks, primers).
It cleans invalid characters and generates reverse complements.

Supports both file paths (str, Path) and file-like objects (e.g. StringIO)
for maximum flexibility in production and testing environments.
"""

from typing import Union, TextIO
from os import PathLike
import os
import pandas as pd
from src.utils.helpers import reverse_complement, clean_sequence


def load_and_process_meta_csv(meta_csv: Union[str, PathLike, TextIO]) -> dict:
    """Loads and processes sample metadata CSV into a structured dictionary.

    Accepts:
        - File path as string or Path-like object
        - File-like object (e.g., StringIO, io.TextIOWrapper)

    Cleans all DNA sequences (removes non-ATCG characters) and generates
    reverse complements for each field.

    Args:
        meta_csv: Input source. Can be:
            - A string path to a CSV file
            - A Path-like object (e.g., pathlib.Path)
            - A readable file-like object (must have .read() method)

    Returns:
        dict: Nested dictionary where:
            - Outer key: sample name (from 'sample' column)
            - Inner dict: sequence fields including original and _rc:
              {
                'target': 'AAACCCGGG',
                'target_rc': 'CCCGGGTTT',
                'up': 'TT',
                'up_rc': 'AA',
                ...
              }

    Raises:
        FileNotFoundError: If a path is provided but the file does not exist.
        ValueError: If required columns are missing or parsing fails.
        TypeError: If input is not a valid path or readable file-like object.

    Example:
        # From string path
        data = load_and_process_meta_csv("data/meta.csv")

        # From pathlib.Path
        from pathlib import Path
        data = load_and_process_meta_csv(Path("data/meta.csv"))

        # From StringIO (for testing)
        from io import StringIO
        csv = "sample,target,up,down,fp,rp\\nsample1,ATCG,TT,AA,AGGT,CTAG"
        data = load_and_process_meta_csv(StringIO(csv))
    """
    # Determine if input is a path type (str or PathLike)
    is_path = isinstance(meta_csv, (str, PathLike))

    # Only check file existence if it's a path
    if is_path:
        path_str = str(meta_csv)
        if not os.path.exists(path_str):
            raise FileNotFoundError(f"Metadata file not found: {path_str}")

    # Let pandas handle the actual reading (supports str, Path, file-like)
    try:
        df = pd.read_csv(meta_csv)
    except Exception as e:
        raise ValueError(f"Failed to read CSV input: {e}")

    # Validate required columns
    required = {'sample', 'target', 'up', 'down', 'fp', 'rp'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Sample meta CSV must have columns: {required}. Missing: {missing}")

    # Clean sequences and generate reverse complements
    for col in ['target', 'up', 'down', 'fp', 'rp']:
        df[col] = df[col].astype(str).apply(clean_sequence)
        df[f"{col}_rc"] = df[col].apply(reverse_complement)

    # Return as nested dictionary indexed by 'sample'
    return df.set_index('sample')[[
        'target', 'target_rc',
        'up', 'up_rc',
        'down', 'down_rc',
        'fp', 'fp_rc',
        'rp', 'rp_rc'
    ]].to_dict('index')
