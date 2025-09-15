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

"""Demultiplexes mixed FASTQ by barcode matching.

Uses sliding window with Hamming distance for barcode detection.
"""

import gzip
from Bio import SeqIO
from pathlib import Path
from typing import TextIO
from src.preprocess.meta_preprocessor import load_and_process_meta_csv
from src.utils.helpers import reverse_complement, hamming_distance


class FASTQDemultiplexer:
    """Demultiplexer for mixed CRISPR sequence FASTQ files."""

    def __init__(
            self, 
            meta_csv: str | Path | TextIO, 
            mismatch: int = 0,
            min_dist: int = 20,
            max_dist: int = 200,
            output_dir: str | Path = "demultiplexed"
        ) -> None:
        """Initializes demultiplexer with flexible barcode matching.

        Args:
            meta_csv: Path to metadata CSV with sample barcodes.
            mismatch: Allowed mismatches in barcode (0, 1, 2).
            min_dist: Minimum distance between fp and rp (to avoid false positives).
            max_dist: Maximum distance between fp and rp.
            output_dir: Output directory for demuxed FASTQ file.
        """
        self.meta_data = load_and_process_meta_csv(meta_csv)
        self.mismatch = mismatch
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.output_dir = Path(output_dir)
        self.writers: dict[str, gzip.GzipFile] = {}
        self.counts: dict[str, int] = {}

    def _get_writer(self, sample: str) -> gzip.GzipFile:
        """Creates a gzip writer for the sample."""
        if sample not in self.writers:
            filepath = self.output_dir / f"{sample}.fq.gz"
            filepath.parent.mkdir(parents=True, exist_ok=True)
            self.writers[sample] = gzip.open(filepath, "wt")
            self.counts[sample] = 0
        return self.writers[sample]

    def _find_subseq_with_mismatch(self, seq: str, pattern: str) -> list[int]:
        """Finds all start positions where pattern matches with <= mismatchs.

        Uses sliding window with Hamming distance.

        Args:
            seq: Query sequence.
            pattern: Subsequence to find.

        Returns:
            List of start indices.
        """
        m = len(seq)
        n = len(pattern)
        if n == 0 or m < n:
            return []
        positions = []

        for i in range(m - n + 1):
            if hamming_distance(seq[i : i + n], pattern) <= self.mismatch:
                positions.append(i)
        return positions

    def _match_read(self, seq: str) -> str:
        """Matches read to sample by locating fp and rp with distance constraints.

        Args:
            seq: Read sequence.

        Returns:
            Sample name or 'unknown'.
        """
        for sample, info in self.meta_data.items():
            fp, rp = info["fp"], info["rp"]
            fp_rc, rp_rc = info["fp_rc"], info["rp_rc"]

            # Forward match: fp before rp
            fp_positions = self._find_subseq_with_mismatch(seq, fp)
            rp_positions = self._find_subseq_with_mismatch(seq, rp)
            for fp_pos in fp_positions:
                for rp_pos in rp_positions:
                    if fp_pos < rp_pos:
                        dist = rp_pos - (fp_pos + len(fp))
                        if self.min_dist <= dist <= self.max_dist:
                            return sample

            # Reverse complement match: rp_rc before fp_rc
            rp_rc_positions = self._find_subseq_with_mismatch(seq, rp_rc)
            fp_rc_positions = self._find_subseq_with_mismatch(seq, fp_rc)
            for rp_rc_pos in rp_rc_positions:
                for fp_rc_pos in fp_rc_positions:
                    if rp_rc_pos < fp_rc_pos:
                        dist = fp_rc_pos - (rp_rc_pos + len(rp_rc))
                        if self.min_dist <= dist <= self.max_dist:
                            return sample
            return "unknown"

    def demultiplex(self, fastq_file: Path | str) -> None:
        """Demultiplexes FASTQ file into samle-specific .fq.gz files.

        Reads are assigned based on presence and spacing of primer sequences.
        Reads without any primers found are put into unknown.fq.gz.

        Args:
            fastq_file: Input FASTQ path (supports .gz).

        Raises:
            FileNotFoundError: If input file does not exist.
        """
        fastq_path = Path(fastq_file)
        if not fastq_path.exists():
            raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")

        opener = gzip.open if fastq_path.suffix.endswith(".gz") else open
        mode = "rt"

        with opener(fastq_path, mode) as f:
            for record in SeqIO.parse(f, "fastq"):
                seq = str(record.seq)
                sample = self._match_read(seq)
                writer = self._get_writer(sample)
                qual = "".join(
                    chr(q + 33) 
                    for q in record.letter_annotations.get("phred_quality", [])
                )
                writer.write(f"@{record.id}\n{seq}\n+\n{qual}\n")
                self.counts[sample] = self.counts.get(sample, 0) + 1

    def close(self) -> None:
        """Closes all output file handles."""
        for f in self.writers.values():
            f.close()

    @property
    def counts_dict(self) -> dict[str, int]:
        """Returns sample-to-read count mapping."""
        return dict(self.counts)

