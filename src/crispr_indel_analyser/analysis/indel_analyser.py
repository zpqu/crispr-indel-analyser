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

"""Indel frequencey and position analyser."""

from pathlib import Path
from typing import Optional, Literal, Any
import json
import pandas as pd
from crispr_indel_analyser.utils.helpers import hamming_distance, reverse_complement
from crispr_indel_analyser.io.fastq import read_fastq
from crispr_indel_analyser.analysis.aligner import ParasailAligner
from crispr_indel_analyser.io.writer import write_table


class IndelAnalyser:
    """Analyses indel frequency and positions in CRISPR-edited reads."""

    def __init__(
        self, 
        meta_data: dict[str, dict], 
        demux_dir: Path | str, 
        result_dir: Path | str, 
        mismatch: int = 1, 
        ) -> None:
        self.meta_data = meta_data
        self.demux_dir = Path(demux_dir)
        self.result_dir = Path(result_dir)
        self.mismatch = mismatch  # allowed mismatches in flanking sequences
        self.aligner = ParasailAligner(match=2, 
            mismatch=-1, 
            gap_open=10, 
            gap_extend=1,
        )
        self.result_summary: dict[str, dict] = {}
        
    def _match_flanks(self, seq: str, up: str, down: str) -> tuple[bool, int, int]:
        """Checks if sequence contains correct flanking sequences."""
        if not up or not down:
            return False, -1, -1

        for i in range(len(seq) - len(up) + 1):
            if hamming_distance(seq[i:i+len(up)], up) <= self.mismatch:
                up_end = i + len(up)
                for j in range(up_end, len(seq) - len(down) + 1):
                    if hamming_distance(seq[j:j+len(down)], down) <= self.mismatch:
                        return True, up_end, j
        return False, -1, -1

    def analyse(
        self, 
        sample_id:str, 
        output_format: Literal["txt", "json"] = "txt", 
    ) -> Optional[dict[str, Any]]:
        """ Analyses one sample and saves result in designated format.

        Args:
            sample_id: Sample name.
            output_format: Output format, 'txt' (tabular) or 'json' (structured).

        Returns:
            Analysis results in dictionary format.
        """
        if sample_id not in self.meta_data:
            print(f"[WARNING] Sample {sample_id} not in metadata.")
            return None

        fastq_path = self.demux_dir / f"{sample_id}.fq.gz"
        if not fastq_path.exists():
            print(f"[WARNING] FASTQ not found: {fastq_path}")
            return None

        info = self.meta_data[sample_id]
        up, down = info["up"], info["down"]
        target = info["target"]

        stats = {
            "num_ins": 0, 
            "num_del": 0,
            "num_other": 0, 
            "num_skip": 0,
        }

        del_pos: dict[int, int] = {}
        ins_pos: dict[tuple[int, str], int] = {}

        for record in read_fastq(fastq_path):
            seq = str(record.seq)
            rc_seq = reverse_complement(seq)
            matched, up_end, down_start = self._match_flanks(
                seq, 
                up, 
                down, 
            )
            rc_matched, rc_up_end, rc_down_start = self._match_flanks(
                rc_seq, 
                up, 
                down, 
            )
            if not matched and not rc_matched:
                stats["num_skip"] += 1
                continue
            if matched:
                window_seq = seq[up_end:down_start]
            else:
                window_seq = rc_seq[rc_up_end:rc_down_start]

            aligned_query, aligned_ref = self.aligner.align_global(window_seq, target)

            if "-" in aligned_query and "-" not in aligned_ref:
                ref_pos = aligned_query.index("-")
                del_pos[ref_pos] = del_pos.get(ref_pos, 0) + 1
                stats["num_del"] += 1
            elif "-" in aligned_ref and "-" not in aligned_query:
                start = aligned_ref.index("-")
                inserted = ""
                i = start
                while i < len(aligned_ref) and aligned_ref[i] == "-":
                    inserted += aligned_query[i]
                    i += 1
                pos_key = f"after:_{start}:{inserted.upper()}"
                ins_pos[pos_key] = ins_pos.get(pos_key, 0) + 1
                stats["num_ins"] += 1
            else:
                stats["num_other"] += 1

            # Calculate percentage of indels
        total = (
            stats["num_ins"] + 
            stats["num_del"] + 
            stats["num_other"]
        )
        stats["per_ins"] = (
            round(stats["num_ins"] / total * 100, 2) 
            if total 
            else 0
        )
        stats["per_del"] = (
            round(stats["num_del"] / total * 100, 2)
            if total
            else 0
        )

        # Add position data
        stats["pos_ins"] = dict(ins_pos)
        stats["pos_del"] = dict(del_pos)

        # Save in requested format
        output_file = self.result_dir / f"{sample_id}_summary.{output_format}"
        self.result_dir.mkdir(parents=True, exist_ok=True)

        try:
            if output_format == "json":
                with open(output_file, "w") as f:
                    json.dump({"sample": sample_id, **stats}, f, indent=2)
            else:
                df = pd.DataFrame([{"sample": sample_id, **stats}])
                write_table(df, output_file)
        except Exception as e:
            print(f"[ERROR] Failed to write {output_file}: {e}")
            return None
        
        self.result_summary[sample_id] = stats
        return stats

    # Summarise results from all samples
    def summarise_all(
        self, 
        output_csv: Path | str = "results/summary.csv", 
    ) -> None:
        """Saves combined summmary table."""
        if not self.result_summary:
            return
        df = pd.DataFrame.from_dict(self.result_summary, orient="index")
        df.index.name = "sample"
        write_table(df.reset_index(), output_csv, sep=",")

