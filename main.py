#!/usr/bin/env python3
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

"""Main entry point for CRISPR Indel Analyser pipeline."""

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path

from src.preprocess.meta_preprocessor import load_and_process_meta_csv
from src.demux.fastq_demux import FASTQDemultiplexer
from src.analysis.indel_analyser import IndelAnalyser


# Configure logger
logging.basicConfig(
    level=logging.INFO, 
    format="[%(levelname)s] %(asctime)s %(message)s", 
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)


def main():
    """Run the full CRISPR analysis pipeline."""
    parser = argparse.ArgumentParser(
        description=(
            "CRISPR Indel Analyser: "
            "Demultiplex and analyse CRISPR editing efficiency"
        ),
    )

    # Input/Output Options
    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument(
        "--fastq",
        type=Path,
        required=True,
        help=(
            "Input FASTQ file path."
            "Supports plain FASTQ (.fq, .fastq) and gzipped (.fq.gz, .fastq.gz)."
        )
    )
    io_group.add_argument(
        "--meta-csv",
        type=Path,
        required=True,
        help=(
            "Path to metadata CSV with columns: "
            "sample, target, up, down, fp, rp. "
            "Defines reference sequences and barcodes for demux and analysis."
        )
    )
    io_group.add_argument(
        "--demux-dir",
        type=Path,
        default="demultiplexed",
        help=(
            "Directory to save demultiplexed FASTQ files. "
            "One .fq.gz per sample will be generated."
        )
    )
    io_group.add_argument(
        "--result-dir",
        type=Path,
        default="results",
        help=(
            "Directory to save analysis results. "
            "Includes per-sample summaries and combined summary csv."
        )
    )

    # Matching Parameters
    match_group = parser.add_argument_group("Matching Parameters")
    match_group.add_argument(
        "--barcode-mismatch", 
        type=int, 
        default=0, 
        choices=[0, 1, 2],
        help=(
            "Allowed mismatches when matching barcodes. "
            "Use 1 if sequencing quality is moderate, 2 only for low-quality data. "
        )
    )
    match_group.add_argument(
        "--min-dist", 
        type=int,
        default=1,
        help=(
            "Minimum distance between forward primer (fp) and reverse primer (rp). "
            "Prevents false positives from random overlaps. "
            "Should be longer than expected UMI/linker length. "
        )
    )
    match_group.add_argument(
        "--max-dist", 
        type=int, 
        default=1000, 
        help=(
            "Maximum allowed distance between fp and rp. "
            "Avoids misassignment due to very long inserts. "
            "Adjust based on amplicon size. "
        )
    )
    match_group.add_argument(
        "--flank-mismatch", 
        type=int, 
        default=0, 
        choices=[0, 1, 2],
        help=(
            "Allowed mismatches when matching upstream/downstream "
            "flanking regions (up/down). Use 1 if sequencing quality "
            "is moderate, 2 only for low-quality data. "
        )
    )

    # Output Format
    output_group = parser.add_argument_group("Output Format")
    output_group.add_argument(
        "--output-format", 
        choices=["txt", "json"],
        default="txt",
        help=(
            "Output format for per-sample summary: "
            "'txt' for tab-separated table (default), "
            "'json' for structured format with indel position statistics. "
        )
    )

    # Verbosity control
    parser.add_argument(
        "--verbose", 
        action="store_true",
        help="Enable debug-level logging."
    )

    args = parser.parse_args()

    # Set log level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)

    logger.info("Program starts")

    # Load metadata
    logger.info("Loading metadata...")
    try:
        meta_data = load_and_process_meta_csv(args.meta_csv)
        logger.debug(f"Loaded metadata for samples: {list(meta_data.keys())}")
    except Exception as e:
        logger.error(f"Failed to read meta CSV: {e}")
        sys.exit(1)

    # Step 1: Demultiplexing
    logger.info("Demultiplexing FASTQ...")
    try:
        demux = FASTQDemultiplexer(
            meta_csv=args.meta_csv, 
            mismatch=args.barcode_mismatch, 
            output_dir=args.demux_dir,
            min_dist=args.min_dist, 
            max_dist=args.max_dist, 
        )
        demux.demultiplex(args.fastq)
        demux.close()
        logger.info("Demultiplexing completed")
        logger.info(f"Sample counts: {dict(demux.counts_dict)}")
    except Exception as e:
        logger.error(f"Demultiplexing failed: {e}")
        sys.exit(1)

    # Step 2: Indel Analysis
    logger.info("Analysing indels...")
    analyser = IndelAnalyser(
        meta_data=meta_data, 
        demux_dir=args.demux_dir, 
        result_dir=args.result_dir, 
        mismatch=args.flank_mismatch, 
    )

    samples = [k for k in demux.counts_dict.keys() if k != "unknown"]
    for sample in samples:
        logger.debug(f"Analysing sample: {sample}")
        result = analyser.analyse(sample, output_format=args.output_format)
        if result is not None:
            logger.info(
                f"{sample}: {result['num_other']} unedited, "
                f"{result['num_ins']} ins ({result['per_ins']:.1f}%)"
                f"{result['num_del']} del ({result['per_del']:.1f}%)"
            )

    # Save combined summary
    try:
        summary_path = args.result_dir / "summary.csv"
        analyser.summarise_all(output_csv=summary_path)
        logger.info(f"Combined results saved to {summary_path}")
    except Exception as e:
        logger.error(f"Failed to save summary: {e}")

    logger.info("Program ends")

if __name__ == "__main__":
    main()

