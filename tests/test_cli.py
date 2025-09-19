# Copyright (C) 2025 Zhipeng Qu
#
# Permission is hereby granted, free of charge, to any person obtaining a copy # of this software and associated documentation files (the "Software"), to deal
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

"""Integration tests for main.py CLI"""

import gzip
import subprocess
import pytest
from pathlib import Path


@pytest.fixture
def temp_dirs(tmp_path):
    """Creates temporary input/output directories with test data."""
    # Input data directory
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    # Output directories
    results_dir = tmp_path / "results"
    demux_dir = tmp_path / "demultiplexed"
    results_dir.mkdir()
    demux_dir.mkdir()

    # Create minimal FASTQ (one forward read)
    fastq_gz = data_dir / "test.fq.gz"
    with gzip.open(fastq_gz, "wt") as f:
        f.write(
            "@read1\nAGGTCATTCAAACCCGGGAAGAGCTAG"
            "\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        )

    # Create metadata CSV
    meta_csv = data_dir / "meta.csv"
    meta_csv.write_text("""sample,target,up,down,fp,rp
sample1,AAACCCGGG,TTC,AAG,AGGTCA,CTAGCT""")

    return {
        "data_dir": data_dir,
        "results_dir": results_dir,
        "demux_dir": demux_dir,
        "fastq": fastq_gz,
        "meta_csv": meta_csv,
    }


def test_main_cli_run_successful(temp_dirs):
    """Test full pipeline runs successfully with valid inputs."""
    cmd = [
        "python", "main.py",
        "--fastq", str(temp_dirs["fastq"]),
        "--meta-csv", str(temp_dirs["meta_csv"]),
        "--demux-dir", str(temp_dirs["demux_dir"]),
        "--result-dir", str(temp_dirs["results_dir"]),
        "--flank-mismatch", "0",
        "--min-dist", "1",
        "--max-dist", "100",
        "--output-format", "txt"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"Command failed: {result.stderr}"
    assert "Demultiplexing completed" in result.stderr
    assert "Analysing indels..." in result.stderr
    assert "Combined results saved" in result.stderr

    # Check output files exist and are non-empty
    demux_file = temp_dirs["demux_dir"] / "sample1.fq.gz"
    assert demux_file.exists()
    assert demux_file.stat().st_size > 0

    summary_txt = temp_dirs["results_dir"] / "sample1_summary.txt"
    assert summary_txt.exists()
    content = summary_txt.read_text()
    assert "num_other" in content
    assert "num_ins" in content

    combined_csv = temp_dirs["results_dir"] / "summary.csv"
    assert combined_csv.exists()
    csv_content = combined_csv.read_text()
    assert "sample1" in csv_content
    assert "per_del" in csv_content


def test_main_cli_json_output(temp_dirs):
    """Test JSON output includes structured position data."""
    json_results = temp_dirs["results_dir"].parent / "json_results"
    json_results.mkdir(exist_ok=True)

    cmd = [
        "python", "main.py",
        "--fastq", str(temp_dirs["fastq"]),
        "--meta-csv", str(temp_dirs["meta_csv"]),
        "--result-dir", str(json_results),
        "--output-format", "json"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0

    json_file = json_results / "sample1_summary.json"
    assert json_file.exists()
    content = json_file.read_text()
    assert '"pos_ins"' in content or '"pos_del"' in content
    assert '"num_other"' in content


def test_main_reverse_strand_read(temp_dirs):
    """Test analyser handles reverse-complement reads."""
    rev_fastq = temp_dirs["data_dir"] / "rev_test.fq.gz"
    # Reverse complement of TTAAACCCGGGAA → TTCCCGGGTTTAA
    with gzip.open(rev_fastq, "wt") as f:
        f.write(
            "@read2\nCTAGCTCTTCCCGGGTTTGAATGACCT"
            "\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        )

    cmd = [
        "python", "main.py",
        "--fastq", str(rev_fastq),
        "--meta-csv", str(temp_dirs["meta_csv"]),
        "--result-dir", str(temp_dirs["results_dir"]),
        "--demux-dir", str(temp_dirs["demux_dir"]),
        "--flank-mismatch", "1"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    assert "sample1" in result.stderr

    final_csv = temp_dirs["results_dir"] / "summary.csv"
    assert final_csv.exists()
    assert "num_other" in final_csv.read_text()


def test_main_invalid_fastq_file(temp_dirs):
    """Test graceful handling of missing FASTQ."""
    result = subprocess.run([
        "python", "main.py",
        "--fastq", "nonexistent.fq.gz",
        "--meta-csv", str(temp_dirs["meta_csv"]),
        "--result-dir", str(temp_dirs["results_dir"])
    ], capture_output=True, text=True)

    assert result.returncode != 0
    assert "FASTQ file not found" in result.stderr

def test_main_invalid_meta_csv(tmp_path):
    """Test invalid CSV structure is caught."""
    bad_csv = tmp_path / "bad.csv"
    bad_csv.write_text("sample,missing_col\ns1,val")

    result = subprocess.run([
        "python", "main.py",
        "--fastq", "data/test.fq.gz",
        "--meta-csv", str(bad_csv),
        "--result-dir", str(tmp_path / "out")
    ], capture_output=True, text=True)

    assert result.returncode != 0
    assert "Missing: " in result.stderr or "ValueError" in result.stderr


def test_main_help_shows_usage():
    """Test --help prints usage information."""
    result = subprocess.run(
        ["python", "main.py", "--help"], 
        capture_output=True, text=True
    )

    assert result.returncode == 0
    assert "--fastq" in result.stdout
    assert "--flank-mismatch" in result.stdout
    assert "--output-format" in result.stdout


def test_main_no_samples_matched(temp_dirs):
    """Test behavior when no reads match any sample."""
    # FASTQ with no matching primers
    bad_fastq = temp_dirs["data_dir"] / "no_match.fq.gz"
    with gzip.open(bad_fastq, "wt") as f:
        f.write("@readX\nATGCATGCATGC\n+\nIIIIIIIIIIII\n")

    cmd = [
        "python", "main.py",
        "--fastq", str(bad_fastq),
        "--meta-csv", str(temp_dirs["meta_csv"]),
        "--result-dir", str(temp_dirs["results_dir"])
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    assert "Sample counts:" in result.stderr
    assert "'unknown'" in result.stderr

