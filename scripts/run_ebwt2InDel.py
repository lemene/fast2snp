#!/usr/bin/env python3
"""
Wrapper for the ebwt2InDel variant calling pipeline (BCR_LCP_GSA → ebwt2InDel → filter_snp).

This is a PRODUCTION wrapper, not a test harness. It runs the three steps
sequentially and records wall‑clock time and maximum memory for each step.

Usage:
  python run_ebwt2indel_pipeline.py \
      --fastq /data/input.fastq.gz \
      --bcr-lcp-gsa /tools/BCR_LCP_GSA/BCR_LCP_GSA \
      --ebwt2indel /tools/ebwt2InDel/build/ebwt2InDel \
      --filter-snp /tools/ebwt2InDel/build/filter_snp \
      --output-base /data/output/ebwt_result \
      --threads 32
"""

import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path


def run_and_monitor(step_name: str, cmd: list, log_file: Path):
    """
    Execute a command with /usr/bin/time -v.
    Prints wall time and max RSS, appends them to log_file.
    Exits on failure.
    """
    print("-" * 60)
    print(f"Step: {step_name}")
    print("Command:", " ".join(cmd))

    time_tmp = "time.tmp"
    full_cmd = ["/usr/bin/time", "-v", "-o", time_tmp] + cmd

    start = time.time()
    proc = subprocess.run(full_cmd, capture_output=True, text=True)
    elapsed = time.time() - start

    if proc.returncode != 0:
        print(f"Error: Step '{step_name}' failed.")
        print(proc.stderr)
        if Path(time_tmp).exists():
            os.remove(time_tmp)
        sys.exit(1)

    # Parse time output
    wall_time = "N/A"
    max_rss = "N/A"
    try:
        with open(time_tmp) as f:
            output = f.read()
        m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)", output)
        if m:
            wall_time = m.group(1).strip()
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", output)
        if m:
            max_rss = m.group(1) + " KB"
    except Exception:
        pass
    finally:
        if Path(time_tmp).exists():
            os.remove(time_tmp)

    print(f"  - Wall Time: {wall_time}")
    print(f"  - Max RSS: {max_rss}")

    with open(log_file, "a") as log:
        log.write(f"[{step_name}]\n")
        log.write(f"  Command: {' '.join(cmd)}\n")
        log.write(f"  Wall Time: {wall_time}\n")
        log.write(f"  Max RSS: {max_rss}\n")
        log.write("-" * 40 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Wrapper for BCR_LCP_GSA → ebwt2InDel → filter_snp pipeline"
    )
    parser.add_argument("--fastq", required=True,
                        help="Input FASTQ file (gzipped)")
    parser.add_argument("--bcr-lcp-gsa", required=True,
                        help="Path to BCR_LCP_GSA executable")
    parser.add_argument("--ebwt2indel", required=True,
                        help="Path to ebwt2InDel executable")
    parser.add_argument("--filter-snp", required=True,
                        help="Path to filter_snp executable")
    parser.add_argument("--output-base", required=True,
                        help="Base name for outputs (e.g. /path/to/result_base)")
    parser.add_argument("--threads", type=int, default=36,
                        help="Threads for ebwt2InDel (default: 36)")
    parser.add_argument("--filter-threshold", type=int, default=5,
                        help="Threshold for filter_snp (default: 5)")
    parser.add_argument("--work-dir", default=".",
                        help="Working directory (default: current directory)")
    parser.add_argument("--force-bcr", action="store_true",
                        help="Force re-running BCR_LCP_GSA even if .ebwt exists")
    parser.add_argument("--log-file", default="pipeline_metrics.log",
                        help="Log file name (default: pipeline_metrics.log)")

    args = parser.parse_args()

    # Change to working directory
    work_dir = Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(work_dir)

    log_file = Path(args.log_file)
    log_file.write_text("")  # clear previous log

    output_base = args.output_base
    bwt_output = output_base + ".bwt"          # BCR_LCP_GSA output prefix
    ebwt_file = output_base + ".ebwt"          # generated ebwt file
    snp_output = output_base + ".snp"
    filtered_snp = output_base + "_filtered.snp"

    print("Starting ebwt2InDel pipeline (wrapper)...")
    total_start = time.time()

    # Step 1: BCR_LCP_GSA
    if not args.force_bcr and Path(ebwt_file).exists():
        print(f"\nStep 1 output ({ebwt_file}) already exists. Skipping BCR_LCP_GSA.")
    else:
        cmd1 = [args.bcr_lcp_gsa, args.fastq, bwt_output]
        run_and_monitor("BCR_LCP_GSA", cmd1, log_file)

    # Step 2: ebwt2InDel
    if not Path(ebwt_file).exists():
        sys.exit(f"Error: Required input for ebwt2InDel not found: {ebwt_file}")
    cmd2 = [args.ebwt2indel, "-t", str(args.threads), "-1", ebwt_file, "-o", snp_output]
    run_and_monitor("ebwt2InDel", cmd2, log_file)

    # Step 3: filter_snp (requires redirection, so use bash -c)
    cmd3 = [args.filter_snp, snp_output, str(args.filter_threshold)]
    # To redirect stdout, wrap with bash -c
    shell_cmd = f"{' '.join(cmd3)} > {filtered_snp}"
    run_and_monitor("filter_snp", ["bash", "-c", shell_cmd], log_file)

    total_elapsed = int(time.time() - total_start)
    print("-" * 60)
    print(f"Pipeline completed successfully.")
    print(f"Total execution time: {total_elapsed} seconds.")
    print("Summary of Metrics:")
    with open(log_file) as f:
        print(f.read())


if __name__ == "__main__":
    main()

# python run_ebwt2indel_pipeline.py \
#     --fastq /data/input.fastq.gz \
#     --bcr-lcp-gsa /tools/BCR_LCP_GSA/BCR_LCP_GSA \
#     --ebwt2indel /tools/ebwt2InDel/build/ebwt2InDel \
#     --filter-snp /tools/ebwt2InDel/build/filter_snp \
#     --output-base /data/output/ebwt_result \
#     --threads 32 \
#     --work-dir /tmp