#!/usr/bin/env python3
"""
Wrapper for running DiscoSNP++ with resource monitoring.

Runs the DiscoSNP++ pipeline and records wall‑clock time and maximum RSS.

Usage:
  python run_discosnp.py \
      --disco-snp /path/to/run_discoSnp++.sh \
      --reads /data/reads.fastq.gz \
      --output-prefix /data/output/disco_result \
      --kmer-size 31 \
      --max-diff 10 \
      --log-file disco_metrics.log
"""

import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Wrapper for DiscoSNP++ with /usr/bin/time -v monitoring"
    )
    parser.add_argument("--disco-snp", required=True,
                        help="Path to run_discoSnp++.sh script")
    parser.add_argument("--reads", required=True,
                        help="Input reads (FASTQ, can be gzipped)")
    parser.add_argument("--output-prefix", required=True,
                        help="Output prefix (-p argument for DiscoSNP++)")
    parser.add_argument("--kmer-size", type=int, default=31,
                        help="k-mer size (default: 31)")
    parser.add_argument("--max-diff", type=int, default=10,
                        help="-u parameter (max difference, default: 10)")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads (if supported, otherwise ignored)")
    parser.add_argument("--log-file", default="discosnp_metrics.log",
                        help="Log file for resource usage (default: discosnp_metrics.log)")
    parser.add_argument("--work-dir", default=".",
                        help="Working directory (default: current directory)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the command without running")

    args = parser.parse_args()

    # Change to work directory
    work_dir = Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(work_dir)

    log_file = Path(args.log_file)
    log_file.write_text("")  # clear

    # Build DiscoSNP++ command
    cmd = [
        args.disco_snp,
        "-u", str(args.max_diff),
        "-r", args.reads,
        "-T",                  # (the original had -T, presumably something)
        "-p", args.output_prefix,
        "-k", str(args.kmer_size),
    ]
    # If threads are supported (DiscoSNP++ can use -t), add:
    if args.threads > 1:
        cmd.extend(["-t", str(args.threads)])

    print("Starting DiscoSNP++ with monitoring...")
    print("-" * 60)
    print("Command:", " ".join(cmd))

    if args.dry_run:
        print("Dry run: not executing.")
        return

    time_tmp = "time.tmp"
    full_cmd = ["/usr/bin/time", "-v", "-o", time_tmp] + cmd

    start = time.time()
    proc = subprocess.run(full_cmd, capture_output=True, text=True)
    elapsed = time.time() - start

    if proc.returncode != 0:
        print("Error: DiscoSNP++ failed.")
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
        log.write("[DiscoSNP++]\n")
        log.write(f"  Command: {' '.join(cmd)}\n")
        log.write(f"  Wall Time: {wall_time}\n")
        log.write(f"  Max RSS: {max_rss}\n")
        log.write("-" * 40 + "\n")

    print("-" * 60)
    print("DiscoSNP++ completed.")
    print("Summary of Metrics:")
    with open(log_file) as f:
        print(f.read())


if __name__ == "__main__":
    main()

# python run_discosnp.py \
#     --disco-snp /tools/DiscoSNP/run_discoSnp++.sh \
#     --reads /data/SRR30860786_1.fastq.gz \
#     --output-prefix /data/output/discotau \
#     --kmer-size 31 \
#     --max-diff 10 \
#     --threads 8