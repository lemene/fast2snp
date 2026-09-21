#!/usr/bin/env python3
"""
run_validation.py — Validate fast2snp calls by mapping to reference and benchmarking.

Steps:
  1. Convert .snp bubbles to VCF via BWA-MEM (using snp_bubble_2_vcf.py)
     Output: <output_base>/validation_mapping/
  2. Evaluate VCF against truth set with hap.py (using run_singularity_happy.py)
     Output: <output_base>/validation_benchmark/

Usage:
  python run_validation.py \
      --snp /path/to/calls.snp \
      --ref /path/to/ref.fa \
      --ref-full /path/to/full_ref.fa \
      --truth /path/to/truth.vcf.gz \
      --bed /path/to/confident.bed \
      --output-base /path/to/results \
      --prefix my_sample
"""

import argparse
import subprocess
import sys
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

SCRIPT_DIR = Path(__file__).resolve().parent
# Assumes snp_bubble_2_vcf.py and run_singularity_happy.py are in the same directory
SNP2VCF_SCRIPT = SCRIPT_DIR / "snp_bubble_2_vcf.py"
HAPPY_SCRIPT = SCRIPT_DIR / "run_singularity_happy.py"

def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate fast2snp variants: map to reference and benchmark."
    )
    parser.add_argument("--snp", required=True,
                        help="Input .snp file from fast2snp")
    parser.add_argument("--ref", required=True,
                        help="Reference FASTA for mapping (e.g., chr22.fa)")
    parser.add_argument("--truth", required=True,
                        help="Truth VCF.gz")
    parser.add_argument("--bed", help="High-confidence BED file")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Base directory for validation outputs")
    parser.add_argument("--prefix", default="variants",
                        help="Sample prefix for output naming")
    parser.add_argument("--threads", type=int, default=8,
                        help="Threads for mapping and hap.py (default: 8)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing")
    parser.add_argument("--happy-extra", default="",
                        help="Extra arguments to pass to run_singularity_happy.py")
    
    return parser.parse_args()

def check_args(args):
    # Check scripts exist
    if not SNP2VCF_SCRIPT.is_file():
        sys.exit(f"Error: {SNP2VCF_SCRIPT} not found.")
    if not HAPPY_SCRIPT.is_file():
        sys.exit(f"Error: {HAPPY_SCRIPT} not found.")

def run_cmd(cmd, dry_run=False):
    """Print and optionally execute a command."""
    print("  Running:", " ".join(cmd))
    if not dry_run:
        subprocess.run(cmd, check=True)

def get_mapping_dir(args):
    """Get the directory for mapping outputs."""
    return Path(args.outdir) / "validation_mapping"

def get_benchmark_dir(args):
    """Get the directory for benchmark outputs."""
    return Path(args.outdir) / "validation_benchmark"

def get_query_vcf(args):
    """Get the query VCF file."""
    mapping_dir = get_mapping_dir(args)
    return mapping_dir / f"{args.prefix}.vcf"

def get_happy_summary_file(args):
    return get_benchmark_dir(args) / f"{args.prefix}.summary.csv"

def run_mapping(args):
  # Create output directories
    out_base = Path(args.outdir)
    mapping_dir = get_mapping_dir(args)
    mapping_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1: map .snp to reference -> VCF
    # ------------------------------------------------------------------
    print("=" * 60)
    print("Step 1: Mapping bubbles to reference genome")
    print("=" * 60)
    if not Path(get_query_vcf(args)).exists() or \
        Path(get_query_vcf(args)).stat().st_atime < Path(args.snp).stat().st_mtime:

        snp2vcf_cmd = [
            sys.executable, str(SNP2VCF_SCRIPT),
            args.snp,
            args.ref,
            str(mapping_dir),
            args.prefix,
            "--threads", str(args.threads),
        ]
        run_cmd(snp2vcf_cmd, dry_run=args.dry_run)
    else:
        print("  Skipping:", get_query_vcf(args))

def run_benchmark(args):
    benchmark_dir = get_benchmark_dir(args)
    mapping_dir = get_mapping_dir(args)

    query_vcf = get_query_vcf(args)

    # ------------------------------------------------------------------
    # Step 2: benchmark with hap.py
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("Step 2: Benchmarking with hap.py")
    print("=" * 60)

    if not Path(get_happy_summary_file(args)).exists() or \
        Path(get_happy_summary_file(args)).stat().st_atime < Path(query_vcf).stat().st_mtime:
        happy_cmd = [
            sys.executable, str(HAPPY_SCRIPT),
            "-t", args.truth,
            "-q", str(query_vcf),
            "-r", args.ref,
            "-o", str(benchmark_dir),
            "--prefix", args.prefix,
            "-p", str(args.threads),
        ]
        if args.bed:
            happy_cmd += ["-b", args.bed]
        print(happy_cmd)
        run_cmd(happy_cmd, dry_run=args.dry_run)
    else:
        print("  Skipping:", get_happy_summary_file(args))

    print("\nValidation finished.")
    print(f"  Mapped VCF: {query_vcf}")
    print(f"  Benchmark results: {benchmark_dir}")

def main():
    args = parse_args()
    check_args(args)

    run_mapping(args)
    run_benchmark(args)


if __name__ == "__main__":
    main()