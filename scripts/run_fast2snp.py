#!/usr/bin/env python3
"""
fast2snp pipeline: KMC counting → heterozygous k-mer identification → variant calling
"""

import argparse
import os
import sys
import subprocess
import glob
from pathlib import Path
import shutil
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# Get the directory where this script is located
SCRIPT_DIR = Path(__file__).resolve().parent

def require_executable(executable_name):
    """Check if an executable is available in the system PATH."""
    if shutil.which(executable_name) is None:
        logging.error(f"{executable_name} is not found in PATH. Please install it and try again.")
        sys.exit(1)

def parse_known_args():
    parser = argparse.ArgumentParser(
        description="fast2snp pipeline: reference-free small variant detection"
    )

    # Required parameters
    parser.add_argument("input", nargs="+", help="Input files (space-separated *.fastq or *.fastq.gz)")
    
    # Optional parameters
    parser.add_argument("-k", "--kmer-size", type=int, default=31,
                        help="k-mer size (default: 31)")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Number of threads (default: 8)")
    parser.add_argument("-o", "--outdir", type=str, default='.',
                        help="Output directory (default: current directory)")

    return parser.parse_known_args()


def get_kmer_dir(args):
    return os.path.join(args.outdir, "kmer_analysis")

def get_kmc_input_list_path(args):
    kmer_dir = get_kmer_dir(args)
    return os.path.join(kmer_dir, "kmc_input.list") 

def get_kmc_db_path(args):
    kmer_dir = get_kmer_dir(args)
    return os.path.join(kmer_dir, "kmc_db")

def get_kmc_hist_file(args):
    return os.path.join(get_kmer_dir(args), "kmc_db.histo")

def get_kmc_png_pattern(args):
    return os.path.join(get_kmer_dir(args), "kmc_db_*.png")
def get_snp_dir(args):
    return os.path.join(args.outdir, "variant_calls")

def get_bubble_snp_file(args):
    return os.path.join(get_snp_dir(args), "bubbles.snp")
def run_kmc(args):
    """Run KMC for k-mer counting."""
    kmc_db = get_kmc_db_path(args)
    kmc_input_list = get_kmc_input_list_path(args)
    kmer_dir = get_kmer_dir(args)
    
    if not (os.path.exists(f"{kmc_db}.kmc_pre") and os.path.exists(f"{kmc_db}.kmc_suf")):
        logging.info(f"Running KMC: counting k-mers with k={args.kmer_size}...")
        cmd = [
            "kmc",
            f"-cs65535",
            f"-k{args.kmer_size}",
            f"-t{args.threads}",
            f"@{kmc_input_list}",
            kmc_db,
            kmer_dir
        ]
        logging.info("  " + " ".join(cmd))
        subprocess.run(cmd, check=True)
    else:
        logging.info(f"KMC database already exists: {kmc_db}")


def run_kmc_histogram(args):
    """Generate k-mer histogram from KMC database."""
    kmc_db = get_kmc_db_path(args)
    hist_file = get_kmc_hist_file(args)
    
    if not os.path.exists(hist_file):
        logging.info("Generating k-mer histogram...")
        cmd = [
            "kmc_tools",
            "transform",
            kmc_db,
            "histogram",
            hist_file
        ]
        logging.info("  " + " ".join(cmd))
        subprocess.run(cmd, check=True)
    else:
        logging.info(f"Histogram file already exists: {hist_file}")
    

def run_find_het_kmers(args):
    """Run heterozygous k-mer identification."""

    hist_file = get_kmc_hist_file(args)
    plot_dir = get_kmer_dir(args)
    existing_pngs = glob.glob(get_kmc_png_pattern(args))
    
    if not existing_pngs:
        logging.info("Running heterozygous k-mer identification...")
        cmd = [
            sys.executable,
            str(SCRIPT_DIR / "find_het_kmers.py"),
            hist_file,
            "-k", str(args.kmer_size),
            "-o", plot_dir
        ]
        
        # # Add optional parameters if specified
        # if args.est_kmercov:
        #     cmd.extend(["-est", args.est_kmercov])
        # if args.cutoff is not None:
        #     cmd.extend(["-cutoff", str(args.cutoff)])
        # if args.het is not None:
        #     cmd.extend(["-het", str(args.het)])
        # if args.rho is not None:
        #     cmd.extend(["-rho", str(args.rho)])
        
        logging.info("  " + " ".join(cmd))
        subprocess.run(cmd, check=True)
    else:
        logging.info(f"Plot files already exist: {existing_pngs}")


def parse_boundaries_from_pngs(args):
    """Extract left and right boundaries from PNG filenames."""

    logging.info(glob.glob(get_kmc_png_pattern(args)))
    png_files = sorted(glob.glob(get_kmc_png_pattern(args)))

    assert len(png_files) == 1
    
    # Extract boundaries from first PNG filename
    first_png = os.path.basename(png_files[0])
    stem = os.path.splitext(first_png)[0]
    parts = stem.split("_")

    assert len(parts) >= 2, "PNG filename does not contain enough parts to extract boundaries."
    left = parts[-2] if len(parts) >= 2 else None
    right = parts[-1] if len(parts) >= 1 else None
    
    return left, right

def run_fast2snp(args, left, right, cxx_args):
    """Run the core fast2snp variant caller."""
    snp_dir = get_snp_dir(args)
    kmc_db = get_kmc_db_path(args)
    bubble_snp = get_bubble_snp_file(args)


    if not os.path.exists(bubble_snp):
        logging.info("Running fast2snp variant calling...")
        logging.info(f"  Using cxx_args: {cxx_args}")
        cmd = [
            "kmc_dump",
            f"-ci{left}",
            kmc_db, 
            "/dev/stdout",
            " | ",
            str(SCRIPT_DIR / "fast2snp"),
            "-i", '-',
            "-l", str(left),
            "-r", str(right),
            "-t", str(args.threads),
            "-o", snp_dir + "/",
            " ".join(cxx_args)
        ]
        logging.info("  " + " ".join(cmd))
        subprocess.run(" ".join(cmd), shell=True, check=True)

    else:
        logging.info(f"SNP output already exists: {bubble_snp}")
    
def step1_count_kmers(args):
    """Step 1: KMC k-mer counting."""

    kmer_dir = get_kmer_dir(args)
    os.makedirs(kmer_dir, exist_ok=True)

    # Write KMC input list file
    kmc_input_list = get_kmc_input_list_path(args)
    with open(kmc_input_list, "w") as f:
        for file in args.input:
            f.write(f"{file}\n")
    
    run_kmc(args)
    run_kmc_histogram(args)
    run_find_het_kmers(args)

def step2_call_variants(args, cxx_args):
    """Step 2: Call variants using fast2snp."""

    snp_dir = get_snp_dir(args)
    os.makedirs(snp_dir, exist_ok=True)

    left, right = parse_boundaries_from_pngs(args)
    logging.info(f"Boundaries: left={left}, right={right}")

    run_fast2snp(args, left, right, cxx_args)
    logging.info(f"\nPipeline complete.")
def main():
    args, cxx_args = parse_known_args()
    
    logging.info(f"Script directory: {SCRIPT_DIR}")
    step1_count_kmers(args)
    step2_call_variants(args, cxx_args)


if __name__ == "__main__":
    main()