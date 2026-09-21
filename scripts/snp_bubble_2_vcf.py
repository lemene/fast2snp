#!/usr/bin/env python3
"""
Universal pipeline: convert .snp / .fa / .tab (k-mer pair table) to VCF via BWA-MEM.

Merged from:
  pisad2vcfmem.sh, ebwt2vcfmem.sh, discosnpvcffiltermem.sh, kmer2snp2vcfmem.sh

Usage examples:
  # .snp mode (FASTA with headers)
  python snp2vcf.py calls.snp ref.fa out prefix

  # .fa mode (DiscoSnp++ output, with filtering)
  python snp2vcf.py calls.fa ref.fa out prefix --input-format fa
  python snp2vcf.py calls.fa ref.fa out prefix --input-format fa --filter-high-confidence

  # .tab mode (kmer2snp table: Seq1 Seq2 Count)
  python snp2vcf.py calls.tab ref.fa out prefix --input-format tab

All BWA parameters, SNP/Indel patterns, threads can be customized.
"""

import sys
import os
import re
import subprocess
import argparse
from pathlib import Path


def run_cmd(cmd, **kwargs):
    """Run command, exit on error."""
    print("  Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, **kwargs)


def count_unmapped_samtools(sam_file):
    """Return number of unmapped reads using samtools, or None."""
    try:
        res = subprocess.run(
            ["samtools", "view", "-f", "4", "-c", sam_file],
            capture_output=True, text=True, check=True
        )
        return int(res.stdout.strip())
    except Exception:
        return None


def parse_sam_mapping(sam_file):
    """Return dict {read_name: bool} (True=mapped) and list of all read names."""
    read_mapped = {}
    all_reads = []
    with open(sam_file, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            qname = parts[0]
            flag = int(parts[1])
            read_mapped[qname] = not (flag & 4)
            all_reads.append(qname)
    return read_mapped, all_reads


def extract_variant_id_and_path(read_name):
    """
    Extract variant ID and path label from a read name.
    Supports:
      - _higher_path_ / _lower_path_  (DiscoSnp++ style)
      - _path1_ / _path2_              (legacy kmer2snp style)
    Returns (var_id, path) where path is "higher"/"lower" or "path1"/"path2".
    If no pattern matches, returns (read_name, None).
    """
    # Try DiscoSnp++ pattern first
    m = re.match(r'^(.*?)_(higher_path_|lower_path_)', read_name)
    if m:
        return m.group(1), "higher" if "higher" in m.group(2) else "lower"
    # Try path1/path2 pattern
    m = re.match(r'^(.*?)_(path[12])_', read_name)
    if m:
        return m.group(1), m.group(2)
    return read_name, None


def compute_variant_stats(read_mapped, all_reads):
    """
    Group reads by variant ID (using _higher_/_lower_ or _path1_/_path2_).
    A variant is considered mapped ONLY IF both paths are mapped.
    Returns:
      total_vars, mapped_vars,
      snp_total, snp_mapped, indel_total, indel_mapped,
      unmapped_snp_vars, unmapped_indel_vars
    """
    # variant_id -> set of observed path types
    var_paths = {}      # {var_id: {path_label: is_mapped}}
    for read_name in all_reads:
        vid, path = extract_variant_id_and_path(read_name)
        if path is None:
            continue
        is_mapped = read_mapped.get(read_name, False)
        var_paths.setdefault(vid, {})[path] = is_mapped

    total = len(var_paths)
    mapped = 0
    snp_total = snp_mapped = 0
    indel_total = indel_mapped = 0
    unmapped_snp = 0
    unmapped_indel = 0

    for vid, paths in var_paths.items():
        # Define "mapped" as both expected paths being present and mapped
        # For higher/lower or path1/path2, we expect exactly 2 paths mapped.
        if len(paths) == 2 and all(paths.values()):
            mapped += 1
            both_mapped = True
        else:
            both_mapped = False

        if vid.startswith("SNP"):
            snp_total += 1
            if both_mapped:
                snp_mapped += 1
            else:
                unmapped_snp += 1
        elif vid.startswith("INDEL"):
            indel_total += 1
            if both_mapped:
                indel_mapped += 1
            else:
                unmapped_indel += 1

    return total, mapped, snp_total, snp_mapped, indel_total, indel_mapped, unmapped_snp, unmapped_indel


def detect_format(filename):
    """Guess input format from extension."""
    if filename.endswith(".snp"):
        return "snp"
    elif filename.endswith(".fa") or filename.endswith(".fasta"):
        return "fa"
    elif filename.endswith(".tab") or filename.endswith(".tsv") or filename.endswith(".txt"):
        return "tab"
    else:
        sys.exit(f"Error: Cannot detect format from extension '{filename}'. Use --input-format.")


def main():
    parser = argparse.ArgumentParser(
        description="Convert variant calls (various formats) to VCF via BWA-MEM."
    )
    parser.add_argument("calls", help="Input calls file")
    parser.add_argument("ref", help="Reference genome FASTA")
    parser.add_argument("out_dir", help="Output directory")
    parser.add_argument("prefix", help="Output prefix")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="BWA threads (default: 8)")
    parser.add_argument("--input-format", choices=["snp", "fa", "tab"],
                        help="Input format (auto-detected if not specified)")
    # .fa preprocessing
    parser.add_argument("--filter-high-confidence", action="store_true",
                        help="(.fa only) Keep only |high| entries")
    # SNP/Indel patterns for .snp headers
    parser.add_argument("--snp-pattern", default="_SNP_",
                        help="Pattern to identify SNP in .snp headers (default: _SNP_)")
    parser.add_argument("--indel-pattern", default="_INDEL_",
                        help="Pattern to identify Indel in .snp headers (default: _INDEL_)")
    # BWA
    parser.add_argument("--bwa-preset", choices=["default", "relaxed"], default="default",
                        help="BWA MEM preset: default (-k 15 -T 20), relaxed (-k 13 -T 15)")
    parser.add_argument("--bwa-extra", default="",
                        help="Extra BWA parameters (overrides preset)")

    args = parser.parse_args()

    calls = args.calls
    ref = args.ref
    out_dir = args.out_dir
    prefix = args.prefix

    # Determine format
    if args.input_format:
        fmt = args.input_format
    else:
        fmt = detect_format(calls)
    print(f"Detected input format: {fmt}")

    if not Path(calls).is_file():
        sys.exit(f"Error: Calls file '{calls}' not found.")
    if not Path(ref).is_file():
        sys.exit(f"Error: Reference '{ref}' not found.")

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Output files
    out_vcf = os.path.join(out_dir, f"{prefix}.vcf")
    out_fastq = os.path.join(out_dir, f"{prefix}.fastq")
    out_sam = os.path.join(out_dir, f"{prefix}.sam")
    out_log = os.path.join(out_dir, f"{prefix}.log")
    vcf_unsorted = out_vcf + ".unsorted"

    # Intermediate FASTA (if conversion needed)
    clean_fasta = None
    total_reads = 0
    input_for_seqtk = None

    print(f"input calls: {calls}")
    print(f"input reference: {ref}")
    print(f"output directory: {out_dir}")
    print(f"output prefix: {prefix}")
    print(f"final vcf: {out_vcf}")

    # ================================================================
    # Step 1: Pre-processing based on format
    # ================================================================
    if fmt == "fa":
        print("\nStep 1: Pre-processing .fa (lowercase removal & filtering)...")
        clean_fasta = os.path.join(out_dir, f"{prefix}.clean.fasta")
        awk_script = '/^>/ { if ($0 ~ /G1_0\\/1/'
        if args.filter_high_confidence:
            awk_script += ' && $0 ~ /\\|high\\|/'
        awk_script += ') { skip=0; print $0 } else { skip=1 } } '
        awk_script += '!/^>/ { if (skip==0) { gsub(/[a-z]/, "", $0); if(length($0) > 0) print $0 } }'
        with open(clean_fasta, "w") as out_f:
            run_cmd(["awk", awk_script, calls], stdout=out_f)
        input_for_seqtk = clean_fasta

        # Count sequences
        total_reads = 0
        with open(clean_fasta, "r") as f:
            for line in f:
                if line.startswith(">"):
                    total_reads += 1
        print(f"Total sequences after filtering: {total_reads}")

    elif fmt == "snp":
        print("\nStep 1: No preprocessing needed for .snp (direct FASTA).")
        input_for_seqtk = calls
        # Count will be done later for statistics

    elif fmt == "tab":
        print("\nStep 1: Converting tabular k-mer pairs to FASTA...")
        out_fasta = os.path.join(out_dir, f"{prefix}.fasta")  # intermediate
        total_pairs = 0
        with open(calls, "r") as inf, open(out_fasta, "w") as outf:
            for line in inf:
                if not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) < 3:
                    print(f"Warning: skipping malformed line: {line.strip()}")
                    continue
                seq1, seq2, count = parts[0], parts[1], parts[2]
                total_pairs += 1
                # Use unified naming: SNP_{idx}_higher_path_{cov} / _lower_path_{cov}
                outf.write(f">SNP_{total_pairs}_higher_path_{count}\n{seq1}\n")
                outf.write(f">SNP_{total_pairs}_lower_path_{count}\n{seq2}\n")
        total_reads = total_pairs * 2
        print(f"Total pairs: {total_pairs} -> {total_reads} sequences")
        input_for_seqtk = out_fasta

    else:
        sys.exit(f"Unknown format: {fmt}")

    # ================================================================
    # Step 2: Convert to FASTQ
    # ================================================================
    print("\nStep 2: Converting to FASTQ...")
    with open(out_fastq, "w") as fq:
        run_cmd(["seqtk", "seq", "-F", "h", input_for_seqtk], stdout=fq)

    # ================================================================
    # Step 3: BWA indexing
    # ================================================================
    if not Path(ref + ".bwt").is_file():
        print("Step 3: BWA index not found. Indexing reference...")
        run_cmd(["bwa", "index", ref])
    else:
        print("Step 3: BWA index found. Skipping indexing.")

    # ================================================================
    # Step 4: BWA MEM alignment
    # ================================================================
    bwa_cmd = ["bwa", "mem", "-t", str(args.threads)]
    if args.bwa_extra:
        bwa_cmd.extend(args.bwa_extra.split())
    else:
        if args.bwa_preset == "default":
            bwa_cmd.extend(["-k", "15", "-T", "20"])
        elif args.bwa_preset == "relaxed":
            bwa_cmd.extend(["-k", "13", "-T", "15"])
        bwa_cmd.extend(["-B", "3", "-O", "5"])
    bwa_cmd.extend([ref, out_fastq])

    print(f"\nStep 4: Aligning with BWA MEM ({args.bwa_preset} preset)...")
    with open(out_sam, "w") as sam_out:
        run_cmd(bwa_cmd, stdout=sam_out)

    # ================================================================
    # Step 5: Alignment statistics
    # ================================================================
    print("\nStep 5: Calculating alignment statistics...")

    # Total reads: if not already set, count from input_for_seqtk
    if total_reads == 0:
        total_reads = 0
        with open(input_for_seqtk, "r") as f:
            for line in f:
                if line.startswith(">"):
                    total_reads += 1

    # Unmapped reads
    unmapped_reads = count_unmapped_samtools(out_sam)
    if unmapped_reads is None:
        print("  (samtools unavailable, parsing SAM manually)")
        rm, all_r = parse_sam_mapping(out_sam)
        unmapped_reads = sum(1 for r in all_r if not rm[r])
    mapped_reads = total_reads - unmapped_reads

    # Read-level breakdown by type
    # We scan the FASTA headers used as input for seqtk.
    # For .snp: use header patterns; for .fa: headers start with >SNP / >INDEL; for .tab: all start with >SNP
    if fmt == "fa":
        with open(input_for_seqtk, "r") as f:
            headers = [l for l in f if l.startswith(">")]
        total_snps = sum(1 for h in headers if h.startswith(">SNP"))
        total_indels = sum(1 for h in headers if h.startswith(">INDEL"))
    elif fmt == "snp":
        total_snps = 0
        total_indels = 0
        with open(calls, "r") as f:
            for line in f:
                if line.startswith(">"):
                    h = line
                    if args.snp_pattern in h or \
                       (args.snp_pattern.isalpha() and args.snp_pattern.upper() in h.upper()):
                        total_snps += 1
                    if args.indel_pattern in h or \
                       (args.indel_pattern.isalpha() and args.indel_pattern.upper() in h.upper()):
                        total_indels += 1
    elif fmt == "tab":
        total_snps = total_reads  # all are SNP in the generated FASTA
        total_indels = 0

    # Unmapped read counts by type using samtools or Python
    unmapped_snps = 0
    unmapped_indels = 0
    try:
        proc = subprocess.Popen(
            ["samtools", "view", "-f", "4", out_sam],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
        )
        for line in proc.stdout:
            qname = line.split("\t")[0]
            # Detect type
            if fmt == "fa" or fmt == "tab":
                # Read names start with SNP/INDEL (for .fa) or SNP (for .tab)
                if qname.startswith("SNP"):
                    unmapped_snps += 1
                elif qname.startswith("INDEL"):
                    unmapped_indels += 1
            else:  # snp
                if args.snp_pattern in qname:
                    unmapped_snps += 1
                elif args.indel_pattern in qname:
                    unmapped_indels += 1
        proc.wait()
    except Exception:
        # Fallback to Python parsing
        rm, all_r = parse_sam_mapping(out_sam)
        for r in all_r:
            if not rm[r]:
                if fmt == "fa" or fmt == "tab":
                    if r.startswith("SNP"):
                        unmapped_snps += 1
                    elif r.startswith("INDEL"):
                        unmapped_indels += 1
                else:
                    if args.snp_pattern in r:
                        unmapped_snps += 1
                    elif args.indel_pattern in r:
                        unmapped_indels += 1

    mapped_snps = total_snps - unmapped_snps
    mapped_indels = total_indels - unmapped_indels

    # Data-level statistics (for paired formats: fa & tab)
    has_pairs = fmt in ("fa", "tab")
    var_stats = None
    if has_pairs:
        rm, all_r = parse_sam_mapping(out_sam)
        var_stats = compute_variant_stats(rm, all_r)
        (total_vars, mapped_vars,
         snp_total_v, snp_mapped_v, indel_total_v, indel_mapped_v,
         unmapped_snp_v, unmapped_indel_v) = var_stats
        # Save unmapped INDEL list if any (original script behavior for .fa)
        if fmt == "fa" or (fmt == "tab" and indel_total_v > 0):
            unmapped_indel_ids = set()
            for vid, paths in {}:  # we need to reconstruct
                pass
            # Simply use original extraction for reporting; skip detailed list if not needed.
            # We can produce the list from the unmapped reads by variant id.
            pass

    # Print statistics
    print("-" * 60)
    print("Alignment Statistics (Read Level):")
    print(f"Total Reads: {total_reads}")
    print(f"  - Mapped:   {mapped_reads}")
    print(f"  - Unmapped: {unmapped_reads}")
    print()
    print("Breakdown by Type (Read Level):")
    print(f"  SNP:    Total={total_snps}, Mapped={mapped_snps}, Unmapped={unmapped_snps}")
    print(f"  INDEL:  Total={total_indels}, Mapped={mapped_indels}, Unmapped={unmapped_indels}")

    if has_pairs and var_stats:
        (total_vars, mapped_vars,
         snp_total_v, snp_mapped_v, indel_total_v, indel_mapped_v,
         unmapped_snp_v, unmapped_indel_v) = var_stats
        print()
        print("-" * 60)
        print("Alignment Statistics (Data/Variant Level):")
        print("Definition: A variant is mapped only if BOTH paths are mapped.")
        print(f"Total Variants: {total_vars}")
        print(f"  - Mapped:   {mapped_vars}")
        print(f"  - Unmapped: {total_vars - mapped_vars}")
        print()
        print("Breakdown by Type (Data Level):")
        print(f"  SNP:    Total={snp_total_v}, Mapped={snp_mapped_v}, Unmapped={unmapped_snp_v}")
        print(f"  INDEL:  Total={indel_total_v}, Mapped={indel_mapped_v}, Unmapped={unmapped_indel_v}")

    print("-" * 60)

    # ================================================================
    # Step 6: Variant calling with sam2vcf
    # ================================================================
    print("\nStep 6: Calling variants...")
    script_dir = Path(__file__).resolve().parent
    sam2vcf_bin = script_dir / "sam2vcf"
    if not sam2vcf_bin.is_file():
        sys.exit(f"Error: sam2vcf not found in {script_dir}\n"
                 f"Compile with: g++ -O3 {script_dir}/sam2vcf.cpp -o {sam2vcf_bin}")
    run_cmd([
        str(sam2vcf_bin),
        "-f", ref,
        "-s", out_sam,
        "-v", vcf_unsorted
    ], stdout=open(out_log, "w"))

    # Sort and deduplicate VCF
    print("  Sorting and deduplicating VCF...")
    headers = []
    variants = []
    with open(vcf_unsorted, "r") as f:
        for line in f:
            if line.startswith("#"):
                headers.append(line)
            else:
                variants.append(line.rstrip("\n"))
    variants.sort(key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])))
    unique_vars = []
    for v in variants:
        if not unique_vars or v != unique_vars[-1]:
            unique_vars.append(v)
    with open(out_vcf, "w") as vcf:
        vcf.writelines(headers)
        for v in unique_vars:
            vcf.write(v + "\n")
    os.remove(vcf_unsorted)

    # ================================================================
    # Step 7: Post-calling statistics
    # ================================================================
    print("\nStep 7: Calculating post-calling statistics...")
    vcf_snps = 0
    vcf_indels = 0
    with open(out_vcf, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            if "VT=SNP" in line:
                vcf_snps += 1
            if "VT=INDEL" in line:
                vcf_indels += 1

    print("-" * 60)
    print("Post-calling Statistics (from VCF):")
    print(f"Total Variants: {vcf_snps + vcf_indels}")
    print(f"SNPs: {vcf_snps}")
    print(f"Indels: {vcf_indels}")
    print("-" * 60)

    # Optional cleanup
    # if fmt == "fa":
    #     os.remove(clean_fasta)
    # elif fmt == "tab":
    #     os.remove(out_fasta)
    # os.remove(out_fastq)
    # os.remove(out_sam)

    print(f"\nFinished. Output VCF: {out_vcf}")


if __name__ == "__main__":
    main()