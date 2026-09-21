#!/usr/bin/env python3
"""
Test fast2snp performance across different k-mer sizes.

For each k, runs fast2snp → hap.py benchmark, collects runtime, memory,
and PASS recall/precision/F1 for SNPs and INDELs.
Outputs two summary tables: kmer_performance_summary.txt and kmer_pass_metrics_summary.csv.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import time
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Get the directory where this script is located
SCRIPT_DIR = Path(__file__).resolve().parent
RUN_FAST2SNP_SCRIPT = SCRIPT_DIR / "run_fast2snp.py"
RUN_VALIDATION_SCRIPT = SCRIPT_DIR / "run_validation.py"

def parse_args():
    parser = argparse.ArgumentParser(description="Batch test fast2snp over multiple k-mer sizes")

    parser.add_argument("input", nargs="+", help="Input files (space-separated *.fastq or *.fastq.gz)")
    parser.add_argument("-k", "--kmer-size", type=int, required=True, default=31)
    parser.add_argument("--covs", required=True, type=str, 
                        help="Coverage values (comma-separated)")
    parser.add_argument("--truth-vcf", required=True,
                        help="Truth VCF (e.g. HG001_chr22_het.vcf.gz)")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--bed", required=True,
                        help="High-confidence BED file")
    parser.add_argument("-o", "--outdir", default=".", help="Root directory for all per‑k outputs")
    parser.add_argument("--threads", type=int, default=8)


    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing")
    args = parser.parse_args()
    
    return parser.parse_args()

def get_output_dir(args, c):
    return Path(args.outdir) / ("cov_" + str(c))

def get_reads_file(args, c):
    return get_output_dir(args, c) / "reads.fasta"


def get_bubble_snp_file(args, k):
    return get_output_dir(args, k) / "variant_calls" / f"bubbles.snp"

def get_happy_prefix(args, k):
    return get_output_dir(args, k) / "validation_benchmark" / "variants"
def get_happy_summary_file(args, k):
    return get_output_dir(args, k) / "validation_benchmark" / f"variants.summary.csv"

def get_results_csv_path(args):
    return Path(args.outdir) / "results.csv"


def get_results_png_path(args):
    return Path(args.outdir) / "results.png"

def run_cmd(cmd, dry_run=False):
    """Print and optionally execute a command."""
    print("  Running:", cmd)
    print("  Running:", " ".join(cmd))
    if not dry_run:
        subprocess.run(cmd, check=True)
def parse_time_mem(time_log: Path):
    """Extract wall‑clock time and maximum resident set size from /usr/bin/time -v output."""
    time_str = mem_str = "NA"
    if time_log.is_file():
        text = time_log.read_text()
        # Elapsed (wall clock) time (h:mm:ss or m:ss)
        m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)", text)
        if m:
            time_str = m.group(1).strip()
        # Maximum resident set size (kbytes)
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
        if m:
            mem_str = m.group(1)
    return time_str, mem_str


def extract_pass_metrics(summary_csv: Path):
    """
    Read hap.py summary CSV, return a dict:
    { "SNP": (recall, precision, f1), "INDEL": (recall, precision, f1) }
    Values are "NA" if not found.
    """
    metrics = {
        "SNP": ("NA", "NA", "NA"),
        "INDEL": ("NA", "NA", "NA"),
    }
    if not summary_csv.is_file():
        return metrics

    with open(summary_csv) as f:
        reader = csv.reader(f)
        rows = list(reader)
    if len(rows) < 2:
        return metrics

    header = [h.strip().upper() for h in rows[0]]
    data = rows[1:]

    # Column indices
    type_idx = filter_idx = recall_idx = precision_idx = f1_idx = None
    for i, h in enumerate(header):
        if h in {"TYPE", "VAR_TYPE", "NENGB"}:
            type_idx = i
        elif h in {"FILTER"}:
            filter_idx = i
        elif h in {"METRIC.RECALL"}:
            recall_idx = i
        elif h in {"METRIC.PRECISION"}:
            precision_idx = i
        elif h in {"METRIC.F1_SCORE"}:
            f1_idx = i

    # Fallback: try to guess columns from content if header is unusual
    if any(idx is None for idx in [type_idx, filter_idx, recall_idx, precision_idx, f1_idx]):
        max_cols = max(len(r) for r in data) if data else 0
        for i in range(max_cols):
            vals = {r[i].strip().upper() for r in data if i < len(r)}
            if type_idx is None and ("SNP" in vals or "INDEL" in vals):
                type_idx = i
            if filter_idx is None and "PASS" in vals:
                filter_idx = i

    if any(idx is None for idx in [type_idx, filter_idx, recall_idx, precision_idx, f1_idx]):
        return metrics

    for row in data:
        vtype = row[type_idx].strip().upper()
        flt = row[filter_idx].strip().upper()
        if (vtype in ("SNP", "INDEL")) and flt == "PASS":
            rec = row[recall_idx].strip() if recall_idx < len(row) else "NA"
            prec = row[precision_idx].strip() if precision_idx < len(row) else "NA"
            f1 = row[f1_idx].strip() if f1_idx < len(row) else "NA"
            metrics[vtype] = (rec, prec, f1)

    return metrics


def get_seq_size(filepath):
    """Return (num_records, total_bases) for FASTA/FASTQ."""
    suffix = Path(filepath).suffix
    fmt = "fastq" if suffix in (".fq", ".fastq", ".fq.gz", ".fastq.gz") else "fasta"
    n = total = 0
    from Bio import SeqIO
    for rec in SeqIO.parse(filepath, fmt):
        n += 1
        total += len(rec)
    return n, total


import subprocess
import csv
import io

def total_bases_seqkit(filepath: str, threads: int = 2) -> int:
    """使用 seqkit stats 快速获取总碱基数"""
    cmd = ["seqkit", "stats", "-T", "-j", str(threads), filepath]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    # 解析 TSV 输出，第二列是 sum_len
    reader = csv.DictReader(io.StringIO(result.stdout), delimiter='\t')
    row = next(reader)
    return int(row['sum_len'])

def total_bases_multi_using_seqkit(files, threads=4):
    """统计多个文件的总碱基数（使用 seqkit）"""
    cmd = ["seqkit", "stats", "-T", "-j", str(threads)] + files
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    reader = csv.DictReader(io.StringIO(result.stdout), delimiter='\t')
    total = sum(int(row['sum_len']) for row in reader)
    return total

def subsample_reads_by_depth(args, cov):
    logging.info(f"Subsampling reads by depth: {cov}")


    reads = get_reads_file(args, cov)
    if not os.path.exists(reads) or os.path.getsize(reads) == 0:
        _, ref_size = get_seq_size(args.ref)
        sample_size = ref_size * cov


        input_size = total_bases_multi_using_seqkit(args.input)
        p = sample_size / input_size
        logging.info(f"ref size: {ref_size}, sample size: {sample_size}, input size: {input_size}, p: {p}")

        cmd = ["seqkit", "sample", f"-s{11}", "-p", str(p)] + args.input
        with open(reads, "w") as f_out:
            subprocess.run(cmd, stdout=f_out, check=True)



def run_fast2snp(args, cov, koutdir):
    print(f"Running fast2snp for Coverage={cov}")

    
    reads = get_reads_file(args, cov)
    run_fast2snp_cmd = [
        sys.executable, str(RUN_FAST2SNP_SCRIPT),
        str(reads),
        "-o", str(koutdir),
        "-k", str(args.kmer_size),
        "-C", "10000",
        "--threads", str(args.threads),
    ]
    print(run_fast2snp_cmd)
    run_cmd(run_fast2snp_cmd, dry_run=args.dry_run)


def run_validation(args, k, koutdir):
    print(f"Running validation for k={k}")

    
    run_validation_cmd = [
        sys.executable, str(RUN_VALIDATION_SCRIPT),
        "--snp", str(get_bubble_snp_file(args, k)),
        "-o", str(get_output_dir(args, k)),
        "--truth", args.truth_vcf,
        "--ref", args.ref,
        "--bed", args.bed,
        "--threads", str(args.threads),
    ]
    run_cmd(run_validation_cmd, dry_run=args.dry_run)

def collect_metrics(args, k, koutdir):
    """Collect wall time, memory, and PASS metrics from hap.py summary."""
    summary_csv = get_happy_summary_file(args, k)

    metrics = extract_pass_metrics(summary_csv)
    return metrics


import pandas as pd
import matplotlib.pyplot as plt
def plot_performance_from_csv(csv_path, png_path):
    """
    从 CSV 读取数据并绘制 SNP/INDEL 性能对比图
    """
    # 1. 读取 CSV 数据
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"错误: 找不到文件 '{csv_path}'，请检查路径。")
        return
    except Exception as e:
        print(f"读取 CSV 时出错: {e}")
        return

    # 2. 设置绘图风格 (支持中文显示，防止乱码)
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans'] 
    plt.rcParams['axes.unicode_minus'] = False
    
    # 创建 1x2 的子图布局，共享 Y 轴以便对比
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))#, sharey=True)
    fig.suptitle('Performance of SNP and INDEL Detection under Different k-mer Lengths', fontsize=16, fontweight='bold')


    ax1.autoscale(axis='y')        # 左轴根据 y1 自动调整
    ax2.autoscale(axis='y')        # 右轴根据 y2 自动调整

    # --- 子图 (a): SNP Performance ---
    # 绘制召回率 (Recall) - 黑色方块
    ax1.plot(df['Coverage'], df['SNP_PASS_Recall'], marker='s', color='black', linestyle='-', label='SNP Recall')
    # 绘制精确率 (Precision) - 红色圆点
    ax1.plot(df['Coverage'], df['SNP_PASS_Precision'], marker='o', color='#d62728', linestyle='-', label='SNP Precision')
    ax1.plot(df['Coverage'], df['SNP_PASS_F1'], marker='^', color='#1f77b4', 
            linestyle='--', label='SNP F1')
    ax1.set_title('(a) SNP Detection Performance', fontsize=14)
    ax1.set_xlabel('k-mer Length', fontsize=12)
    ax1.set_ylabel('Recall or Precision', fontsize=12)
    ax1.legend(loc='lower right', fontsize=10)
    ax1.grid(True, linestyle='--', alpha=0.5)

    # --- 子图 (b): INDEL Performance ---
    # 绘制召回率 (Recall) - 黑色方块
    ax2.plot(df['Coverage'], df['INDEL_PASS_Recall'], marker='s', color='black', linestyle='-', label='INDEL Recall')
    # 绘制精确率 (Precision) - 红色圆点
    ax2.plot(df['Coverage'], df['INDEL_PASS_Precision'], marker='o', color='#d62728', linestyle='-', label='INDEL Precision')
    ax2.plot(df['Coverage'], df['INDEL_PASS_F1'], marker='^', color='#1f77b4', 
        linestyle='--', label='INDEL F1')
    ax2.set_title('(b) INDEL Detection Performance', fontsize=14)
    ax2.set_xlabel('k-mer Length', fontsize=12)
    ax2.legend(loc='lower right', fontsize=10)
    ax2.grid(True, linestyle='--', alpha=0.5)

    # 调整布局并保存
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"图片已保存至: {png_path}")
    plt.show()

def main():
    args = parse_args()

    mf = open(get_results_csv_path(args), "w", newline="")
    writer = csv.writer(mf)
    writer.writerow(["Coverage",
                        "SNP_PASS_Recall", "SNP_PASS_Precision", "SNP_PASS_F1",
                        "INDEL_PASS_Recall", "INDEL_PASS_Precision", "INDEL_PASS_F1"])

    covs = [int(c) for c in args.covs.split(",")]
    for c in covs:
        coutdir = get_output_dir(args, c)
        coutdir.mkdir(parents=True, exist_ok=True)
        subsample_reads_by_depth(args, c)
        run_fast2snp(args, c, coutdir)
        run_validation(args, c, coutdir)
        metrics = collect_metrics(args, c, coutdir)

        snp_rec, snp_prec, snp_f1 = metrics["SNP"]
        indel_rec, indel_prec, indel_f1 = metrics["INDEL"]
        print(f"Coverage={c}, SNP: {snp_rec}, {snp_prec}, {snp_f1}, INDEL: {indel_rec}, {indel_prec}, {indel_f1}")
        writer.writerow([c,
                            snp_rec, snp_prec, snp_f1,
                            indel_rec, indel_prec, indel_f1])

    mf.close()
    plot_performance_from_csv(get_results_csv_path(args), get_results_png_path(args))

if __name__ == "__main__":
    main()