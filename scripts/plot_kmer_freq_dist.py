#!/usr/bin/env python3
"""
kmer_freq_count.py

统计指定 k-mer 集合中不同频次的 k-mer 数量。
输出两列：frequency  count
只输出频次 > 0 的项，按频次升序排列。

用法:
  python kmer_freq_count.py \
      --freq kmer_counts.txt \
      --kmers het_kmers.fasta \
      --output freq_count.txt

可选:
  --include-missing-as-zero  将频次文件中缺失的 k-mer 频次计为 0（默认忽略缺失）
"""

import argparse
import sys
from pathlib import Path
from collections import Counter

def read_freq_file(freq_file):
    """读取频次文件，返回字典 {kmer: count}"""
    freq = {}
    with open(freq_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                kmer, count = parts[0], int(parts[1])
                freq[kmer] = count
    return freq

def read_kmer_set(kmer_file):
    """
    读取 k-mer 集合，返回列表。
    自动识别 FASTA 格式（首行以 > 开头）或纯文本（每行一个 k-mer）。
    """
    kmers = []
    with open(kmer_file, 'r') as f:
        first_line = f.readline().strip()
        f.seek(0)
        if first_line.startswith('>'):
            # FASTA 格式
            seq = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if seq:
                        kmers.append(''.join(seq))
                        seq = []
                else:
                    seq.append(line)
            if seq:
                kmers.append(''.join(seq))
        else:
            # 纯文本
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                kmers.append(parts[0])
    return kmers

def main():
    parser = argparse.ArgumentParser(description="统计 k-mer 集合的频次分布")
    parser.add_argument("--freq", required=True, help="k-mer 频次文件")
    parser.add_argument("--kmers", required=True, help="k-mer 集合文件（纯文本或 FASTA）")
    parser.add_argument("--output", default=None, help="输出文件（默认 stdout）")
    parser.add_argument("--include-missing-as-zero", action="store_true",
                        help="将频次文件中缺失的 k-mer 频次计为 0（默认忽略缺失）")
    args = parser.parse_args()

    if not Path(args.freq).is_file():
        sys.exit(f"错误: 频次文件不存在: {args.freq}")
    if not Path(args.kmers).is_file():
        sys.exit(f"错误: k-mer 集合文件不存在: {args.kmers}")

    # 读取频次字典
    print("读取频次文件...", file=sys.stderr)
    freq = read_freq_file(args.freq)
    print(f"  频次文件中包含 {len(freq)} 个唯一 k-mer", file=sys.stderr)

    # 读取 k-mer 集合
    print("读取 k-mer 集合...", file=sys.stderr)
    kmer_list = read_kmer_set(args.kmers)
    print(f"  集合中共有 {len(kmer_list)} 个 k-mer", file=sys.stderr)

    # 去重
    unique_kmers = list(dict.fromkeys(kmer_list))
    print(f"  去重后 {len(unique_kmers)} 个", file=sys.stderr)

    # 统计频次分布
    freq_counter = Counter()
    missing = 0
    for kmer in unique_kmers:
        if kmer in freq:
            freq_counter[freq[kmer]] += 1
        else:
            if args.include_missing_as_zero:
                freq_counter[0] += 1
            else:
                missing += 1

    if missing > 0 and not args.include_missing_as_zero:
        print(f"  忽略 {missing} 个频次文件中缺失的 k-mer", file=sys.stderr)

    # 按频次升序输出，只输出频次 > 0 的项
    out_fh = open(args.output, 'w') if args.output else sys.stdout
    out_fh.write("frequency\tcount\n")
    for freq_value in sorted(freq_counter.keys()):
        if freq_value == 0 and not args.include_missing_as_zero:
            continue  # 如果没有包含0，则跳过
        if freq_value <= 0:
            continue
        out_fh.write(f"{freq_value}\t{freq_counter[freq_value]}\n")

    if args.output:
        out_fh.close()

    print("完成。", file=sys.stderr)

if __name__ == "__main__":
    main()