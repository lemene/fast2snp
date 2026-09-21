#!/usr/bin/env python3
"""
query_read_kmers.py

查询 reads 文件中每个 k-mer 在频率文件中的频次。
频率文件中的 k-mer 为标准形式（canonical，即与反向互补序列取字典序较小者）。
查询时也会对 read 中的 k-mer 进行同样的标准化。

用法:
  python query_read_kmers.py \
      --freq kmer_counts.txt \
      --reads reads.fastq.gz \
      --k 31 \
      --output read_kmer_freq.txt

参数:
  --freq      k-mer 频率文件，每行格式: kmer count (空格或制表符分隔)
  --reads     输入 reads 文件（FASTA/FASTQ，支持 .gz）
  --k         k-mer 长度（若不指定，自动从频率文件推断）
  --output    输出文件（默认 stdout）
  --missing   设置 k-mer 缺失时的频次值（默认 0）
"""

import argparse
import sys
import gzip
from pathlib import Path

def revcomp(seq):
    """返回反向互补序列"""
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(comp.get(base, 'N') for base in reversed(seq))

def canonical(kmer):
    """返回 k-mer 的标准形式（与反向互补序列取字典序较小者）"""
    rc = revcomp(kmer)
    return min(kmer, rc)

def read_freq_file(freq_file):
    """读取频率文件，返回 (dict, k)，键为标准 k-mer"""
    freq = {}
    k = None
    with open(freq_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                kmer = parts[0].upper()
                count = int(parts[1])
                # 标准化
                canon = canonical(kmer)
                freq[canon] = freq.get(canon, 0) + count  # 若重复累加
                if k is None:
                    k = len(kmer)
    return freq, k

def parse_reads(reads_file):
    """解析 FASTA/FASTQ，yield (read_name, sequence)"""
    opener = gzip.open if reads_file.endswith('.gz') else open
    with opener(reads_file, 'rt') as f:
        first_line = f.readline()
        if not first_line:
            return
        if first_line.startswith('>'):
            # FASTA
            f.seek(0)
            name = None
            seq_parts = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if name is not None:
                        yield name, ''.join(seq_parts)
                    name = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if name is not None:
                yield name, ''.join(seq_parts)
        elif first_line.startswith('@'):
            # FASTQ
            f.seek(0)
            while True:
                name = f.readline().strip()
                if not name:
                    break
                if not name.startswith('@'):
                    continue
                seq = f.readline().strip()
                plus = f.readline().strip()
                qual = f.readline().strip()
                yield name[1:].split()[0], seq
        else:
            raise ValueError("无法识别 reads 文件格式（需为 FASTA 或 FASTQ）")

def main():
    parser = argparse.ArgumentParser(description="查询 reads 中每个 k-mer 的频次")
    parser.add_argument("--freq", required=True, help="k-mer 频率文件")
    parser.add_argument("--reads", required=True, help="reads 文件（FASTA/FASTQ）")
    parser.add_argument("--k", type=int, default=None, help="k-mer 长度（默认自动推断）")
    parser.add_argument("--output", default=None, help="输出文件（默认 stdout）")
    parser.add_argument("--missing", type=int, default=0,
                        help="k-mer 缺失时的频次值（默认 0）")
    args = parser.parse_args()

    if not Path(args.freq).is_file():
        sys.exit(f"错误: 频率文件不存在: {args.freq}")
    if not Path(args.reads).is_file():
        sys.exit(f"错误: reads 文件不存在: {args.reads}")

    # 读取频率文件（键已标准化）
    print("加载频率文件...", file=sys.stderr)
    freq, inferred_k = read_freq_file(args.freq)
    if args.k is None:
        if inferred_k is None:
            sys.exit("无法从频率文件推断 k，请用 --k 指定")
        args.k = inferred_k
        print(f"  推断 k = {args.k}", file=sys.stderr)
    else:
        print(f"  使用 k = {args.k}", file=sys.stderr)
    print(f"  频率文件包含 {len(freq)} 个唯一标准 k-mer", file=sys.stderr)

    # 输出文件
    out_fh = open(args.output, 'w') if args.output else sys.stdout
    out_fh.write("read_name\tposition\tkmer\tfrequency\n")

    total_kmers = 0
    missing_kmers = 0

    # 遍历 reads
    for read_name, seq in parse_reads(args.reads):
        seq = seq.upper()
        if len(seq) < args.k:
            continue
        for i in range(len(seq) - args.k + 1):
            raw_kmer = seq[i:i+args.k]
            # 跳过含 N 的 k-mer
            if 'N' in raw_kmer:
                continue
            # 标准化
            canon = canonical(raw_kmer)
            count = freq.get(canon, args.missing)
            if canon not in freq:
                missing_kmers += 1
            out_fh.write(f"{read_name}\t{i+1}\t{canon}\t{count}\n")
            total_kmers += 1

    if args.output:
        out_fh.close()

    print(f"完成：查询了 {total_kmers} 个 k-mer，其中 {missing_kmers} 个在频率文件中缺失", file=sys.stderr)

if __name__ == "__main__":
    main()