#!/usr/bin/env python3
"""
compare_first_kmers.py

比较两个 FASTA 文件中每个序列前 k 个碱基构成的 k-mer 集合。

用法:
  python compare_first_kmers.py --fasta1 sample1.fa --fasta2 sample2.fa --k 31 [--output result.txt]

输出:
  集合大小、交集大小、仅在文件1中的数量、仅在文件2中的数量，并可选保存详细列表。
"""

import argparse
import gzip
from pathlib import Path

def read_first_kmers(fasta_file, k):
    """
    读取 FASTA 文件，提取每个序列的前 k 个碱基作为 k-mer。
    返回集合。
    支持 .gz 压缩文件。
    """
    kmers = set()
    opener = gzip.open if fasta_file.endswith('.gz') else open

    with opener(fasta_file, 'rt') as f:
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    # 处理上一个序列
                    full_seq = ''.join(seq)
                    if len(full_seq) >= k:
                        kmers.add(full_seq[:k].upper())
                    seq = []
            else:
                seq.append(line)
        # 最后一个序列
        if seq:
            full_seq = ''.join(seq)
            if len(full_seq) >= k:
                kmers.add(full_seq[:k].upper())
    return kmers

def main():
    parser = argparse.ArgumentParser(description="比较两个 FASTA 文件的第一个 k-mer 集合")
    parser.add_argument("--fasta1", required=True, help="第一个 FASTA 文件（支持 .gz）")
    parser.add_argument("--fasta2", required=True, help="第二个 FASTA 文件（支持 .gz）")
    parser.add_argument("--k", type=int, required=True, help="k-mer 长度")
    parser.add_argument("--output", default=None, help="输出文件（默认 stdout）")
    args = parser.parse_args()

    if not Path(args.fasta1).is_file():
        raise FileNotFoundError(f"文件不存在: {args.fasta1}")
    if not Path(args.fasta2).is_file():
        raise FileNotFoundError(f"文件不存在: {args.fasta2}")

    print("读取第一个文件...")
    set1 = read_first_kmers(args.fasta1, args.k)
    print(f"  文件1: {len(set1)} 个唯一 k-mer")

    print("读取第二个文件...")
    set2 = read_first_kmers(args.fasta2, args.k)
    print(f"  文件2: {len(set2)} 个唯一 k-mer")

    only_in_1 = set1 - set2
    only_in_2 = set2 - set1
    common = set1 & set2

    lines = [
        f"k-mer 长度: {args.k}",
        f"文件1 ({args.fasta1}): {len(set1)} 个唯一 k-mer",
        f"文件2 ({args.fasta2}): {len(set2)} 个唯一 k-mer",
        f"仅在文件1中: {len(only_in_1)}",
        f"仅在文件2中: {len(only_in_2)}",
        f"共有: {len(common)}",
    ]

    if args.output:
        with open(args.output, 'w') as f:
            f.write('\n'.join(lines) + '\n')
            f.write("\n# 仅在文件1中的 k-mer:\n")
            f.write('\n'.join(sorted(only_in_1)) + '\n')
            f.write("\n# 仅在文件2中的 k-mer:\n")
            f.write('\n'.join(sorted(only_in_2)) + '\n')
        print(f"结果已保存至 {args.output}")
    else:
        print('\n'.join(lines))
        print("\n# 仅在文件1中的 k-mer (前 10 个):")
        for kmer in sorted(only_in_1)[:10]:
            print(kmer)
        print("\n# 仅在文件2中的 k-mer (前 10 个):")
        for kmer in sorted(only_in_2)[:10]:
            print(kmer)

if __name__ == "__main__":
    main()