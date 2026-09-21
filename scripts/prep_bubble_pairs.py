#!/usr/bin/env python3
"""
prep_bubble_pairs.py

从 fast2snp 输出的 FASTA 文件中提取 bubble 序列名称，
按照 Bubble_<ID>_<TYPE>_Path<1|2>_Ext<X> 模式配对，
将每对序列名输出到指定文件（每行两个名称，空格分隔）。

用法:
  python prep_bubble_pairs.py \
      --fasta input.fasta \
      --output bubble_pairs.txt \
      --pattern "Bubble_{id}_{type}_Path{path}_Ext{ext}"
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict


def parse_bubble_name(seq_name, pattern="Bubble_{id}_{type}_Path{path}_Ext{ext}"):
    """
    从序列名称中解析 bubble ID、类型（SNP/INDEL）、路径编号（1/2）和 Ext 值。
    
    默认模式: Bubble_{id}_{type}_Path{path}_Ext{ext}
    例如: Bubble_52482_INDEL_Path1_Ext8
          -> id=52482, type=INDEL, path=1, ext=8
    
    返回 dict 或 None（如果不匹配）
    """
    # 构建正则表达式
    # 将模式转换为带命名组的正则
    regex_pattern = pattern
    regex_pattern = regex_pattern.replace("{id}", "(?P<id>[^_]+)")
    regex_pattern = regex_pattern.replace("{type}", "(?P<type>SNP|INDEL|Snp|Indel|IsolatedSNP|NonIsolatedSNP)")
    regex_pattern = regex_pattern.replace("{path}", "(?P<path>[12])")
    regex_pattern = regex_pattern.replace("{ext}", "(?P<ext>[^_]+)")
    
    match = re.match(regex_pattern, seq_name)
    if match:
        return {
            'id': match.group('id'),
            'type': 'SNP' if match.group('type').upper().find("SNP") >= 0 else "INDEL",
            'path': int(match.group('path')),
            'ext': match.group('ext')
        }
    return None


def extract_bubble_pairs(fasta_file, pattern="Bubble_{id}_{type}_Path{path}_Ext{ext}"):
    """
    从 FASTA 文件中提取所有序列名，并按 bubble 配对。
    
    返回:
      pairs: list of (name1, name2, type)
      unpaired: list of 未能配对的序列名
    """
    seq_names = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seq_names.append(line.split()[0][1:])
    
    bubbles = defaultdict(dict)  # {(id, type): {path: seq_name}}
    unpaired = []
    
    for name in seq_names:
        info = parse_bubble_name(name, pattern)
        if info is None:
            unpaired.append(name)
            continue
        
        key = (info['id'], info['type'])
        path = info['path']
        
        if path in bubbles[key]:
            print(f"警告: 重复的路径 {name}")
            unpaired.append(name)
        else:
            bubbles[key][path] = name
    
    # 配对
    pairs = []
    for (bid, btype), paths in sorted(bubbles.items(), key=lambda x: int(x[0][0]) if x[0][0].isdigit() else x[0][0]):
        if 1 in paths and 2 in paths:
            pairs.append((paths[1], paths[2], btype))
        else:
            missing = [p for p in [1, 2] if p not in paths]
            present = [paths[p] for p in [1, 2] if p in paths]
            print(f"警告: Bubble {bid} ({btype}) 缺少 Path{missing}，现有: {present}")
            for name in present:
                unpaired.append(name)
    
    return pairs, unpaired


def main():
    parser = argparse.ArgumentParser(
        description="从 FASTA 提取 bubble 序列名并配对"
    )
    parser.add_argument("--fasta", required=True,
                        help="输入 FASTA 文件（fast2snp 输出）")
    parser.add_argument("--output", required=True,
                        help="输出配对文件（每行两个序列名和类型，空格分隔）")
    parser.add_argument("--pattern", default="Bubble_{id}_{type}_Path{path}_Ext{ext}",
                        help="序列命名模式（默认: Bubble_{id}_{type}_Path{path}_Ext{ext}）")
    parser.add_argument("--unpaired-output", default=None,
                        help="输出未配对序列名到指定文件（可选）")
    args = parser.parse_args()
    
    if not Path(args.fasta).is_file():
        sys.exit(f"错误: 找不到文件 '{args.fasta}'")
    
    print("解析 FASTA 文件...")
    pairs, unpaired = extract_bubble_pairs(args.fasta, args.pattern)
    
    print(f"  共找到 {len(pairs)} 对 bubble")
    print(f"  未配对序列: {len(unpaired)}")
    
    # 输出配对（三列：path1_name path2_name type）
    with open(args.output, 'w') as f:
        for name1, name2, btype in pairs:
            f.write(f"{name1} {name2} {btype}\n")
    
    print(f"配对结果已保存至: {args.output}")
    
    if args.unpaired_output:
        with open(args.unpaired_output, 'w') as f:
            for name in unpaired:
                f.write(f"{name}\n")
        print(f"未配对序列已保存至: {args.unpaired_output}")
    elif unpaired:
        print("\n警告: 存在未配对序列，前几个如下:")
        for name in unpaired[:10]:
            print(f"  {name}")
        if len(unpaired) > 10:
            print(f"  ... 共 {len(unpaired)} 个")


if __name__ == "__main__":
    main()