#!/usr/bin/env python3
"""
extract_vcf_from_mpileup.py

从 bcftools mpileup 生成的 VCF/BCF 中提取候选变异位点。
规则：只要某个位点存在至少一个非参考等位基因（且其深度 ≥ 1），就输出。
输出为简单的 VCF 文件，供 hap.py 或其他评估工具使用。

用法:
  python extract_vcf_from_mpileup.py \
      --pileup all_sites.bcf \
      --output candidates.vcf \
      --sample-name SAMPLE
"""

import argparse
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="从 mpileup 输出中提取候选变异"
    )
    parser.add_argument("--pileup", required=True,
                        help="bcftools mpileup 输出文件（.bcf 或 .vcf）")
    parser.add_argument("--output", required=True,
                        help="输出 VCF 文件路径")
    parser.add_argument("--sample-name", default="SAMPLE",
                        help="样本名（写入 VCF 头部，默认 SAMPLE）")
    args = parser.parse_args()

    pileup_path = Path(args.pileup)
    if not pileup_path.is_file():
        sys.exit(f"错误: 找不到输入文件 '{pileup_path}'")

    # 1. 使用 bcftools query 提取位点信息：CHROM POS REF ALT AD
    query_cmd = [
        "bcftools", "query",
        "-f", "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n",
        str(pileup_path)
    ]
    result = subprocess.run(query_cmd, capture_output=True, text=True, check=True)

    lines = result.stdout.strip().splitlines()
    if not lines:
        print("警告: pileup 文件中没有位点。")
        return

    # 2. 解析并筛选候选变异
    candidates = []
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        chrom, pos, ref, alt = fields[0], fields[1], fields[2], fields[3]
        ad = fields[4] if len(fields) >= 5 else ""

        # 跳过参考位点（ALT 为空或 '.'）
        if alt in ("", "."):
            continue

        # 解析 AD（等位基因深度），格式如 "2,1,0" 对应 REF,ALT1,ALT2...
        has_alt = False
        alt_depth = 1  # 默认值，用于输出
        if ad:
            try:
                ad_parts = [int(x) for x in ad.split(",")]
                if len(ad_parts) > 1:
                    # 检查是否有任何非参考等位基因深度 >= 1
                    non_ref_depths = ad_parts[1:]
                    if max(non_ref_depths) >= 1:
                        has_alt = True
                        alt_depth = max(non_ref_depths)
                else:
                    # AD 只有一个值，但 ALT 非空，保守认为存在
                    has_alt = True
            except ValueError:
                # AD 解析失败，保守认为存在
                has_alt = True
        else:
            # 没有 AD 字段，但 ALT 非空，保留
            has_alt = True

        if has_alt:
            # 只取第一个 ALT 等位基因（简化）
            alt_allele = alt.split(",")[0]
            candidates.append((chrom, pos, ref, alt_allele, alt_depth))

    # 3. 写 VCF 文件
    vcf_lines = [
        "##fileformat=VCFv4.2",
        "##source=extract_candidates_from_mpileup.py",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample_name}"
    ]

    for chrom, pos, ref, alt, alt_depth in candidates:
        vcf_lines.append(
            f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t30\tPASS\t.\tGT:AD\t0/1:{alt_depth},{alt_depth}"
        )

    with open(args.output, "w") as f:
        f.write("\n".join(vcf_lines) + "\n")

    print(f"完成: 共提取 {len(candidates)} 个候选变异位点，已写入 {args.output}")


if __name__ == "__main__":
    main()