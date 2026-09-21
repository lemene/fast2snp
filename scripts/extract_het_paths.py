#!/usr/bin/env python3
"""
extract_variant_paths.py

从 VCF 文件中提取变异的两条路径（参考和替代），支持自动合并临近变异。
路径包含变异位点上下游指定长度的侧翼序列，保证除变异区域外完全相同。
对于 INDEL，自动去除锚定碱基（首碱基相同）以保证侧翼一致。

合并规则：
  默认将距离小于 flank 的变异合并到同一条路径中，避免侧翼重叠。
  可通过 --merge-distance 覆盖合并距离（设 0 表示不合并）。

用法:
  python extract_variant_paths.py \
      --ref reference.fasta \
      --vcf variants.vcf.gz \
      --flank 100 \
      --output variant_paths.fasta

可选:
  --pass-only        只处理 FILTER 为 PASS 或 '.' 的变异
  --merge-distance N 合并距离阈值，默认等于 flank（设 0 表示不合并）
"""

import argparse
import pysam

def normalize_variant(record, alt):
    """
    处理单个变异，返回 (actual_pos, actual_ref, actual_alt, var_type)
    对于 INDEL，若 REF 和 ALT 首碱基相同，去除锚定碱基，位置后移一位。
    """
    chrom = record.chrom
    pos = record.pos
    ref_allele = record.ref
    is_snp = (len(ref_allele) == 1 and len(alt) == 1)
    var_type = "SNP" if is_snp else "INDEL"

    if var_type == "INDEL" and ref_allele[0] == alt[0]:
        actual_ref = ref_allele[1:]
        actual_alt = alt[1:]
        actual_pos = pos + 1
    else:
        actual_ref = ref_allele
        actual_alt = alt
        actual_pos = pos
    return actual_pos, actual_ref, actual_alt, var_type

def main():
    parser = argparse.ArgumentParser(description="提取变异路径序列，支持合并临近变异")
    parser.add_argument("--ref", required=True, help="参考基因组 FASTA")
    parser.add_argument("--vcf", required=True, help="输入 VCF 文件")
    parser.add_argument("--flank", type=int, default=100, help="侧翼长度（默认100）")
    parser.add_argument("--output", required=True, help="输出 FASTA 文件")
    parser.add_argument("--pass-only", action="store_true", help="只处理 PASS 或 '.' 记录")
    parser.add_argument("--merge-distance", type=int, default=None,
                        help="合并距离阈值，默认等于 flank，设 0 表示不合并")
    args = parser.parse_args()

    # 合并距离：默认等于 flank
    if args.merge_distance is None:
        merge_dist = args.flank
    else:
        merge_dist = args.merge_distance

    ref_fasta = pysam.FastaFile(args.ref)
    vcf = pysam.VariantFile(args.vcf)

    out_fasta = open(args.output, "w")
    written = 0
    skipped = 0

    variants = []
    for record in vcf:
        if args.pass_only:
            if record.filter.keys() and "PASS" not in record.filter.keys() and "." not in record.filter.keys():
                continue
        chrom = record.chrom
        for alt in record.alts:
            if alt is None or alt in ("<NON_REF>", "<*>"):
                continue
            actual_pos, actual_ref, actual_alt, var_type = normalize_variant(record, alt)
            if not actual_ref:
                continue
            variants.append((chrom, actual_pos, actual_ref, actual_alt, var_type, record.pos, alt))

    # 按位置排序
    variants.sort(key=lambda x: (x[0], x[1]))

    # 分组：距离小于 merge_dist 的变异合并
    groups = []
    if variants:
        current_group = [variants[0]]
        last_end = variants[0][1] + len(variants[0][2])  # 实际变异末端（1-based）
        for var in variants[1:]:
            if var[1] - last_end < merge_dist:
                current_group.append(var)
                last_end = max(last_end, var[1] + len(var[2]))
            else:
                groups.append(current_group)
                current_group = [var]
                last_end = var[1] + len(var[2])
        groups.append(current_group)

    # 对每个组生成路径
    for group in groups:
        chrom = group[0][0]
        min_pos = min(v[1] for v in group)
        max_pos = max(v[1] + len(v[2]) - 1 for v in group)
        start = min_pos - args.flank
        end = max_pos + args.flank

        chrom_len = ref_fasta.get_reference_length(chrom)
        if start < 1 or end > chrom_len:
            skipped += len(group)
            continue

        ref_window = ref_fasta.fetch(chrom, start - 1, end).upper()

        # 校验所有变异在窗口内的 REF 一致
        valid = True
        for var in group:
            pos_in_window = var[1] - start
            if ref_window[pos_in_window:pos_in_window+len(var[2])] != var[2].upper():
                valid = False
                break
        if not valid:
            skipped += len(group)
            continue

        # 生成替代序列：在参考窗口上依次替换每个变异
        alt_seq = list(ref_window)
        for var in group:
            pos_in_window = var[1] - start
            alt_seq[pos_in_window:pos_in_window+len(var[2])] = var[3].upper()
        alt_seq = ''.join(alt_seq)

        # 构建 header 中的位置和类型信息
        if len(group) == 1:
            # 单个变异：保留原格式
            pos_str = str(group[0][5])
            type_str = group[0][4]
            out_fasta.write(f">{chrom}:{pos_str}_TYPE_{type_str}_REF\n{ref_window}\n")
            out_fasta.write(f">{chrom}:{pos_str}_TYPE_{type_str}_ALT\n{alt_seq}\n")
        else:
            # 合并多个变异：位置列表 + 类型列表，无范围
            pos_list = "_".join(str(v[5]) for v in group)
            type_list = "_".join(v[4] for v in group)
            out_fasta.write(f">{chrom}:{pos_list}_VARIANTS_{type_list}_REF\n{ref_window}\n")
            out_fasta.write(f">{chrom}:{pos_list}_VARIANTS_{type_list}_ALT\n{alt_seq}\n")
        written += 2

    out_fasta.close()
    vcf.close()
    ref_fasta.close()

    print(f"完成：输出了 {written} 条路径序列，跳过 {skipped} 个变异。")

if __name__ == "__main__":
    main()