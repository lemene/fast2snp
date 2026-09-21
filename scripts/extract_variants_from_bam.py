#!/usr/bin/env python3
"""
extract_variants_from_bam.py

直接从 BAM 文件中提取所有与参考基因组不同的位点（SNP 和 indel），
输出符合 hap.py 格式要求的 VCF。
每个替代等位基因单独一行，包含 FORMAT/SAMPLE 和 VT 标签。
跳过参考等位基因为 N 的位点。

用法:
  python extract_variants_from_bam.py \
      --bam input.bam \
      --ref reference.fasta \
      --output variants.vcf \
      [--min-mapq 0] \
      [--min-baseq 0]
"""

import argparse
import sys
from collections import defaultdict
import pysam

def main():
    parser = argparse.ArgumentParser(description="从 BAM 提取所有 SNP 和 indel，输出标准 VCF")
    parser.add_argument("--bam", required=True, help="输入 BAM 文件（需已排序和索引）")
    parser.add_argument("--ref", required=True, help="参考基因组 FASTA（需已建立索引）")
    parser.add_argument("--output", required=True, help="输出 VCF 文件")
    parser.add_argument("--min-mapq", type=int, default=0, help="最低比对质量（默认 0 不过滤）")
    parser.add_argument("--min-baseq", type=int, default=0, help="最低碱基质量（默认 0 不过滤）")
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    ref_fasta = pysam.FastaFile(args.ref)

    # 从 BAM 头部提取所有参考序列名称（用于 contig 行）
    contigs = bam.references

    # 存储变异信息: {(chrom, pos): {"ref": ref_base, "alts": set(), "types": {alt: vt}}}
    variant_dict = defaultdict(lambda: {"ref": "", "alts": set(), "types": {}})

    def add_variant(chrom, pos1, ref_allele, alt_allele, vt):
        # 跳过参考为 N 或含 N 的情况
        if ref_allele == "N" or "N" in ref_allele:
            return
        key = (chrom, pos1)
        variant_dict[key]["ref"] = ref_allele
        variant_dict[key]["alts"].add(alt_allele)
        variant_dict[key]["types"][alt_allele] = vt

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        if read.mapping_quality < args.min_mapq:
            continue

        if read.query_sequence is None:
            continue
        if read.query_qualities is None:
            # 若质量缺失，用空列表代替
            read.query_qualities = []

        chrom = read.reference_name
        ref_pos = read.reference_start  # 0-based
        query_pos = 0
        query_seq = read.query_sequence
        query_qual = read.query_qualities if read.query_qualities else []

        for op, length in read.cigartuples:
            if op == 0:  # M (match or mismatch)
                for i in range(length):
                    ref_base = ref_fasta.fetch(chrom, ref_pos + i, ref_pos + i + 1).upper()
                    # 跳过参考为 N 的位点
                    if not ref_base or ref_base == "N":
                        continue
                    q_base = query_seq[query_pos + i].upper()
                    q_qual = query_qual[query_pos + i] if query_qual else 0
                    if q_qual < args.min_baseq:
                        continue
                    if ref_base and ref_base != q_base:
                        # SNP
                        pos1 = ref_pos + i + 1  # 1-based
                        add_variant(chrom, pos1, ref_base, q_base, "SNP")
                ref_pos += length
                query_pos += length

            elif op == 1:  # I (insertion)
                ins_seq = query_seq[query_pos:query_pos+length]
                anchor_pos = ref_pos  # 0-based, insertion after anchor
                anchor_base = ref_fasta.fetch(chrom, anchor_pos - 1, anchor_pos) if anchor_pos > 0 else ref_fasta.fetch(chrom, anchor_pos, anchor_pos + 1)
                # 跳过锚点碱基为 N 的插入
                if not anchor_base or anchor_base == "N":
                    query_pos += length
                    continue
                ref_allele = anchor_base
                alt_allele = anchor_base + ins_seq
                pos1 = anchor_pos
                add_variant(chrom, pos1, ref_allele, alt_allele, "INDEL")
                query_pos += length

            elif op == 2:  # D (deletion)
                del_seq = ref_fasta.fetch(chrom, ref_pos, ref_pos + length)
                anchor_pos = ref_pos  # 0-based, deletion after anchor
                if anchor_pos > 0:
                    anchor_base = ref_fasta.fetch(chrom, anchor_pos - 1, anchor_pos)
                else:
                    anchor_base = ref_fasta.fetch(chrom, anchor_pos + length, anchor_pos + length + 1)
                # 跳过锚点或缺失序列含 N 的情况
                if not anchor_base or anchor_base == "N" or not del_seq or "N" in del_seq:
                    ref_pos += length
                    continue
                ref_allele = anchor_base + del_seq
                alt_allele = anchor_base
                pos1 = anchor_pos
                add_variant(chrom, pos1, ref_allele, alt_allele, "INDEL")
                ref_pos += length

            elif op == 4:  # S (soft clip)
                query_pos += length
            elif op == 5:  # H (hard clip)
                pass
            else:
                if op == 3:  # N
                    ref_pos += length
                elif op == 6:  # P
                    pass

    bam.close()
    ref_fasta.close()

    # 输出 VCF
    with open(args.output, "w") as vcf:
        # 头部
        vcf.write("##fileformat=VCFv4.1\n")
        for contig in contigs:
            vcf.write(f"##contig=<ID={contig}>\n")
        vcf.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        vcf.write('##INFO=<ID=VT,Number=1,Type=String,Description="Variant type">\n')
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        # 变异记录
        for (chrom, pos), info in sorted(variant_dict.items()):
            ref = info["ref"]
            alts = info["alts"]
            if not alts or not ref:
                continue
            # 输出前再次确保 REF 和 ALT 不含 N
            if "N" in ref:
                continue
            for alt in sorted(alts):
                if "N" in alt:
                    continue
                vt = info["types"].get(alt, "SNP" if len(ref) == len(alt) == 1 else "INDEL")
                vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\tVT={vt}\tGT\t0|1\n")

    print(f"完成，共输出 {len(variant_dict)} 个变异位点")

if __name__ == "__main__":
    main()