#!/usr/bin/env python3
import sys
import csv
from cyvcf2 import VCF
from pyfaidx import Fasta
import argparse
from collections import defaultdict

def range_in_set(chrom, start, end, s):
    for i in range(start,end):
        if i in s[chrom]:
            return True
    return False

def extract_flanking_kmers(ref_path, vcf_path, output_path, k=31):
    """
    提取变异位点紧挨着的左右 K-mer 并保存到 CSV 文件。
    """
    # 1. 初始化参考基因组索引 (自动处理 1-based 坐标)
    print("正在加载参考基因组索引...")
    ref = Fasta(ref_path)
    #print(ref['chr22'])
    # 2. 初始化 VCF 读取器
    print("正在读取 VCF 文件...")
    vcf = VCF(vcf_path)
    varpos = defaultdict(set)
    for variant in vcf:
        for b in range(variant.POS-1 , variant.POS-1+len(variant.REF)):
            varpos[variant.CHROM].add(b)

    vcf = VCF(vcf_path)

    # 3. 准备写入 CSV 文件
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # 写入表头
        writer.writerow(["Chrom", "Variant_Pos", "Ref_Base", "Side", "Kmer"])

        count = 0
        for variant in vcf:
            chrom = variant.CHROM
            pos = variant.POS - 1     # 1-based
            ref_base = variant.REF
            alt_base = variant.ALT[0]

            is_snp = len(ref_base) == len(alt_base)

            # 检查染色体是否在参考基因组中
            if chrom not in ref:
                print(f"警告: 染色体 {chrom} 未在参考基因组中找到，跳过。")
                continue

            chrom_len = len(ref[chrom])
            old_count = count
            # --- 提取左侧 K-mer ---
            # 左侧 K-mer 的 1-based 范围: [pos-k, pos-1]
            left_end = pos  if is_snp else pos + 1
            left_start = left_end - k
            if left_start >= 1 and not range_in_set(chrom, left_start, left_end, varpos):  # 确保不越界
                left_kmer = ref[chrom][left_start:left_end].seq.upper()
                writer.writerow([chrom, pos+1, ref_base, "Left", left_kmer])
                count += 1

            # --- 提取右侧 K-mer ---
            # 右侧 K-mer 的 1-based 范围: [pos+1, pos+k]
            right_start = pos + 1 if is_snp else pos + len(ref_base)
            right_end = right_start + k
            if right_end <= chrom_len and not range_in_set(chrom, right_start, right_end, varpos):  # 确保不越界
                right_kmer = ref[chrom][right_start:right_end].seq.upper()
                writer.writerow([chrom, pos+1, ref_base, "Right", right_kmer])
                count += 1

            #if old_count == count: print(chrom , pos+1, is_snp, ref_base, alt_base)
    vcf.close()
    print(f"提取完成！共生成 {count} 条 K-mer 记录，已保存至: {output_path}")

def main(argv):
    parser = argparse.ArgumentParser(description='Extract flanking kmers from VCF file.')
    parser.add_argument('ref', help='Path to reference genome FASTA file.')
    parser.add_argument('vcf', help='Path to VCF file.')
    parser.add_argument('out', help='Path to output CSV file.')
    parser.add_argument('-k', type=int, default=16, help='K-mer size (default: 31).')

    #REF_FILE = "../../chr22.fa"
    #VCF_FILE = "../../../../data/raw/人类chr22-1/HG001_chr22_het.vcf.gz"
    #OUT_FILE = "flanking_kmers.csv"
    #K_SIZE = 16  # 设定你的 K-mer 长度

    try:
        args = parser.parse_args()
        extract_flanking_kmers(args.ref, args.vcf, args.out, k=args.k)
    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)
    
if __name__ == "__main__":
    main(sys.argv)
