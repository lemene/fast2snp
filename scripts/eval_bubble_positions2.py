#!/usr/bin/env python3
"""
eval_bubble_positions.py

通用 bubble 位置评估框架，支持多种 reference-free SNP 检测工具。
完整流程：转换 bubble → 提取序列 → BWA 比对 → 构建位置索引 → 与真值比较。

用法:
  python eval_bubble_positions.py \
      --tool fast2snp \
      --bubble-file calls.fasta \
      --ref reference.fa \
      --truth truth.vcf.gz \
      --bed confident.bed \
      --output-prefix eval_result \
      --mapq 20 \
      --threads 8
"""

import argparse
import sys
import os
import re
import logging
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
from abc import ABC, abstractmethod
import pysam
import shutil

# ============================================================
# 日志配置
# ============================================================
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


# ============================================================
# 基类：BubbleEvaluator
# ============================================================
class BubbleEvaluator(ABC):
    """所有工具评估器的基类"""
    
    def __init__(self, bubble_file, ref_fasta, truth_vcf, bed_file,
                 min_mapq=20, threads=8, work_dir=None, keep_bam=False):
        self.bubble_file = bubble_file
        self.ref_fasta = ref_fasta
        self.truth_vcf = truth_vcf
        self.bed_file = bed_file
        self.min_mapq = min_mapq
        self.threads = threads
        self.pairs = []  # [(name1, name2, type)]
        
        # 工作目录
        if work_dir:
            self.work_dir = Path(work_dir)
            self.work_dir.mkdir(parents=True, exist_ok=True)
            self._temp_dir = None
        else:
            self._temp_dir = tempfile.TemporaryDirectory(prefix="bubble_eval_")
            self.work_dir = Path(self._temp_dir.name)
        
        self.keep_bam = keep_bam
        self.bam_path = None
        
        # 生成的文件路径
        self.fasta_path = self.work_dir / "bubble_sequences.fasta"
        self.pairs_path = self.work_dir / "bubble_pairs.txt"
    
    @abstractmethod
    def _convert_to_fasta(self, output_fasta, output_pairs):
        """
        将 bubble 文件转换为 FASTA 文件和 pairs 文件。
        
        参数:
          output_fasta: 输出 FASTA 文件路径（用于 BWA 比对）
          output_pairs: 输出 pairs 文件路径（每行: path1_name path2_name type）
        
        子类必须实现此方法。
        """
        pass
    
    def load_bubbles(self):
        """从生成的 pairs 文件中加载 bubbles"""
        pairs = []

        with open(self.pairs_path, 'r') as f:
            for line in f:
                parts = line.split()
                assert len(parts) >= 3, f"行格式错误: {line}"

                pairs.append((parts[0:-1], parts[-1]))
        return pairs
    
    def align_bubbles_to_reference(self):
        """将 bubble 序列比对到参考基因组，返回排序后 BAM 路径"""
        logger.info("检查 BWA 索引...")
        if not Path(self.ref_fasta + ".bwt").exists():
            logger.info("BWA 索引不存在，开始构建...")
            subprocess.run(
                ["bwa", "index", self.ref_fasta],
                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            logger.info("BWA 索引构建完成")
        else:
            logger.info("BWA 索引已存在")
        
        # BWA MEM 比对
        logger.info("执行 BWA MEM 比对...")
        sam_path = self.work_dir / "bubble_alignments.sam"
        bwa_cmd = [
            "bwa", "mem",
            "-t", str(self.threads),
            "-k", "13",
            "-T", "10",
            "-a",
            self.ref_fasta,
            str(self.fasta_path)
        ]
        
        with open(sam_path, 'w') as sam_out:
            subprocess.run(bwa_cmd, stdout=sam_out, check=True)
        
        # SAM → BAM → 排序 → 索引
        logger.info("转换并排序 BAM...")
        unsorted_bam = self.work_dir / "bubble_alignments.unsorted.bam"
        sorted_bam = self.work_dir / "bubble_alignments.sorted.bam"
        
        subprocess.run(
            ["samtools", "view", "-bS", str(sam_path), "-o", str(unsorted_bam)],
            check=True
        )
        subprocess.run(
            ["samtools", "sort", str(unsorted_bam), "-o", str(sorted_bam)],
            check=True
        )
        subprocess.run(
            ["samtools", "index", str(sorted_bam)],
            check=True
        )
        
        self.bam_path = str(sorted_bam)
        logger.info(f"排序后的 BAM: {self.bam_path}")
        
        if not self.keep_bam:
            os.remove(sam_path)
            os.remove(unsorted_bam)
        
        return self.bam_path
    
    def load_bed_regions(self):
        """加载 BED 区域"""
        regions = defaultdict(list)
        with open(self.bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split("\t")
                if len(parts) >= 3:
                    regions[parts[0]].append((int(parts[1]), int(parts[2])))
        return regions
    
    def is_in_bed(self, chrom, pos1, bed_regions):
        """检查 1-based VCF 位置是否在 BED 内"""
        if chrom not in bed_regions:
            return False
        pos0 = pos1 - 1
        for start, end in bed_regions[chrom]:
            if start <= pos0 < end:
                return True
        return False
    
    def load_truth_variants(self, bed_regions):
        """加载真值变异，按类型分类"""
        snp_positions = set()
        indel_positions = set()
        vcf = pysam.VariantFile(self.truth_vcf)
        for record in vcf:
            if not self.is_in_bed(record.chrom, record.pos, bed_regions):
                continue
            for alt in record.alts:
                if alt is None or alt == '<NON_REF>' or alt == '<*>':
                    continue
                if len(record.ref) == 1 and len(alt) == 1:
                    snp_positions.add((record.chrom, record.pos))
                else:
                    indel_positions.add((record.chrom, record.pos))
        vcf.close()
        return {"SNP": snp_positions, "INDEL": indel_positions}
    
    def build_position_index(self):
        """构建位置索引：{(chrom, pos0): set(bubble_idx)}，按类型分开"""
        pos_index_by_type = {
            "SNP": defaultdict(set),
            "INDEL": defaultdict(set)
        }
        
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        query_to_alignments = defaultdict(list)
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            if read.mapping_quality < self.min_mapq:
                continue
            query_to_alignments[read.query_name].append(
                (read.reference_name, read.reference_start, read.reference_end)
            )
        bam.close()

        for bubble_idx, (names, btype) in enumerate(self.pairs):
            
            all_alignments = []
            for name in names:
                if name in query_to_alignments:
                    all_alignments.extend(query_to_alignments[name])

            target_index = pos_index_by_type[btype]
            for ref_name, ref_start, ref_end in all_alignments:
                for pos0 in range(ref_start, ref_end):
                    target_index[(ref_name, pos0)].add(bubble_idx)
        
        return pos_index_by_type

    def write_pairs(self, pairs_path, pairs):
        """pairs: [[], ""]"""

        type_counts = defaultdict(int)
        with open(pairs_path, 'w') as f:
            for names, btype in pairs.values():
                if len(names) > 1 and btype:
                    type_counts[btype] += 1
                    f.write(f"{' '.join(names)} {btype}\n")
                else:
                    logger.warning(f"Bubble 配对不完整: size={len(names)}, paths='{names}', type='{btype}'")
        
        logger.info(f"  共转换 {len(pairs)} 对 bubble")
        logger.info(f"    SNP: {type_counts.get('SNP', 0)}")
        logger.info(f"    INDEL: {type_counts.get('INDEL', 0)}")

    def evaluate(self):
        """执行完整评估流程"""
        logger.info(f"转换 bubble 文件（工具: {self.tool_name}）...")
        self._convert_to_fasta(self.fasta_path, self.pairs_path)
        logger.info(f"  FASTA: {self.fasta_path}")
        logger.info(f"  Pairs: {self.pairs_path}")
        
        logger.info("加载 bubble pairs...")
        self.pairs = self.load_bubbles()
        logger.info(f"  共 {len(self.pairs)} 个 bubbles")
        
        type_counts = defaultdict(int)
        pairs_by_type = {"SNP": [], "INDEL": []}
        for i, (_, bt) in enumerate(self.pairs):
            type_counts[bt] += 1
            pairs_by_type[bt].append(i)
        
        logger.info(f"  SNP bubbles: {type_counts.get('SNP', 0)}")
        logger.info(f"  INDEL bubbles: {type_counts.get('INDEL', 0)}")
        
        self.align_bubbles_to_reference()
        
        logger.info("加载 BED 区域...")
        bed_regions = self.load_bed_regions()
        
        logger.info("加载真值变异...")
        truth_by_type = self.load_truth_variants(bed_regions)
        logger.info(f"  真值 SNP: {len(truth_by_type['SNP'])}")
        logger.info(f"  真值 INDEL: {len(truth_by_type['INDEL'])}")
        
        logger.info("构建位置索引...")
        pos_index_by_type = self.build_position_index()
        for vtype in ["SNP", "INDEL"]:
            logger.info(f"  {vtype} 索引包含 {len(pos_index_by_type[vtype])} 个唯一位置")
        
        results = {}
        for vtype in ["SNP", "INDEL"]:
            type_indices = set(pairs_by_type[vtype])
            pos_index = pos_index_by_type[vtype]
            truth_positions = truth_by_type[vtype]
            
            tp_bubbles = set()
            covered_truth = set()
            
            for (chrom, pos1) in truth_positions:
                pos0 = pos1 - 1
                bubbles_here = pos_index.get((chrom, pos0), set())
                if bubbles_here:
                    covered_truth.add((chrom, pos1))
                    tp_bubbles.update(bubbles_here)
            
            fp_bubbles = type_indices - tp_bubbles
            fn_positions = truth_positions - covered_truth
            
            results[vtype] = {
                'tp': len(tp_bubbles),
                'fp': len(fp_bubbles),
                'fn': len(fn_positions),
                'tp_bubbles': tp_bubbles,
                'fp_bubbles': fp_bubbles,
                'fn_positions': fn_positions,
                'covered_truth': covered_truth,
                'total_bubbles': len(type_indices),
                'total_truth': len(truth_positions)
            }
        
        return results
    
    def cleanup(self):
        """清理临时文件"""
        if self._temp_dir and not self.keep_bam:
            self._temp_dir.cleanup()


# ============================================================
# Fast2snp 评估器
# ============================================================
class Fast2snpEvaluator(BubbleEvaluator):
    """
    fast2snp 评估器
    输入: FASTA，header 格式 Bubble_<id>_<type>_Path<1|2>_Ext<X>
    """
    
    tool_name = "fast2snp"
    
    def _convert_to_fasta(self, output_fasta, output_pairs):
        """
        从 fast2snp FASTA 文件生成:
          - FASTA 文件（保持原样）
          - pairs 文件（每行: path1_name path2_name type）
        """
        shutil.copy(self.bubble_file, output_fasta)
        
        pairs = defaultdict(lambda: [[], ""])

        with open(self.bubble_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    name = line.split()[0][1:]
                    match = re.match(r'Bubble_(\d+)_(INDEL|IsolatedSNP|NonIsolatedSNP)_Path([12])_Ext(\S+)', name)
                    if match:
                        bid = match.group(1)
                        btype = match.group(2).upper()
                        path = int(match.group(3))
                        assert path in (1, 2)
                        pairs[bid][0].append(name)
                        pairs[bid][1] = 'SNP' if btype.upper().find("SNP") >= 0 else "INDEL";

                    else:
                        logger.warning(f"无法解析 header: {name}")
        
        self.write_pairs(output_pairs, pairs)

# ============================================================
# Kmer2snp 评估器
# ============================================================
class Kmer2snpEvaluator(BubbleEvaluator):
    """
    kmer2snp 评估器
    输入: 表格格式 Seq1 Seq2 Count
    """
    
    tool_name = "kmer2snp"
    
    def _convert_to_fasta(self, output_fasta, output_pairs):

        pairs = defaultdict(lambda: [[], ""])

        """将表格格式转换为 FASTA + pairs 文件"""
        with open(self.bubble_file, 'r') as inf, \
             open(output_fasta, 'w') as fasta_out:
            
            for i, line in enumerate(inf):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    logger.warning(f"第{i}行格式错误，跳过")
                    continue
                
                seq1 = parts[0]
                seq2 = parts[1]
                
                name1 = f"kmer2snp_{i}_path1"
                name2 = f"kmer2snp_{i}_path2"

                fasta_out.write(f">{name1}\n{seq1}\n")
                fasta_out.write(f">{name2}\n{seq2}\n")

                pairs[i][0] = [name1, name2]
                pairs[i][1] = "SNP"
                
        self.write_pairs(output_pairs, pairs)

# ============================================================
# DiscoSnp 评估器
# ============================================================
class DiscoSnpEvaluator(BubbleEvaluator):
    """
    discoSnp++ 评估器
    输入: FASTA，header 格式:
      >SNP_<id>_higher_path_<info>
      >SNP_<id>_lower_path_<info>
      >INDEL_<id>_higher_path_<info>
      >INDEL_<id>_lower_path_<info>
    其中 <info> 可能包含 |P_... 等附加信息
    """
    
    tool_name = "discoSnp"
    
    def _convert_to_fasta(self, output_fasta, output_pairs):
        """
        从 discoSnp++ FASTA 文件生成:
          - FASTA 文件（保持原样）
          - pairs 文件（每行: higher_path_name lower_path_name type）
        """

        # 复制原始 FASTA
        shutil.copy(self.bubble_file, output_fasta)
        
        # 解析 header 并配对
        pairs = defaultdict(lambda: [[], ""])
        
        with open(self.bubble_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    name = line.split()[0][1:]
                    segs = name.split('|')[0].split('_')

                    if len(segs) == 4 and (segs[0] == "SNP" or segs[0] == "INDEL") and segs[2] == "path":
                        id = int(segs[3])
                        path = 1 if segs[1] == "higher" else 2
                        pairs[id][0].append(name)
                        pairs[id][-1] = segs[0]
                    else:
                        logger.warning(f"无法解析 header: {name}")
                        continue
                        
        self.write_pairs(output_pairs, pairs)

# ============================================================
# Ebwt2Indel 评估器
# ============================================================
class Ebwt2IndelEvaluator(BubbleEvaluator):
    """
    ebwt2InDel 评估器
    输入: FASTA，header 格式:
      >cluster:<id>_id:<d>_right:<len>_cov:<cov>_type:_<SNP|INDEL>_event:<REF/ALT>
    
    示例:
      >cluster:2252_id:1_right:30_cov:4_type:_INDEL_event:A/
      GAAAAACTTGGTTTTAATTAAAAACAAACCAAAAAAAAAAAAAAAAAAACACCTCGGGGAC
      >cluster:2252_id:2_right:30_cov:3_type:_INDEL_event:A/
      GGAAAAACTTGGTTTTAATTAAAAACAAACCAAAAAAAAAAAAAAAAAACACCTCGGGGAC
    """
    
    tool_name = "ebwt2indel"
    
    def _convert_to_fasta(self, output_fasta, output_pairs):
        """
        从 ebwt2InDel FASTA 文件生成:
          - FASTA 文件（保持原样）
          - pairs 文件（每行: path1_name path2_name type）
        """
        # 复制原始 FASTA
        shutil.copy(self.bubble_file, output_fasta)
        
        pairs = defaultdict(lambda: [[], ""])
        
        header_pattern = re.compile(
            r'^cluster:(\d+)_id:(\d+)_right:\d+_cov:\d+_type:_(SNP|INDEL)_event:'
        )
        
        with open(self.bubble_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    name = line.split()[0][1:]
                    match = header_pattern.match(name)
                    
                    if match:
                        cluster_id = int(match.group(1))
                        btype = match.group(3)
                        
                        pairs[cluster_id][0].append(name)
                        pairs[cluster_id][1] = btype
                    else:
                        logger.warning(f"无法解析 header: {name}")
                        continue

        self.write_pairs(output_pairs, pairs)

# ============================================================
# 工具注册表
# ============================================================
TOOL_REGISTRY = {
    "fast2snp": Fast2snpEvaluator,
    "kmer2snp": Kmer2snpEvaluator,
    "discosnp": DiscoSnpEvaluator,
    "discoSnp": DiscoSnpEvaluator,
    "ebwt2indel": Ebwt2IndelEvaluator,
}


# ============================================================
# 主函数
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="通用 bubble 位置评估框架"
    )
    parser.add_argument("--tool", default="fast2snp",
                        choices=list(TOOL_REGISTRY.keys()),
                        help="工具名称")
    parser.add_argument("--input", required=True,
                        help="工具输出的 bubble 文件")
    parser.add_argument("--ref", required=True,
                        help="参考基因组 FASTA")
    parser.add_argument("--truth", required=True,
                        help="标准 VCF 真值")
    parser.add_argument("--bed", required=True,
                        help="高置信区间 BED")
    parser.add_argument("--output-prefix", default="eval_result",
                        help="输出文件前缀")
    parser.add_argument("--mapq", type=int, default=20,
                        help="最低 MAPQ 阈值")
    parser.add_argument("--threads", type=int, default=8,
                        help="BWA 线程数")
    parser.add_argument("--work-dir", default=None,
                        help="工作目录（保存中间文件）")
    parser.add_argument("--keep-bam", action="store_true",
                        help="保留中间 BAM 文件")
    args = parser.parse_args()
    
    for f in [args.input, args.ref, args.truth, args.bed]:
        if not Path(f).is_file():
            logger.error(f"找不到文件: {f}")
            sys.exit(1)
    
    evaluator_class = TOOL_REGISTRY[args.tool]
    evaluator = evaluator_class(
        bubble_file=args.input,
        ref_fasta=args.ref,
        truth_vcf=args.truth,
        bed_file=args.bed,
        min_mapq=args.mapq,
        threads=args.threads,
        work_dir=args.work_dir,
        keep_bam=args.keep_bam
    )
    
    try:
        results = evaluator.evaluate()
        
        logger.info("="*60)
        logger.info(f"评估结果（工具: {args.tool}）")
        logger.info("="*60)
        
        for vtype in ["SNP", "INDEL"]:
            r = results[vtype]
            precision = r['tp'] / (r['tp'] + r['fp']) if (r['tp'] + r['fp']) > 0 else 0
            recall = r['tp'] / (r['tp'] + r['fn']) if (r['tp'] + r['fn']) > 0 else 0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
            
            logger.info(f"--- {vtype} ---")
            logger.info(f"总 {vtype} Bubbles: {r['total_bubbles']}")
            logger.info(f"正确 Bubbles (TP): {r['tp']}")
            logger.info(f"错误 Bubbles (FP): {r['fp']}")
            logger.info(f"未覆盖真值变异 (FN): {r['fn']}")
            logger.info(f"覆盖的真值变异: {len(r['covered_truth'])}")
            logger.info(f"真值变异总数: {r['total_truth']}")
            logger.info(f"Precision: {precision:.4f}")
            logger.info(f"Recall:    {recall:.4f}")
            logger.info(f"F1 Score:  {f1:.4f}")
        
        logger.info("="*60)
        
        for vtype in ["SNP", "INDEL"]:
            r = results[vtype]
            summary_file = f"{args.output_prefix}_{args.tool}_{vtype}_summary.txt"
            with open(summary_file, 'w') as f:
                f.write("Bubble_Idx\tPair\tType\tCorrect\n")
                for i in range(len(evaluator.pairs)):
                    names, bt = evaluator.pairs[i]
                    if bt != vtype:
                        continue
                    correct = i in r['tp_bubbles']
                    f.write(f"{i}\t{','.join(names)}\t{bt}\t{correct}\n")
            
            fn_file = f"{args.output_prefix}_{args.tool}_{vtype}_fn_variants.txt"
            with open(fn_file, 'w') as f:
                f.write("Chrom\tPos\n")
                for chrom, pos in sorted(r['fn_positions']):
                    f.write(f"{chrom}\t{pos}\n")
        
        logger.info("详细结果已保存:")
        for vtype in ["SNP", "INDEL"]:
            logger.info(f"  {args.output_prefix}_{args.tool}_{vtype}_summary.txt")
            logger.info(f"  {args.output_prefix}_{args.tool}_{vtype}_fn_variants.txt")
    
    finally:
        evaluator.cleanup()


if __name__ == "__main__":
    main()