#!/usr/bin/env python3
"""
eval_common.py

公共基类与工具提取器，供 position 评估和 hap.py 评估共用。
包含：日志配置、BaseEvaluator、ToolExtractor 及其子类。
"""

import logging
import subprocess
import tempfile
import shutil
import os
from pathlib import Path
from collections import defaultdict
from abc import ABC, abstractmethod
import re
import pysam

# ---------- 日志 ----------
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


# ============================================================
# 工具提取器
# ============================================================
class ToolExtractor(ABC):
    """各工具只需实现 extract()，将自身输出转换为 FASTA + pairs"""
    def __init__(self, bubble_file):
        self.bubble_file = bubble_file

    @abstractmethod
    def extract(self, output_fasta, output_pairs):
        pass


def write_pairs_to_file(pairs, output_pairs):
    """pairs: dict {id: [[names], type]}，将完整配对写入文件"""
    type_counts = defaultdict(int)
    with open(output_pairs, 'w') as f:
        for names, btype in pairs.values():
            if len(names) > 1 and btype:
                type_counts[btype] += 1
                f.write(f"{' '.join(names)} {btype}\n")
            else:
                pass
                # logger.warning(f"Bubble 配对不完整: size={len(names)}, paths='{names}', type='{btype}'")
    logger.info(f"  共转换 {len(pairs)} 对 bubble")
    logger.info(f"    SNP: {type_counts.get('SNP', 0)}")
    logger.info(f"    INDEL: {type_counts.get('INDEL', 0)}")


class Fast2snpExtractor(ToolExtractor):
    tool_name = "fast2snp"
    def extract(self, output_fasta, output_pairs):
        shutil.copy(self.bubble_file, output_fasta)
        pairs = defaultdict(lambda: [[], ""])
        with open(self.bubble_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    name = line.split()[0][1:]
                    m = re.match(r'Bubble_(\d+)_(INDEL|IsolatedSNP|NonIsolatedSNP)_Path([12])_Ext(\S+)', name)
                    if m:
                        bid = m.group(1)
                        btype = 'SNP' if 'SNP' in m.group(2).upper() else 'INDEL'
                        path = int(m.group(3))
                        pairs[bid][0].append(name)
                        pairs[bid][1] = btype
                    else:
                        logger.warning(f"无法解析 header: {name}")
        write_pairs_to_file(pairs, output_pairs)


class Kmer2snpExtractor(ToolExtractor):
    tool_name = "kmer2snp"
    def extract(self, output_fasta, output_pairs):
        pairs = defaultdict(lambda: [[], ""])
        with open(self.bubble_file) as inf, open(output_fasta, 'w') as fo:
            for i, line in enumerate(inf):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                name1 = f"kmer2snp_{i}_path1"
                name2 = f"kmer2snp_{i}_path2"
                fo.write(f">{name1}\n{parts[0]}\n>{name2}\n{parts[1]}\n")
                pairs[i][0] = [name1, name2]
                pairs[i][1] = "SNP"
        write_pairs_to_file(pairs, output_pairs)


class DiscoSnpExtractor(ToolExtractor):
    tool_name = "discoSnp"
    def extract(self, output_fasta, output_pairs):
        pairs = defaultdict(lambda: [[], ""])

        # 打开输入和输出文件
        with open(self.bubble_file, 'r') as fin, open(output_fasta, 'w') as fout:
            current_valid = False   # 当前 header 是否保留
            current_name = None
            current_bid = None
            current_path = None
            current_type = None

            for line in fin:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    # 检查是否包含杂合标记 G1_0/1
                    if 'G1_0/1' in line:
                        current_valid = True
                        # 解析 header 信息
                        name = line.split()[0][1:]
                        segs = name.split('|')[0].split('_')
                        if len(segs) == 4 and segs[0] in ("SNP", "INDEL") and segs[2] == "path":
                            current_bid = segs[3]
                            current_path = 1 if segs[1] == "higher" else 2
                            current_type = segs[0]
                            current_name = name.split('|')[0]
                            # 写入 header 到输出 fasta
                            fout.write(">" + name.split('|')[0] + '\n')
                        else:
                            logger.warning(f"无法解析 header: {name}")
                            current_valid = False
                    else:
                        # 非杂合，跳过该记录
                        current_valid = False
                else:
                    # 序列行
                    if current_valid:
                        # 去除小写碱基，只保留大写
                        cleaned = ''.join(c for c in line if c.isupper())
                        if cleaned and current_name is not None:
                            fout.write(cleaned + '\n')
                            # 记录到 pairs
                            if current_bid is not None and current_path is not None:
                                pairs[current_bid][0].append(current_name)
                                pairs[current_bid][1] = current_type

        # 写 pairs 文件
        write_pairs_to_file(pairs, output_pairs)

class Ebwt2IndelExtractor(ToolExtractor):
    tool_name = "ebwt2indel"
    def extract(self, output_fasta, output_pairs):
        shutil.copy(self.bubble_file, output_fasta)
        pairs = defaultdict(lambda: [[], ""])
        pat = re.compile(r'^cluster:(\d+)_id:(\d+)_right:\d+_cov:\d+_type:_(SNP|INDEL)_event:')
        with open(self.bubble_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    name = line.split()[0][1:]
                    m = pat.match(name)
                    if m:
                        bid = m.group(1)
                        pairs[bid][0].append(name)
                        pairs[bid][1] = m.group(3)
                    else:
                        logger.warning(f"无法解析 header: {name}")
        write_pairs_to_file(pairs, output_pairs)


TOOL_EXTRACTORS = {
    "fast2snp": Fast2snpExtractor,
    "kmer2snp": Kmer2snpExtractor,
    "discosnp": DiscoSnpExtractor,
    "discoSnp": DiscoSnpExtractor,
    "ebwt2indel": Ebwt2IndelExtractor,
}


# ============================================================
# 公共评估基类
# ============================================================
class BaseEvaluator(ABC):
    def __init__(self, tool_name, bubble_file, ref_fasta, truth_vcf, bed_file,
                 min_mapq=20, threads=8, work_dir=None, keep_bam=False):
        self.tool_name = tool_name
        self.extractor = TOOL_EXTRACTORS[tool_name](bubble_file)
        self.ref_fasta = ref_fasta
        self.truth_vcf = truth_vcf
        self.bed_file = bed_file
        self.min_mapq = min_mapq
        self.threads = threads
        self.pairs = []

        if work_dir:
            self.work_dir = Path(work_dir)
            self.work_dir.mkdir(parents=True, exist_ok=True)
            self._temp_dir = None
        else:
            self._temp_dir = tempfile.TemporaryDirectory(prefix="eval_")
            self.work_dir = Path(self._temp_dir.name)

        self.keep_bam = keep_bam
        self.bam_path = None
        self.fasta_path = self.work_dir / "bubble_sequences.fasta"
        self.pairs_path = self.work_dir / "bubble_pairs.txt"

    # ---------- 抽象方法 ----------
    @abstractmethod
    def _evaluate(self, bed_regions, truth_variants):
        pass

    @abstractmethod
    def _print_results(self, results):
        pass

    # ---------- 公共流程 ----------
    def prepare(self):
        logger.info(f"[{self.tool_name}] 提取 bubble 文件...")
        self.extractor.extract(self.fasta_path, self.pairs_path)

        logger.info("加载 bubble pairs...")
        self.pairs = self._load_pairs()
        logger.info(f"  共 {len(self.pairs)} 个 bubbles")

        self.bam_path = self._align_to_reference()

        bed_regions = self._load_bed() if self.bed_file else None
        truth_variants = self._load_truth(bed_regions)
        return bed_regions, truth_variants

    def run(self):
        bed_regions, truth_variants = self.prepare()
        results = self._evaluate(bed_regions, truth_variants)
        self._print_results(results)
        return results

    def cleanup(self):
        if self._temp_dir and not self.keep_bam:
            self._temp_dir.cleanup()

    # ---------- 辅助方法 ----------
    def _load_pairs(self):
        pairs = []
        with open(self.pairs_path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    pairs.append((parts[:-1], parts[-1]))
        return pairs

    def _align_to_reference(self):
        logger.info("检查 BWA 索引...")
        if not Path(self.ref_fasta + ".bwt").exists():
            subprocess.run(["bwa", "index", self.ref_fasta],
                           check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sam = self.work_dir / "alignments.sam"
        with open(sam, 'w') as f:
            subprocess.run(["bwa", "mem", "-t", str(self.threads), "-k", "13",
                            "-T 15 -B 3 -O 5 ", self.ref_fasta, str(self.fasta_path)],
                           stdout=f, check=True)
        unsorted_bam = self.work_dir / "alignments.unsorted.bam"
        sorted_bam = self.work_dir / "alignments.sorted.bam"
        subprocess.run(["samtools", "view", "-bS", str(sam), "-o", str(unsorted_bam)], check=True)
        subprocess.run(["samtools", "sort", str(unsorted_bam), "-o", str(sorted_bam)], check=True)
        subprocess.run(["samtools", "index", str(sorted_bam)], check=True)
        if not self.keep_bam:
            os.remove(sam)
            os.remove(unsorted_bam)
        self.bam_path = str(sorted_bam)
        logger.info(f"  BAM: {self.bam_path}")
        return self.bam_path

    def _load_bed(self):
        regions = defaultdict(list)
        with open(self.bed_file) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    regions[parts[0]].append((int(parts[1]), int(parts[2])))
        return regions

    def _is_in_bed(self, chrom, pos1, bed_regions):
        if chrom not in bed_regions:
            return False
        pos0 = pos1 - 1
        return any(start <= pos0 < end for start, end in bed_regions[chrom])

    def _load_truth(self, bed_regions):
        snp_pos = set()
        indel_pos = set()
        vcf = pysam.VariantFile(self.truth_vcf)
        for rec in vcf:
            if bed_regions and not self._is_in_bed(rec.chrom, rec.pos, bed_regions):
                continue
            for alt in rec.alts:
                if alt is None or alt in ('<NON_REF>', '<*>'):
                    continue
                if len(rec.ref) == 1 and len(alt) == 1:
                    snp_pos.add((rec.chrom, rec.pos))
                else:
                    indel_pos.add((rec.chrom, rec.pos))
        vcf.close()
        return {"SNP": snp_pos, "INDEL": indel_pos}

    def _get_query_alignments(self):
        query_to_align = defaultdict(list)
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.mapping_quality < self.min_mapq:
                continue
            query_to_align[read.query_name].append(
                (read.reference_name, read.reference_start, read.reference_end)
            )
        bam.close()
        return query_to_align