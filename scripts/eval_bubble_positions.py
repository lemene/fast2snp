#!/usr/bin/env python3
"""
eval_bubble_positions.py

位置评估：检查 bubble 在参考基因组上的比对位置是否覆盖真实变异。
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
from eval_common import BaseEvaluator, TOOL_EXTRACTORS, logger


class PositionEvaluator(BaseEvaluator):
    def _evaluate(self, bed_regions, truth_variants):
        logger.info("构建位置索引...")
        query_to_align = self._get_query_alignments()
        pos_index = {"SNP": defaultdict(set), "INDEL": defaultdict(set)}

        for bubble_idx, (names, btype) in enumerate(self.pairs):
            aligns = []
            for name in names:
                if name in query_to_align:
                    aligns.extend(query_to_align[name])
            target = pos_index[btype]
            for ref_name, start, end in aligns:
                for pos0 in range(start, end):
                    target[(ref_name, pos0)].add(bubble_idx)

        results = {}
        for vtype in ["SNP", "INDEL"]:
            type_indices = {i for i, (_, bt) in enumerate(self.pairs) if bt == vtype}
            truth_positions = truth_variants[vtype]
            tp_bubbles = set()
            covered_truth = set()
            for (chrom, pos1) in truth_positions:
                pos0 = pos1 - 1
                hits = pos_index[vtype].get((chrom, pos0), set())
                if hits:
                    covered_truth.add((chrom, pos1))
                    tp_bubbles.update(hits)
            results[vtype] = {
                'tp': len(tp_bubbles),
                'fp': len(type_indices - tp_bubbles),
                'fn': len(truth_positions - covered_truth),
                'tp_bubbles': tp_bubbles,
                'fn_positions': truth_positions - covered_truth,
                'covered_truth': covered_truth,
                'total_bubbles': len(type_indices),
                'total_truth': len(truth_positions),
            }
        return results

    def _print_results(self, results):
        logger.info("=" * 60)
        logger.info(f"位置评估结果（工具: {self.tool_name}）")
        logger.info("=" * 60)
        for vtype in ["SNP", "INDEL"]:
            r = results[vtype]
            prec = r['tp'] / (r['tp'] + r['fp']) if r['tp'] + r['fp'] else 0
            rec = r['tp'] / (r['tp'] + r['fn']) if r['tp'] + r['fn'] else 0
            f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0
            logger.info(f"--- {vtype} ---")
            logger.info(f"  Total bubbles: {r['total_bubbles']}")
            logger.info(f"  TP={r['tp']}, FP={r['fp']}, FN={r['fn']}")
            logger.info(f"  Precision={prec:.4f}, Recall={rec:.4f}, F1={f1:.4f}")


def main():
    parser = argparse.ArgumentParser(description="位置评估")
    parser.add_argument("--tool", default="fast2snp", choices=list(TOOL_EXTRACTORS.keys()))
    parser.add_argument("--input", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--bed", default="")
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--mapq", type=int, default=20)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--keep-bam", action="store_true")
    args = parser.parse_args()

    for f in [args.input, args.ref, args.truth]:
        if not Path(f).is_file():
            sys.exit(f"错误: 找不到文件 '{f}'")

    ev = PositionEvaluator(
        tool_name=args.tool,
        bubble_file=args.input,
        ref_fasta=args.ref,
        truth_vcf=args.truth,
        bed_file=args.bed,
        min_mapq=args.mapq,
        threads=args.threads,
        work_dir=args.outdir,
        keep_bam=args.keep_bam,
    )
    try:
        results = ev.run()
        # 保存结果
        for vtype in ["SNP", "INDEL"]:
            r = results[vtype]

            summary = Path(args.outdir) / f"{args.tool}_{vtype}_summary.txt"
            with open(summary, 'w') as f:
                f.write("Bubble_Idx\tPair\tType\tCorrect\n")
                for i, (names, bt) in enumerate(ev.pairs):
                    if bt == vtype:
                        f.write(f"{i}\t{','.join(names)}\t{bt}\t{i in r['tp_bubbles']}\n")
            fn_file = Path(args.outdir) / f"{args.tool}_{vtype}_fn_variants.txt"
            with open(fn_file, 'w') as f:
                f.write("Chrom\tPos\n")
                for chrom, pos in sorted(r['fn_positions']):
                    f.write(f"{chrom}\t{pos}\n")
        logger.info("结果文件已保存。")
    finally:
        ev.cleanup()


if __name__ == "__main__":
    main()