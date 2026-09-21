#!/usr/bin/env python3
"""
eval_with_happy.py

使用 hap.py 评估：从 BAM 提取变异生成 VCF，运行 hap.py。
"""

import argparse
import sys
import subprocess
from pathlib import Path
from eval_common import BaseEvaluator, TOOL_EXTRACTORS, logger


class HappyEvaluator(BaseEvaluator):
    def __init__(self, *args, hap_py_bin="hap.py", **kwargs):
        super().__init__(*args, **kwargs)
        self.hap_py_bin = hap_py_bin
        self.query_vcf = self.work_dir / "query.vcf"

    def _evaluate(self, bed_regions, truth_variants):
        logger.info("从 BAM 生成 query VCF...")
        self._generate_vcf_from_bam()
        logger.info("运行 hap.py ...")
        summary = self._run_happy()
        return {"summary_csv": summary}

    def _generate_vcf_from_bam(self):
        # 使用 extract_variants_from_bam.py 从 BAM 直接提取所有差异
        script = Path(__file__).parent / "extract_variants_from_bam.py"
        subprocess.run([
            sys.executable, str(script),
            "--bam", self.bam_path,
            "--ref", self.ref_fasta,
            "--output", str(self.query_vcf),
        ], check=True)


    def _run_happy(self):
        script = Path(__file__).parent / "run_singularity_happy.py"
        out_prefix = self.work_dir / "happy_result"
        cmd = [
            str(script),
            "-t", self.truth_vcf,
            "-q", str(self.query_vcf),
            "-r", self.ref_fasta,
            "-o", str(out_prefix),
            "-p", f"{self.threads}",
        ]
        if self.bed_file:
            cmd += ["-b", self.bed_file]

        logger.info(f"  { ' '.join(cmd) }")
        subprocess.run(cmd, check=True)
        return Path(str(out_prefix) + ".summary.csv")

    def _print_results(self, results):
        logger.info(f"hap.py 评估完成，summary: {results['summary_csv']}")


def main():
    parser = argparse.ArgumentParser(description="hap.py 评估")
    parser.add_argument("--tool", default="fast2snp", choices=list(TOOL_EXTRACTORS.keys()))
    parser.add_argument("--input", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--bed")
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--mapq", type=int, default=20)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--keep-bam", action="store_true")
    parser.add_argument("--hap-py-bin", default="hap.py")
    args = parser.parse_args()

    for f in [args.input, args.ref, args.truth, args.bed]:
        if f and not Path(f).is_file():
            sys.exit(f"错误: 找不到文件 '{f}'")

    ev = HappyEvaluator(
        tool_name=args.tool,
        bubble_file=args.input,
        ref_fasta=args.ref,
        truth_vcf=args.truth,
        bed_file=args.bed,
        min_mapq=args.mapq,
        threads=args.threads,
        work_dir=args.outdir,
        keep_bam=args.keep_bam,
        hap_py_bin=args.hap_py_bin,
    )
    try:
        results = ev.run()
        # 可复制 summary 到输出前缀
        import shutil
        if Path(results["summary_csv"]).is_file():
            shutil.copy(results["summary_csv"], f"{args.output_prefix}_{args.tool}_happy_summary.csv")
            logger.info(f"summary 已复制到 {args.output_prefix}_{args.tool}_happy_summary.csv")
    finally:
        ev.cleanup()


if __name__ == "__main__":
    main()