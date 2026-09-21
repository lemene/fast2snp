#!/usr/bin/env python3
"""
Automatically map local paths into a Singularity container for hap.py.

Usage example:
  python run_happy.py \
      -t /path/to/HG001_truth.vcf.gz \
      -q /path/to/query.vcf \
      -r /path/to/ref.fa \
      -o /path/to/output \
      --prefix my_sample \
      -b /path/to/confident.bed \
      -p 4
"""

import argparse
import subprocess
import sys
from pathlib import Path

IMAGE = Path.home() / "tool" / "happy_docker_image.sif"
HAPPY_BIN = "/opt/hap.py/bin/hap.py"


def main():
    parser = argparse.ArgumentParser(
        description="Run hap.py inside Singularity with automatic path mapping."
    )
    parser.add_argument("-t", "--truth", required=True)
    parser.add_argument("-q", "--query", required=True)
    parser.add_argument("-r", "--ref", required=True)
    parser.add_argument("-o", "--outdir", required=True,
                        help="Output directory")
    parser.add_argument("--prefix", default=None,
                        help="Output file prefix (default: derived from query filename)")
    parser.add_argument("-b", "--bed")
    parser.add_argument("-p", "--threads", type=int, default=1)
    parser.add_argument("--image", default=IMAGE)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("extra", nargs=argparse.REMAINDER)

    args = parser.parse_args()

    truth = Path(args.truth).resolve()
    query = Path(args.query).resolve()
    ref = Path(args.ref).resolve()
    outdir = Path(args.outdir).resolve()
    bed = Path(args.bed).resolve() if args.bed else None
    image = Path(args.image).resolve()

    outdir.mkdir(parents=True, exist_ok=True)

    # 自动生成输出前缀
    if args.prefix:
        prefix = args.prefix
    else:
        # 去掉 query 文件的已知扩展名
        qname = query.name
        for ext in [".vcf.gz", ".vcf", ".bcf"]:
            if qname.endswith(ext):
                qname = qname[:-len(ext)]
                break
        prefix = qname

    # 绑定管理
    bind_flags = []
    used_dirs = {}

    def add_bind(path, label):
        nonlocal bind_flags
        path = path.resolve()
        if path.is_file():
            host_dir = path.parent
            filename = path.name
        else:
            host_dir = path
            filename = None
        if host_dir not in used_dirs:
            mount_point = f"/data_{label}"
            bind_flags.append(f"{host_dir}:{mount_point}")
            used_dirs[host_dir] = mount_point
        else:
            mount_point = used_dirs[host_dir]
        if filename:
            return f"{mount_point}/{filename}"
        else:
            return mount_point

    cont_truth = add_bind(truth, "truth")
    cont_query = add_bind(query, "query")
    cont_ref = add_bind(ref, "ref")
    cont_out = add_bind(outdir, "out")
    cont_bed = add_bind(bed, "bed") if bed else None

    # 输出前缀：容器内路径
    cont_out_prefix = f"{cont_out}/{prefix}"

    # 构建命令
    cmd = ["singularity", "exec"]
    for b in bind_flags:
        cmd.extend(["-B", b])
    cmd.append(str(image))
    cmd.append(HAPPY_BIN)
    cmd.append(cont_truth)
    cmd.append(cont_query)
    cmd.extend(["-r", cont_ref])
    cmd.extend(["-o", cont_out_prefix])
    if cont_bed:
        cmd.extend(["-f", cont_bed])
    cmd.append("--engine=vcfeval")
    cmd.append(f"--threads={args.threads}")
    cmd.append("--pass-only")
    if args.extra:
        cmd.extend(args.extra)

    print("Executing:")
    print(" ".join(cmd))
    print()

    if not args.dry_run:
        subprocess.run(cmd, check=True)
    else:
        print("Dry run: command not executed.")


if __name__ == "__main__":
    main()