#!/usr/bin/env python3
"""
merge_coverage_pngs.py

将不同覆盖度目录（cov5, cov10, ..., cov60）下的 PNG 图片合并为一张多子图。
每个子图以覆盖度作为标题，可通过参数控制子图大小和字体大小。

用法:
  python merge_coverage_pngs.py --root .. --output merged_coverage.png --cols 3 --scale 8 --title-size 28
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def parse_args():
    parser = argparse.ArgumentParser(description="合并覆盖度目录下的 PNG 图片")
    parser.add_argument("--root", default="..", help="覆盖度目录的父目录（默认 ..）")
    parser.add_argument("--png", default="*.png", help="每个目录下的 PNG 文件名模式（默认 *.png）")
    parser.add_argument("--output", default="merged_coverage.png", help="输出图片路径")
    parser.add_argument("--cols", type=int, default=3, help="子图列数（默认 3）")
    parser.add_argument("--scale", type=float, default=8.0,
                        help="每个子图的边长（英寸，默认 8.0，越大图片越大）")
    parser.add_argument("--dpi", type=int, default=300, help="输出分辨率")
    parser.add_argument("--title-size", type=int, default=24, help="子图标题字号")
    return parser.parse_args()

def main():
    args = parse_args()

    coverages = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    root = Path(args.root)

    images = []
    for cov in coverages:
        dir_path = root / f"{cov}"
        if not dir_path.is_dir():
            print(f"警告: 目录不存在，跳过: {dir_path}")
            continue
        png_files = sorted(dir_path.glob(args.png))
        if not png_files:
            print(f"警告: 目录中没有匹配的 PNG，跳过: {dir_path}")
            continue
        images.append((cov, png_files[0]))
        print(f"cov{cov}: {png_files[0].name}")

    if not images:
        raise SystemExit("错误: 没有找到任何图片。")

    n = len(images)
    cols = args.cols
    rows = (n + cols - 1) // cols

    # 每个子图的尺寸由 --scale 控制，整体图也随之放大
    fig, axes = plt.subplots(rows, cols,
                             figsize=(args.scale * cols, args.scale * rows))
    axes = axes.flatten() if n > 1 else [axes]

    for idx, (cov, img_path) in enumerate(images):
        ax = axes[idx]
        img = mpimg.imread(img_path)
        ax.imshow(img)
        ax.set_title(f"Coverage {cov}×",
                     fontsize=args.title_size, pad=0)
                     #fontweight='bold')  # 标题加粗更醒目
        ax.axis("off")

    for j in range(n, len(axes)):
        axes[j].axis("off")

    plt.tight_layout(pad=0)  # pad 控制子图间距，值越大间距越大
    #plt.tight_layout(pad=0.1, h_pad=0.5, w_pad=0.1)
    plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    print(f"合并完成，已保存至 {args.output}")

if __name__ == "__main__":
    main()