#!/usr/bin/env python3
"""
plot_combined_kmer_effect_3x2.py

将三个数据集的 k-mer 长度影响曲线绘制为 3 行 2 列的子图。
行：ONT R10, Illumina, PacBio HiFi
列：SNP, INDEL
每个子图显示 Recall、Precision、F1 随 k 的变化。
每个子图有小标题，整图无总标题，字体较大。

用法:
  python plot_combined_kmer_effect_3x2.py \
      --r10 r10_results.csv \
      --illumina illumina_results.csv \
      --hifi hifi_results.csv \
      --output combined_kmer_effect_3x2.png
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description="绘制 3x2 的 k-mer 影响合并图")
    parser.add_argument("--r10", required=True, help="R10 数据 CSV 文件")
    parser.add_argument("--illumina", required=True, help="Illumina 数据 CSV 文件")
    parser.add_argument("--hifi", required=True, help="HiFi 数据 CSV 文件")
    parser.add_argument("--output", default="combined_kmer_effect_3x2.png", help="输出图片路径")
    return parser.parse_args()

def main():
    args = parse_args()

    df_r10 = pd.read_csv(args.r10)
    df_illumina = pd.read_csv(args.illumina)
    df_hifi = pd.read_csv(args.hifi)

    datasets = [
        ("Illumina", df_illumina),
        ("PacBio HiFi", df_hifi),
        ("ONT R10", df_r10),
    ]

    # 全局字体设置：进一步增大字体
    plt.rcParams.update({
        'font.size': 16,
        'axes.titlesize': 18,      # 子图标题
        'axes.labelsize': 16,      # 坐标轴标签
        'xtick.labelsize': 14,     # x轴刻度
        'ytick.labelsize': 14,     # y轴刻度
        'legend.fontsize': 16,     # 图例
        'font.sans-serif': ['SimHei', 'Arial Unicode MS', 'DejaVu Sans'],
        'axes.unicode_minus': False
    })

    # 创建 3 行 2 列子图，无总标题
    fig, axes = plt.subplots(3, 2, figsize=(18, 16))

    metric_styles = {
        "Recall": {"linestyle": "-", "marker": "s", "color": "black"},
        "Precision": {"linestyle": "--", "marker": "o", "color": "#d62728"},
        "F1": {"linestyle": "-.", "marker": "^", "color": "#1f77b4"}
    }

    # 绘制每个子图
    for row_idx, (dataset_name, df) in enumerate(datasets):
        # SNP 子图（左列）
        ax_snp = axes[row_idx, 0]
        for metric, style in metric_styles.items():
            col_name = f"SNP_PASS_{metric}"
            ax_snp.plot(df['K'], df[col_name],
                        marker=style['marker'],
                        linestyle=style['linestyle'],
                        color=style['color'],
                        label=metric)
        ax_snp.set_xlabel('k-mer Length', fontsize=16)
        ax_snp.set_ylabel('Metric Value', fontsize=16)
        ax_snp.grid(True, linestyle='--', alpha=0.5)
        ax_snp.autoscale(axis='y')
        ax_snp.set_title(f"{dataset_name} - SNP", fontsize=18)

        # INDEL 子图（右列）
        ax_indel = axes[row_idx, 1]
        for metric, style in metric_styles.items():
            col_name = f"INDEL_PASS_{metric}"
            ax_indel.plot(df['K'], df[col_name],
                          marker=style['marker'],
                          linestyle=style['linestyle'],
                          color=style['color'],
                          label=metric)
        ax_indel.set_xlabel('k-mer Length', fontsize=16)
        ax_indel.set_ylabel('Metric Value', fontsize=16)
        ax_indel.grid(True, linestyle='--', alpha=0.5)
        ax_indel.autoscale(axis='y')
        ax_indel.set_title(f"{dataset_name} - INDEL", fontsize=18)

    # 统一图例（放在图下方，字体已加大）
    handles = []
    labels = []
    for metric, style in metric_styles.items():
        line, = plt.plot([], [], marker=style['marker'], linestyle=style['linestyle'],
                         color=style['color'], label=metric)
        handles.append(line)
        labels.append(metric)
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.02),
               ncol=3, fontsize=16)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"图片已保存至: {args.output}")

if __name__ == "__main__":
    main()