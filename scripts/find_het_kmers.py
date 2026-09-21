#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import scipy
from scipy.stats import nbinom
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from functools import partial
import math
from pathlib import Path
import os
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="GenomeScope-like model fitting for k-mer frequency histograms.")
    parser.add_argument("input", help="Path to the k-mer histogram file.")

    # Optional parameters
    parser.add_argument("-o", "--output", help="Path to save the output plot.")
    parser.add_argument("-est", "--est_kmercov", type=int, help="Estimated k-mer coverage (optional).")
    parser.add_argument("-het", "--het_ratio", type=float, help="Initial heterozygosity(0/0.12) (optional).")
    parser.add_argument("-rho", "--rho", type=float, default=0.2, help="Initial rho value (0.2)(optional).")
    parser.add_argument("-k", "--kmer", type=int, default=21, help="kmer(21)(optional)")
    parser.add_argument("--only_plot_histo", action="store_true", help="Only plot the histogram (optional).")

    return parser.parse_args()


def load_kmer_histogram_all(hist_file):
    """Load k-mer frequency histogram file."""
    hist_df = pd.read_csv(hist_file, sep='\\s+', header=None, names=['coverage', 'frequency'], skiprows=0, nrows=1000)
    coverage = hist_df['coverage'].values
    frequency = hist_df['frequency'].values
    return coverage, frequency


def negative_binomial_pmf(x, mu, size):
    """Negative binomial probability mass function (PMF)."""
    return nbinom.pmf(x, n=size, p=size / (size + mu))


def genomescope_diploid_model_nb(x, r, d, lambda_val, rho, G, kmer_length=21):
    """GenomeScope diploid mixture model function."""
    alpha_val = 2 * (1 - d) * (1 - (1 - r) ** kmer_length) + 2 * d * (1 - (1 - r) ** kmer_length) ** 2 + 2 * d * (
                (1 - r) ** kmer_length) * (1 - (1 - r) ** kmer_length)
    beta_val = (1 - d) * ((1 - r) ** kmer_length) + d * (1 - (1 - r) ** kmer_length) ** 2
    term1 = alpha_val * negative_binomial_pmf(x, lambda_val, lambda_val / rho)
    term2 = beta_val * negative_binomial_pmf(x, 2 * lambda_val, 2 * lambda_val / rho)
    return G * (term1 + term2)

def genomescope_diploid_model_nb1(x, r, d, lambda_val, rho, G, kmer_length=21):
    """GenomeScope diploid mixture model function (revised)."""
    alpha_val = 2 * (1 - d) * (1 - (1 - r) ** kmer_length) + 2 * d * (1 - (1 - r) ** kmer_length) ** 2 + 2 * d * (
                (1 - r) ** kmer_length) * (1 - (1 - r) ** kmer_length)
    term1 = alpha_val * negative_binomial_pmf(x, lambda_val, lambda_val / rho)
    return G * term1


def genomescope_diploid_model_nb2(x, r, d, lambda_val, rho, G, kmer_length=21):
    """GenomeScope diploid mixture model function (revised)."""
    beta_val = (1 - d) * ((1 - r) ** kmer_length) + d * (1 - (1 - r) ** kmer_length) ** 2
    term2 = beta_val * negative_binomial_pmf(x, 2 * lambda_val, 2 * lambda_val / rho)
    return G * term2


def estimate_het_kmer_cov(coverage_data, frequency_data):
    """Estimate estKmerCov and estGenomeSize."""

    cumulative_sum = np.cumsum(frequency_data)
    total_sum = cumulative_sum[-1]

    
    from scipy.signal import find_peaks
    peak_indices, _ = find_peaks(frequency_data, prominence=0.1 * np.mean(frequency_data), distance=5)
    max_peak = max(frequency_data[peak_indices]) if len(peak_indices) > 0 else 0
    logging.info(f"Peak indices found: {[coverage_data[idx] for idx in peak_indices]} {[frequency_data[idx] for idx in peak_indices]}")
    peak_indices = [idx for idx in peak_indices if frequency_data[idx] >= 0.1 * max_peak]


    if len(peak_indices) == 1:
        peak2 = coverage_data[peak_indices[0]]
        peak1 = coverage_data[peak_indices[0]] / 2
        est_kmercov = int(peak1)
    elif len(peak_indices) == 2:
        peak1, peak2 = coverage_data[peak_indices[0]], coverage_data[peak_indices[1]]
        ratio = peak2 / peak1
        est_kmercov = peak1
    elif len(peak_indices) > 2:
        top_two_peak_indices = sorted(
            sorted(peak_indices, key=lambda idx: frequency_data[idx], reverse=True)[:2]
        )
        
        peak1, peak2 = coverage_data[top_two_peak_indices[0]], coverage_data[top_two_peak_indices[1]]
        est_kmercov = peak1
        ratio = peak2 / peak1

    else:
        est_kmercov = 1 
        print(f"The coverage is too low, the selection of hybrid peak may be wrong !!")

    logging.info(f"Estimated k-mer coverage: {est_kmercov}")
    return est_kmercov, peak1, peak2


def plot_kmer_histogram(coverage_all, frequency_all, output_file):
    """Plot k-mer frequency histogram."""
    plt.figure(figsize=(10, 6))
    plt.bar(coverage_all[:int(3 * 5)], frequency_all[:int(3 * 5)], label='k-mer spectrum',
            width=0.9, alpha=0.7)
    plt.gca().set_facecolor('lightgray')
    plt.ylim(0, max(frequency_all[2*1]*2, frequency_all[2*2]*2))
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(output_file)

def main():
# Parse arguments
    args = parse_arguments()
    input_file = args.input
    file_name = Path(input_file).stem

    # Heterozygosity bounds
    min_het, max_het = 0, 1

    # Load k-mer data
    coverage_all, frequency_all = load_kmer_histogram_all(input_file)

    # 全局字体设置：进一步增大字体
    plt.rcParams.update({
        'font.size': 18,
        'axes.titlesize': 20,      # 子图标题
        'axes.labelsize': 18,      # 坐标轴标签
        'xtick.labelsize': 16,     # x轴刻度
        'ytick.labelsize': 16,     # y轴刻度
        'legend.fontsize': 18,     # 图例
        'font.sans-serif': ['SimHei', 'Arial Unicode MS', 'DejaVu Sans'],
        'axes.unicode_minus': False
    })

    if args.only_plot_histo:
        plot_kmer_histogram(coverage_all, frequency_all, f"{file_name}_histo.png")
        return
    # Estimate initial parameters

    # Use provided parameters or defaults
    est_kmercov, peak1, peak2 = (args.est_kmercov, args.est_kmercov, args.est_kmercov*2) if args.est_kmercov else estimate_het_kmer_cov(coverage_all, frequency_all)
    print(f"Estimated k-mer coverage: {est_kmercov}, Peak1: {peak1}, Peak2: {peak2}")
    if est_kmercov == 1:
        logging.error("no peaks, maybe too high error rate")
        raise SystemExit(1)
    
    initial_het = args.het_ratio if args.het_ratio is not None else (min_het + 1e-3)
    initial_rho = args.rho 
    k = args.kmer 
    est_genomesize = np.sum(coverage_all * frequency_all) / (2 * est_kmercov)

    print(f"Using parameters: est_kmercov={est_kmercov}, initial_het={initial_het}, initial_rho={initial_rho}")


    # Fit the model
    model_func = partial(genomescope_diploid_model_nb, kmer_length=k)
    optimal_params, _ = curve_fit(
        model_func,
        coverage_all[math.floor(est_kmercov / 2):],
        frequency_all[math.floor(est_kmercov / 2):],
        p0=[initial_het, 0, est_kmercov, initial_rho, est_genomesize],
        bounds=([0, 0, 0, 0.15, 0], [1, 1, np.inf, 0.65, np.inf])
    )
    r_opt, d_opt, lambda_opt, rho_opt, G_opt = optimal_params

    print("\nFitting successful!")
    print(f"r (heterozygosity): {r_opt:.4f}")
    print(f"d (repeat fraction): {d_opt:.4f}")
    print(f"lambda (haploid coverage): {lambda_opt:.2f}")
    print(f"G (genome scaling factor): {G_opt:.2e}")
    print(f"rho: {rho_opt:.2e}")

    # Central 90% interval (5%-95%) for the heterozygous peak NB component.
    het_mu = lambda_opt
    het_size = lambda_opt / rho_opt
    het_p = het_size / (het_size + het_mu)
    het_left = int(np.floor(nbinom.ppf(0.05, n=het_size, p=het_p)))
    het_right = int(np.ceil(nbinom.ppf(0.95, n=het_size, p=het_p)))
    het_left = max(0, het_left)
    het_right = max(het_left, het_right)
    print(f"heterozygous peak 90% interval: {het_left} - {het_right}")

    # Plot results
    x_plot = np.linspace(1, int(lambda_opt * 5), int(lambda_opt * 5))
    y_predicted1 = genomescope_diploid_model_nb1(x_plot, *optimal_params)
    y_predicted2 = genomescope_diploid_model_nb2(x_plot, *optimal_params)
    y_predicted = genomescope_diploid_model_nb(x_plot, *optimal_params)

    # Compute residuals in the range (0, est_kmercov)
    mask_diff = (coverage_all >= 0) & (coverage_all <= est_kmercov)
    x_diff = coverage_all[mask_diff]
    observed_diff = frequency_all[mask_diff]
    predicted_diff = genomescope_diploid_model_nb(x_diff, *optimal_params)
    residual = observed_diff - predicted_diff
    residual = np.clip(residual, 0, None)  # Replace all negative values with 0

    total_error_kmers = np.sum(residual * x_diff)
    total_kmers = np.sum(coverage_all * frequency_all)
    error_rate = 1 - (1 - (total_error_kmers / total_kmers)) ** (1 / k)
    genome_length = (total_kmers - total_error_kmers) / (lambda_opt * 2)

    length_diff = len(y_predicted) - len(residual)
    residual_padded = np.pad(residual, (0, length_diff), mode='constant', constant_values=0)
    r2 = r2_score(frequency_all[int(lambda_opt):int(lambda_opt * 5)], y_predicted[int(lambda_opt):])

    left = max(2, int(np.floor(-est_kmercov/2 + est_kmercov)))
    right = int(np.ceil(est_kmercov/2 + est_kmercov))

    print(f"\nleft,right: {left} - {right}")
    print(f"error_rate: {error_rate:.4f}")
    print(f"genome_length: {genome_length:.2e}")

    plt.figure(figsize=(10, 6))
    plt.bar(coverage_all[:int(lambda_opt * 5)], frequency_all[:int(lambda_opt * 5)], label='k-mer spectrum',
            width=0.9, alpha=0.7)
    #plt.axvline(x=max_coverage_95, color='r', linestyle='--', linewidth=2, label='cut-off threshold')
    plt.axvline(x=peak1, color='black', linestyle='--', linewidth=1, label='peaks')
    plt.axvline(x=peak2, color='black', linestyle='--', linewidth=1)
#    plt.plot(x_plot, y_predicted1, '-', label='heterozygous', linewidth=2, color='orange')
#    plt.plot(x_plot, y_predicted2, '-', label='homozygous', linewidth=2, color='yellow')
#    if len(x_diff) > 1:
#        plt.plot(x_diff, residual, label='erroneous', linewidth=2, color='green')
#    elif len(x_diff) == 1:
#        plt.vlines(x=x_diff[0], ymin=0, ymax=residual[0], colors='green', linestyles='--', linewidth=2,
#                   label='Residual')
#    plt.plot(x_plot, y_predicted+residual_padded, '-', label='full model', linewidth=3, color='black')

    # title = (f'Diploid Model Fitting\n'
    #          '\n'
    #          f'len: {genome_length:.2e}    '
    #          f'r: {r_opt*100:.3f}%    '
    #          f'err: {error_rate*100:.3f}%    '
    #          f'kcov: {lambda_opt:.1f}    '
    #          f'R^2: {r2:.2f}    '
    #          f'range: {left}-{right}\n'
    #          f'het90: {het_left}-{het_right}'
    #          )

    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    #plt.title(title)
    plt.ylim(0, max(frequency_all[est_kmercov*2]*2,frequency_all[est_kmercov*1]*2))
    plt.gca().set_facecolor('lightgray')
    plt.axvspan(left, right, color='purple', alpha=0.3, label='[$lowCov$, $highCov$]')
    plt.legend()
    #plt.axvspan(het_left, het_right, color='purple', alpha=0.15, label='Het 90% interval')
    output_file = f"{file_name}_{left}_{right}.png"
    if args.output:
        output_file = os.path.join(args.output, output_file)
    plt.savefig(output_file)
    print(f"Graph saved to {output_file}")

if __name__ == "__main__":
    main()
    