#!/usr/bin/env python3

import os
from collections import defaultdict
import glob

sss = defaultdict(list)

for i, line in enumerate(open("read_kmer_freq.txt")):
    if i == 0: continue

    its = line.split()
    sss[its[0]].append(int(its[-1]))

# get low high parameter
png_fname = glob.glob("../fast2snp/kmer_analysis/*.png")
assert (len(png_fname) == 1)
first_png = os.path.basename(png_fname[0])
stem = os.path.splitext(first_png)[0]
parts = stem.split("_")
low, high = int(parts[2]), int(parts[3])
def in_range(x):
    return x >= low and x <= high

snp_counts = [0, 0, 0, 0]
indel_counts = [0, 0, 0, 0]
mix_counts = [0, 0, 0, 0]
snp_het_radio = []
indel_het_radio = []
mix_het_radio = []
for name, freqs in sss.items():
    if name.find("SNP") < 0:
        indel_counts [0] += 2
        if in_range(freqs[0]):
            indel_counts[1] += 1
        if in_range(freqs[-1]):
            indel_counts[1] += 1

        c = sum(1 for i in freqs if in_range(i))
        indel_counts[2] += len(freqs)
        indel_counts[3] += c
        indel_het_radio.append(c/len(freqs))

    elif name.find("INDEL") < 0:
        snp_counts[0] += 2
        if in_range(freqs[0]):
            snp_counts[1] += 1
        if in_range(freqs[-1]):
            snp_counts[1] += 1

        c = sum(1 for i in freqs if in_range(i))
        snp_counts[2] += len(freqs)
        snp_counts[3] += c
        snp_het_radio.append(c/len(freqs))

    else:
        mix_counts [0] += 2
        if in_range(freqs[0]):
            mix_counts[1] += 1
        if in_range(freqs[-1]):
            mix_counts[1] += 1

        c = sum(1 for i in freqs if in_range(i))
        mix_counts[2] += len(freqs)
        mix_counts[3] += c
        mix_het_radio.append(c/len(freqs))

def print_counts(ccc):
    print(ccc, ccc[1]/ccc[0], ccc[3]/ccc[2])

print_counts(snp_counts)
print_counts(indel_counts)
print_counts(mix_counts)


import matplotlib.pyplot as plt

plt.hist(snp_het_radio, label="snp", bins=100)
plt.hist(indel_het_radio, label="indel", bins=100)
plt.hist(mix_het_radio, label="mixed")
plt.legend()
plt.show()