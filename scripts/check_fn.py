#!/usr/bin/env python3
import sys
import csv
from collections import defaultdict
from Bio import SeqIO

fn_vcf_file = sys.argv[1]
std_seed_file = sys.argv[2]
seed_file = sys.argv[3]
def load_std_seed(std_seed_file):
    std_seed = defaultdict(lambda: [None, None])
    dict_reader = csv.DictReader(open(std_seed_file))
    for row in dict_reader:
        chrom = row['Chrom']
        pos = int(row['Variant_Pos'])
        if row['Side'] == 'Left':
            std_seed[(chrom, pos)][0] = row['Kmer']
        else:
            assert row['Side'] == 'Right'
            std_seed[(chrom, pos)][1] = row['Kmer']
    return std_seed


def iter_fn_vcf(fn_vcf_file):
    for line in open(fn_vcf_file):
        if line.startswith("#"):  continue
        its = line.split()
        chrom = its[0]
        pos = int(its[1])
        ref = its[3]
        alt = its[4]
        yield (chrom, pos, ref, alt)

def load_seed(seed_file, k):
    seeds = set()
    for record in SeqIO.parse(seed_file, "fastq"):
        seeds.add(record.seq[:k])
    return seeds

std_seeds = load_std_seed(std_seed_file)
seeds = load_seed(seed_file, 16)

def reverse_complement(kmer):
    m = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

    rv = [m[k] for k in kmer]
    return "".join(rv[::-1])


for chrom, pos, ref, alt in iter_fn_vcf(fn_vcf_file):
    if (chrom, pos) not in std_seeds:
        print("nopos", chrom, pos, ref, alt)
    else:
        left_seed, right_seed = std_seeds[(chrom, pos)]
        if left_seed and left_seed in seeds or right_seed and reverse_complement(right_seed) in seeds:
            print("xxseed", chrom, pos, ref, alt)
        else:
            print("noseed", chrom, pos, ref, alt)



