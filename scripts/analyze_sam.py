#!/usr/bin/env python3
from collections import defaultdict
import sys 
import pysam
import argparse

def collect_reads_with_overhang(fname):

    CIGAR_S = 4  # Soft clip (软剪切)
    CIGAR_H = 5  # Hard clip (硬剪切)
    def has_overhang(cigar_tuples):
        if not cigar_tuples:
            return False

        first_op, _ = cigar_tuples[0]
        last_op, _ = cigar_tuples[-1]

        # 检查第一个和最后一个操作符是否是 S 或 H
        return first_op in (CIGAR_S, CIGAR_H) or \
            last_op in (CIGAR_S, CIGAR_H)

    bamfile = pysam.AlignmentFile(fname, "rb")
    for read in bamfile:
        if has_overhang(read.cigartuples):
            print(f"Read {read.query_name} has an overhang.")

def parird_reads_not_in_same_position(fname):
    mmm = {}
    for line in open(fname):
        if line.startswith("@"): continue
        its = line.split()
        _, n, _, p, _ = its[0].split('_')
        if n not in mmm:
            mmm[n] = int(its[3])
        else:
            if mmm[n] != int(its[3]):
                print(line, end="")

def collect_reads_at_position(fname, pos, offset):
    for line in open(fname):
        if line.startswith("@"): continue

        its = line.split()
        if abs(pos-int(its[3])) < offset:
            print(line, end="")


def test():
    import numpy as np  
    sss = "86,90,129,88,72,89,29,138,135,146,131,124,132"
    sss = np.array([int(x) for x in sss.split(",")])

    def detect_iqr(data):
        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        print(q1, q3, iqr)
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        outliers = [i for i, x in enumerate(data) if x < lower_bound or x > upper_bound]
        return outliers, lower_bound, upper_bound

    print(detect_iqr(sss))

if __name__ == "__main__":
    #parird_reads_not_in_same_position(sys.argv[1])
    
    #test()
    #collect_reads_with_overhang(sys.argv[1])

    if sys.argv[1] == "location":
        collect_reads_at_position(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
    elif sys.argv[1] == "overhang":
        collect_reads_with_overhang(sys.argv[2])
