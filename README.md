fast2snp
========

fast2snp is a fast and accurate reference-free tool for detecting heterozygous small variants (SNPs and indels) directly from sequencing reads. It constructs a de Bruijn graph from heterozygous k-mers and performs coverage-aware localized graph traversal to identify variant-containing bubbles, avoiding global graph traversal.

Features
--------

Reference-free: detects variants directly from reads without a reference genome.
Heterozygous SNP and indel detection: identifies both variant types in a single run.
Broad platform compatibility: works with Illumina short reads, PacBio HiFi, and ONT R10 long reads.
Fast: achieves 2.3-16.6x speedup in CPU time compared to DiscoSnp++ across diverse datasets.
Robust: maintains stable performance across varying sequencing depths and genome complexities.

Dependencies
------------

C++17 compiler (e.g., g++)
OpenMP
KMC (https://github.com/refresh-bio/KMC)
parallel-hashmap (https://github.com/greg7mdp/parallel-hashmap)
Python 3
SciPy, NumPy, pandas, pysam

Installation
------------

git clone --recurse-submodules https://github.com/lemene/fast2snp.git
cd fast2snp
make

The compiled binary will be placed in ../build/bin/fast2snp. Add it to your PATH or call it with the full path.

Quick Start
-----------

Step 1: Count k-mers with KMC

kmc -k31 -m64 -t8 reads.fastq.gz kmc_db tmp_dir

Step 2: Estimate heterozygous k-mer range

python scripts/find_het_kmers.py kmc_db.histo -k 31 -o output_dir

This produces a plot and prints the estimated lowCov and highCov values.

Step 3: Run fast2snp

fast2snp -i kmc_db -l <lowCov> -r <highCov> -k 31 -t 8 -o output_dir

Alternatively, use the wrapper script that performs all steps automatically:

python run_fast2snp.py -i reads.fastq.gz -k 31 -t 8 -o output_dir

Input
-----

Reads: FASTQ file(s) (single-end or paired-end) or a KMC database prefix.
k-mer size: default 31.
Heterozygous coverage range: estimated from the k-mer spectrum or provided manually.

Output
------

A FASTA file containing paired allele sequences for each detected variant.
Each record contains two sequences: the reference allele and the alternative allele.

Example output:

>Bubble_1_SNP_Path1_Ext8
ACGTACGTACGT...
>Bubble_1_SNP_Path2_Ext8
ACGTTCGTACGT...

Parameters
----------

-i  KMC database prefix (or - for stdin)
-l  Lower bound of heterozygous coverage
-r  Upper bound of heterozygous coverage
-k  k-mer size (default: 31)
-t  Number of threads (default: 8)
-o  Output directory (default: current directory)
-e  Extension length (default: -1, auto)

Evaluation
----------

Scripts for benchmarking and evaluation are provided in the scripts/ directory:

extract_variants_from_bam.py: extracts variants from aligned BAM files for VCF generation.
eval_with_happy.py: runs hap.py to benchmark fast2snp against truth sets.
extract_vcf_from_mpileup.py: extracts candidate variants from mpileup output.

See the individual script headers for usage details.

License
-------

This project is released under the MIT License. See LICENSE for details.

Contact
-------

For questions or issues, please open an issue on GitHub or contact the corresponding author.
