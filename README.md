
# fast2snp

**fast2snp** is a fast and accurate reference-free tool for detecting heterozygous small variants (SNPs and indels) directly from sequencing reads. It constructs a coverage-filtered de Bruijn graph and performs heterozygosity-guided localized graph traversal to identify variant-containing bubbles, avoiding exhaustive global graph traversal.

## Features

- **Reference-free**: Detects variants directly from reads without requiring a reference genome.
- **Heterozygous SNP and Indel Detection**: Identifies both variant types in a single run.
- **Broad Platform Compatibility**: Evaluated on Illumina short reads, PacBio HiFi, and ONT R10 long reads.
- **Fast and Efficient**: Employs a localized processing strategy with aggressive in-memory indexing, significantly reducing CPU time while maintaining a highly efficient workflow (at the cost of moderately higher peak memory usage).
- **Robust**: Maintains stable performance across varying sequencing coverages (optimized for $\ge$ 10x) and diverse genome complexities.

## Dependencies

- C++17 compiler (e.g., `g++`)
- OpenMP
- [KMC](https://github.com/refresh-bio/KMC) (k-mer counting)
- [parallel-hashmap](https://github.com/greg7mdp/parallel-hashmap) (included as a git submodule)
- Python 3 (for wrapper scripts and evaluation)
- Python packages: `scipy`, `pysam`, `pandas`, `matplotlib`

## Installation

```bash
# Clone the repository along with the parallel-hashmap submodule
git clone --recurse-submodules https://github.com/lemene/fast2snp.git
cd fast2snp

# Compile the C++ core
make
```
All scripts and compiled binaries will be placed in the `./build/bin/` directory.

## Quick Start

fast2snp can be run in a single step directly from raw sequencing reads. The provided wrapper script `run_fast2snp.py` automatically performs k-mer counting, heterozygous range estimation, and localized graph traversal to extract candidate variant bubbles.

```bash
python run_fast2snp.py <reads.fastq.gz> -o ./out_dir -t 24
```
The final candidate variant bubbles can be found in `variant_calls/bubbles.snp`.

## Step-by-step execution

If you prefer to run the pipeline manually to inspect intermediate files or customize specific steps, follow the instructions below:

**Step 1: K-mer counting**
```bash
mkdir wrkdir
kmc -cs65535 -k31 -t8 reads.fastq.gz kmc_db .
```
This step uses KMC to count the k-mers and generates the KMC database (`kmc_db`).

**Step 2: Heterozygous k-mer identification**
```bash
kmc_tools transform kmc_db histogram kmc_db.histo
find_het_kmers.py kmc_db.histo -k31 -o .
```
This step generates a histogram and runs the estimation script. It will produce a plot named like `kmc_db_xx_yy.png`, where `xx` and `yy` represent the estimated lower and upper boundaries of the heterozygous k-mer range, respectively.

**Step 3: Variant calling**
```bash
kmc_dump -ci <xx> kmc_db /dev/stdout | fast2snp -i - -l <xx> -r <yy> -t 8 -o .
```
In this step, the k-mers are streamed directly from the KMC database to `fast2snp` via a pipe. You can optionally specify advanced parameters such as `-N`, `-L`, and `-S` in this command. The final output (paired bubble sequences) will be saved in `bubbles.snp`.


## Output

fast2snp outputs a FASTA file located at `bubbles.snp`, containing paired allele sequences for each detected variant. Each bubble consists of two sequences representing the reference and alternative alleles.

Example output:
```fasta
>Bubble_1_SNP_Path1_Ext30
ACGTACGTACGT...
>Bubble_1_SNP_Path2_Ext30
ACGTTCGTACGT...
```

## Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-k` | k-mer size | 31 |
| `-N` | Number of top paths retained per group during branch pruning | 20 |
| `-L` | Maximum additional extension length (bounded for small variants) | 50 |
| `-S` | Multiplier for standard deviation when computing the coverage fluctuation threshold ($\tau = \text{mean} + S \times \text{SD}$) | 3 |
| `-t` | Number of threads | 8 |


## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

## Contact

For questions or issues, please open an issue on GitHub or contact the author (lemene@sina.com).