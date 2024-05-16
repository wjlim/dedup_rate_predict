# Deduplication Rate Prediction v.1.0 Package 🧬

This package provides powerful tools for predicting deduplication rates in FASTQ files, optimized for genomics data like the human genome. Leveraging Cython, it handles large datasets efficiently, making it an essential tool for genomic researchers.

## 🌟 Features

- **Efficient Processing**: Predicts deduplication rates in gzipped FASTQ files.
- **Multiprocessing Support**: Utilizes multiprocessing for enhanced performance on multi-core systems.
- **Customizability**: Offers adjustable block and chunk sizes for sequence reading and processing.
- **Optimized with Cython**: Includes critical computational path optimizations.

## 📋 Requirements

- Python 3.6 or higher
- Cython
- psutil

## 🛠 Installation

### Clone and Install

```bash
git clone https://github.com/wjlim/dedup_rate_prediction.git
cd kmer_analysis
pip install .
```
## 🚀 Usage
The package comes with a command-line tool dedup_predictor for easy use:
```bash
kmer_processor -f /path/to/your/forward_reads.fastq.gz -o /path/to/output_prefix
```
## 📊 Default Configuration Performance
With the default configuration (block_size of 100MB and chunk_size of 300,000), here's what you can expect:

Raw data: 50X Whole Genome Sequencing data
Memory Usage: Approximately 41GB
Processors Used: 16
Time Taken: About ~3000 seconds

This benchmark is based on processing the forward reads raw data of a human whole genome at 50X coverage on a system with adequate memory and CPU resources.

⚙️ Customization
Enhance performance by customizing block, chunk sizes and max block with -b ,-c and -m flags:

```bash
kmer_processor --help
usage: kmer_processor [-h] -f IFILE -o OUTPUT [-p THREADS] [-c CHUNK_SIZE] [-a ARRAY_SIZE] [-b BLOCK_SIZE] [-m MAX_BLOCK_SIZE] [-u CUTOFF]

Analyze FASTQ files to calculate PCR duplicates and mean depth.

options:
  -h, --help            show this help message and exit
  -f IFILE, --ifile IFILE
                        Input gzip compressed FASTQ file
  -o OUTPUT, --output OUTPUT
                        Output prefix for files
  -p THREADS, --threads THREADS
                        Number of threads for parallel processing (default: 8)
  -c CHUNK_SIZE, --chunk_size CHUNK_SIZE
                        Chunk size for processing sequences
  -a ARRAY_SIZE, --array_size ARRAY_SIZE
                        Array size for hashing kmer (default: 10e9)
  -b BLOCK_SIZE, --block_size BLOCK_SIZE
                        Block size for reading FASTQ files (default: 500MB)
  -m MAX_BLOCK_SIZE, --max_block_size MAX_BLOCK_SIZE
                        Limitation of read file size (default: 20Gb)
  -u CUTOFF, --cutoff CUTOFF
                        cutoff for determining PCR duplicates (default: 2)
```

📄 License
This project is under the MIT License. For more details, see the LICENSE file.
