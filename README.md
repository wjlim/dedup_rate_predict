# K-mer Analysis Package 🧬

This package provides powerful tools for analyzing k-mers in FASTQ files with non-overlapping method, optimized for genomics data like the human genome. Leveraging Cython, it handles large datasets efficiently, making it an essential tool for genomic researchers.

## 🌟 Features

- **Efficient Processing**: Analyzes k-mer frequency in gzipped FASTQ files.
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
git clone https://github.com/wjlim/kmer_analysis.git
cd kmer_analysis
pip install .
```

## 🚀 Usage
The package comes with a command-line tool kmer_processor for easy use:
kmer_processor -f /path/to/your/forward_reads.fastq.gz -k 21 -o /path/to/output_prefix

## 📊 Default Configuration Performance
With the default configuration (block_size of 100MB and chunk_size of 300,000), here's what you can expect:

Raw data: 50X Whole Genome Sequencing data
Memory Usage: Approximately 41GB
Processors Used: 16
Time Taken: About ~3000 seconds
This benchmark is based on processing the forward reads raw data of a human whole genome at 50X coverage on a system with adequate memory and CPU resources.

## ⚙️ Customization
Enhance performance by customizing block and chunk sizes with -b and -c flags:
kmer_processor --help

## 📄 License
This project is under the MIT License. For more details, see the LICENSE file.

Happy k-mer analyzing! 🧪