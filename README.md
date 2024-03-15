# K-mer Analysis Package

This package provides tools for analyzing k-mers in FASTQ files, particularly useful for genomics data like the human genome. It includes optimizations with Cython to handle large datasets efficiently.

## Features

- Processing of gzipped FASTQ files to analyze k-mer frequency.
- Use of multiprocessing for improved performance on multi-core systems.
- Customizable block and chunk sizes for reading and processing sequences.
- Cython optimizations for critical computational paths.

## Requirements

- Python 3.6 or higher
- Cython
- psutil

## Installation

Clone this repository and navigate into the package directory:

```bash
git clone https://github.com/wjlim/kmer_analysis.git
cd kmer_analysis

Install the package:
pip install .

Usage
The package includes a command-line tool kmer_processor which can be used as follows:
kmer_processor -f /path/to/your/forward_reads.fastq.gz -k 21 -o /path/to/output_prefix

Default Configuration Performance
Using the default configuration (block_size of 100MB and chunk_size of 300,000), processing the forward reads raw data of a human whole genome at 30X coverage:

Memory Usage: Approximately 80GB
Processors Used: 8
Time Taken: About 30 minutes
These metrics are based on processing data on a system configured with sufficient memory and CPU resources to handle high-throughput genomic data.

Customization
You can customize the block and chunk sizes used for reading and processing the FASTQ files by using the -b and -c flags, respectively. This can be useful for adjusting the performance based on your system's hardware capabilities or specific requirements of your dataset.

For more details on the available options, use:
kmer_processor --help

License
This project is licensed under the MIT License - see the LICENSE file for details.
