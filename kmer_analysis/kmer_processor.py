import argparse
import gzip
import multiprocessing
import os
import time
import gc
import signal
import matplotlib.pyplot as plt
import psutil
import pickle

from kmer_hashing.kmer_hashing import hash_kmer
from collections import Counter

class KmerProcessor:
    def __init__(self, fastq_file, output_prefix, k, cutoff, array_size=10e9, threads=8, block_size=500*1024*1024, chunk_size=300000):
        self.fastq_file = fastq_file
        self.output_prefix = output_prefix
        self.k = k
        self.cutoff = cutoff
        self.array_size = array_size
        self.threads = threads
        self.block_size = block_size
        self.chunk_size = chunk_size

    def sequence_generator(self, file_path, block_size=500*1024*1024):
        with gzip.open(file_path, 'rt') as f:
            block = ""
            while True:
                new_block = f.read(block_size)
                if not new_block and not block:
                    break
                block += new_block
                lines = block.split('\n')
                pos = block.rfind('\n')
                nlines = len(lines)

                if pos != nlines - 1:
                    remainder = f.readline()
                    last_incomplete_line = lines.pop(-1) + remainder
                    lines.append(last_incomplete_line)

                complete_lines_count = nlines - (nlines % 4)
                
                for i in range(1, complete_lines_count, 4):
                    seq = lines[i].strip()
                    if 'N' not in seq:
                        yield seq
                block = '\n'.join(lines[complete_lines_count:])
                
    def process_sequence(self, sequence, array_size = 10e9):
        local_counter = Counter()
        for i in range(0, len(sequence) - self.k + 1, self.k):
            kmer = sequence[i:i+self.k]
            index = hash_kmer(kmer.encode(), array_size)
            local_counter[index] += 1
        return local_counter
    
    def process_chunk(self, sequences):
        chunk_counter = Counter()
        for sequence in sequences:
            sequence_dict = self.process_sequence(sequence, self.array_size)
            chunk_counter.update(sequence_dict)
        return chunk_counter

    def chunked_sequences_generator(self):
        chunk = []
        for sequence in self.sequence_generator(self.fastq_file, self.block_size):
            chunk.append(sequence)
            if len(chunk) == self.chunk_size:
                yield chunk
                chunk = []
                gc.collect()
        if chunk:
            yield chunk


def plot_frequency(frequency_dict, output_file):
    labels, values = zip(*sorted(frequency_dict.items()))
    plt.figure(figsize=(10, 6))
    plt.plot(labels, values, marker='o', linestyle='-', color='purple')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Count Frequency (log scale)')
    plt.ylabel('Number of Occurrences (log scale)')
    plt.title('Frequency of Count Frequencies')
    plt.grid(True)
    plt.xticks(labels, labels=[str(label) for label in labels])
    plt.savefig(output_file, dpi=300)

def calculate_statistics(total_counter, start_time, initial_memory, args):
    plot_frequency(Counter(total_counter.values()), f'{args.output}_frequency_plot.png')
    cutoff = args.cutoff
    
    total_kmers = sum(total_counter.values())
    uniq_kmers = sum(count if count <= cutoff else cutoff for count in total_counter.values())
    duplicate_rate = (total_kmers - uniq_kmers) / total_kmers

    elapsed_time = time.time() - start_time
    final_memory = psutil.Process(os.getpid()).memory_info().rss / (1024 ** 2)
    memory_usage = final_memory - initial_memory

    header = (
        "Forward FASTQ"
        ",Deduplicated Rate"
        ",Duplicate Rate"
        ",Unique K-mers"
        ",Total K-mers"
        ",Elapsed Time (seconds)"
        ",Memory Usage (MB)"
    )

    data = (
        f"{args.ifile}"
        f",{1-duplicate_rate:.2%}"
        f",{duplicate_rate:.2%}"
        f",{uniq_kmers}"
        f",{total_kmers}"
        f",{elapsed_time:.2f}"
        f",{memory_usage:.2f}"
    )

    with open(f"{args.output}", "w") as file:
        file.write(header + '\n' + data + '\n')
        
    with open(f"{args.output}_kmer_counter.pkl", "wb") as f:
        pickle.dump(total_counter, f)
        
def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze FASTQ files to calculate PCR duplicates and mean depth.")
    parser.add_argument('-f', '--ifile', required=True, help='Input gzip compressed FASTQ file')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for files')
    parser.add_argument('-p', '--threads', type=int, default=16, help='Number of threads for parallel processing (default: 16)')
    parser.add_argument('-c', '--chunk_size', type=int, default=300000, help='Chunk size for processing sequences (default: 300000)')
    parser.add_argument('-a', '--array_size', type=int, default=int(10e9), help='Array size for hashing kmer (default: 10^9)')
    parser.add_argument('-t', '--cutoff', type=int, default=10, help='cutoff for determining PCR duplicates (default: 10)')
    parser.add_argument('-b', '--block_size', type=int, default=500*1024*1024, help='Block size for reading FASTQ files (default: 500Mb)')
    parser.add_argument('-k', '--kmer_size', type=int, default=13, help='kmer size (default: 13)')
    return parser.parse_args()

def main():
    start_time = time.time()
    args = parse_arguments()
    multiprocessing.set_start_method('spawn')
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    
    kmer_processor = KmerProcessor(args.ifile, args.output, args.kmer_size, args.cutoff, array_size=args.array_size, threads=args.threads,
                                    block_size=args.block_size, chunk_size=args.chunk_size)
    initial_memory = psutil.Process(os.getpid()).memory_info().rss / (1024 ** 2)

    with multiprocessing.Pool(processes=args.threads) as pool:
        total_counter = Counter()
        for chunk_dict in pool.imap_unordered(kmer_processor.process_chunk, kmer_processor.chunked_sequences_generator()):
            total_counter.update(chunk_dict)
            
    plot_frequency(total_counter, f'{args.output}.png')
    calculate_statistics(total_counter, start_time, initial_memory, args)

if __name__ == '__main__':
    main()
