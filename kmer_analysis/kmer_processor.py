#!/usr/bin/env python
import multiprocessing
from multiprocessing.pool import Pool
import argparse
import signal
import gzip
import time
import psutil
import sys
import os
from collections import Counter
from kmer_hashing.kmer_hashing import hash_kmer

class KmerProcessor:
    def __init__(self, fastq_file_forward, k, genome_size, output_prefix, array_size=10**9, threads=8, block_size=10*1024*1024, chunk_size=300000):
        self.fastq_file_forward = fastq_file_forward
        self.k = k
        self.genome_size = genome_size
        self.output_prefix = output_prefix
        self.array_size = array_size
        self.threads = threads
        self.block_size = block_size
        self.chunk_size = chunk_size
            
    def sequence_generator(self, file_path, block_size=10*1024*1024):
        with gzip.open(file_path, 'rt') as f:
            eof = False
            while not eof:
                lines = f.readlines(block_size)
                if not lines:
                    eof = True
                sequences = (line.strip() for i, line in enumerate(lines) if (i % 4) == 1)
                for sequence in sequences:
                    yield sequence

    def process_sequence(self, sequence):
        local_counter = Counter()
        for i in range(0, len(sequence) - self.k + 1, self.k):
            kmer = sequence[i:i+self.k]
            index = hash_kmer(kmer.encode(), self.array_size)
            local_counter[index] += 1
        return local_counter

    def process_chunk(self, sequences):
        chunk_counter = Counter()
        for sequence in sequences:
            sequence_dict = self.process_sequence(sequence)
            chunk_counter.update(sequence_dict)
        return chunk_counter

    def chunked_sequences_generator(self):
        chunk = []
        for sequence in self.sequence_generator(self.fastq_file_forward):
            chunk.append(sequence)
            if len(chunk) == self.chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk
            
    def process_file(self):
        start_time = time.time()
        start_memory = psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)

        with Pool(processes=self.threads) as pool:
            total_counter = Counter()
            for chunk_counter in pool.imap_unordered(self.process_chunk, self.chunked_sequences_generator()):
                total_counter.update(chunk_counter)

        elapsed_time = time.time() - start_time
        final_memory = psutil.Process().memory_info().rss / (1024 ** 2)  # Memory in MB
        memory_usage = final_memory - start_memory

        unique_kmers = len(total_counter)
        total_kmers = sum(total_counter.values())
        dup_kmers = (total_kmers - unique_kmers)
        
        duplicate_rate = dup_kmers / total_kmers
        mean_depth = total_kmers * self.k * 2 / self.genome_size if self.genome_size else 0
        mappable_mean_depth = unique_kmers * self.k * 2 / self.genome_size if self.genome_size else 0

        with open(f"{self.output_prefix}_stats.csv", "w") as f:
            header = "forward_fastq,k,genome_size,duplicate_rate,deduplicate_rate,\
predicted_mean_depth,predicted_mappable_mean_depth,elapsed_time_seconds,memory_usage_mb"
            body = f"{self.fastq_file_forward},{self.k},{self.genome_size},{duplicate_rate:.2%},\
{1-duplicate_rate:.2%},{mean_depth:.2f},{mappable_mean_depth:.2f},{elapsed_time},{memory_usage}"
            f.write(header + '\n' + body + '\n')
            
def main():
    multiprocessing.set_start_method('spawn')
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    parser = argparse.ArgumentParser(description="Analyze FASTQ files to achieve a target coverage of the human genome, given the total number of reads, and calculate PCR duplicates and mean depth.")
    parser.add_argument('-f', '--forward', required=True, help='Input forward gzip compressed FASTQ file')
    parser.add_argument('-k', '--kmer_size', type=int, default=18, help='Size of k-mers (default: 18)')
    parser.add_argument('-g', '--genome_size', type=int, default=3000000000, help='Size of the genome in base pairs (default: 3Gb)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for statistics')
    parser.add_argument('-p', '--threads', type=int, default=8, help='Number of threads for parallel processing')
    parser.add_argument('-b', '--block_size', type=int, default=10*1024*1024, help='Block size for reading FASTQ files (default: 100MB)')
    parser.add_argument('-c', '--chunk_size', type=int, default=300000, help='Chunk size for processing sequences (default: 300000)')
    parser.add_argument('-a', '--array_size', type=int, default=10**9, help='array_size for hashing kmer (default: 10**9)')
    args = parser.parse_args()

    kmer_processor = KmerProcessor(args.forward, args.kmer_size, args.genome_size,
                                args.output, array_size=args.array_size, threads=args.threads,
                                block_size=args.block_size, chunk_size=args.chunk_size)
    kmer_processor.process_file()
    
if __name__ == "__main__":
    main()