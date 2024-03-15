#!/usr/bin/env python
import multiprocessing
from multiprocessing.pool import Pool
import argparse
import signal
import gzip
import time
import psutil
import sys
import kmer_hashing

class KmerProcessor:
    def __init__(self, fastq_file_forward, k, genome_size, output_prefix, array_size=10**9, threads=8, block_size=100*1024*1024, chunk_size=300000):
        self.fastq_file_forward = fastq_file_forward
        self.k = k
        self.genome_size = genome_size
        self.output_prefix = output_prefix
        self.array_size = array_size
        self.threads = threads
        self.block_size = block_size
        self.chunk_size = chunk_size

    def sequence_generator(self, file_path):
        with gzip.open(file_path, 'rb') as f:
            eof = False
            while not eof:
                lines = f.readlines(self.block_size)
                if not lines:
                    eof = True
                sequences = (line.strip().decode('UTF-8') for i, line in enumerate(lines) if (i % 4) == 1)
                for sequence in sequences:
                    yield sequence

    def process_sequence(self, sequence):
        local_dict = {}
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i+self.k]
            index = kmer_hashing.hash_kmer(kmer.encode(), self.array_size)
            local_dict[kmer] = local_dict.get(kmer, 0)  + 1
        return local_dict

    def process_chunk(self, sequences):
        chunk_dict = {}
        for sequence in sequences:
            sequence_dict = self.process_sequence(sequence)
            for kmer, count in sequence_dict.items():
                chunk_dict[kmer] = chunk_dict.get(kmer, 0) + count
        return chunk_dict

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
            total_dict = {}
            for chunk_dict in pool.imap_unordered(self.process_chunk, self.chunked_sequences_generator()):
                for kmer in chunk_dict:
                    total_dict[kmer] = total_dict.get(kmer, 0) + chunk_dict[kmer]

        elapsed_time = time.time() - start_time
        memory_usage = psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024) - start_memory

        unique_kmers = len(total_dict)
        total_kmers = sum(total_dict.values())
        duplicate_rate = (total_kmers - unique_kmers) / total_kmers if total_kmers else 0
        mean_depth = (total_kmers * 2 / self.genome_size) if self.genome_size else 0

        with open(f"{self.output_prefix}_stats.csv", "w") as f:
            header = "forward_fastq,k,genome_size,duplicate_rate,deduplicate_rate,\
predicted_mean_depth,elapsed_time_seconds,memory_usage_mb"
            body = f"{self.fastq_file_forward},{self.k},{self.genome_size},{duplicate_rate:.2%},\
{1-duplicate_rate:.2%},{mean_depth:.2f},{elapsed_time},{memory_usage}"
            f.write(header + '\n' + body + '\n')
            
def main():
    parser = argparse.ArgumentParser(description="Analyze FASTQ files to achieve a target coverage of the human genome, given the total number of reads, and calculate PCR duplicates and mean depth.")
    parser.add_argument('-f', '--forward', required=True, help='Input forward gzip compressed FASTQ file')
    parser.add_argument('-k', '--kmer_size', type=int, default=18, help='Size of k-mers (default: 18)')
    parser.add_argument('-g', '--genome_size', type=int, default=3000000000, help='Size of the genome in base pairs (default: 3Gb)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for statistics')
    parser.add_argument('-p', '--threads', type=int, default=8, help='Number of threads for parallel processing')
    parser.add_argument('-b', '--block_size', type=int, default=100*1024*1024, help='Block size for reading FASTQ files (default: 100MB)')
    parser.add_argument('-c', '--chunk_size', type=int, default=300000, help='Chunk size for processing sequences (default: 300000)')
    args = parser.parse_args()

    kmer_processor = KmerProcessor(args.forward, args.kmer_size, args.genome_size,
                                args.output, array_size=10**9, threads=args.threads,
                                block_size=args.block_size, chunk_size=args.chunk_size)
    kmer_processor.process_file()
    
if __name__ == "__main__":
    main()
