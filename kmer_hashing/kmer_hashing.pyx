# cython: language_level=3

from libc.stdint cimport int64_t
import cython

@cython.boundscheck(False)  # Disable bounds-checking for the entire function
@cython.wraparound(False)  # Disable negative index wrapping for the entire function
cpdef int64_t hash_kmer(bytes kmer, int64_t array_size):
    cdef int64_t hash_integer = 0
    cdef int i
    cdef int length = len(kmer)
    
    for i in range(length):
        hash_integer = hash_integer * 37 + kmer[i]
    
    return hash_integer % array_size
