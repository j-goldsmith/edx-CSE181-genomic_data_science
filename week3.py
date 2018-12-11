from week2 import hamming_ball, ApproximatePatternCount, HammingDistance
import sys
from itertools import product

def MotifEnumeration(Dna, k, d):
    '''
    Patterns ← an empty set
    for each k-mer Pattern in Dna
        for each k-mer Pattern’ differing from Pattern by at most d mismatches
            if Pattern' appears in each string from Dna with at most d mismatches
                add Pattern' to Patterns
    remove duplicates from Patterns
    return Patterns
    '''
    patterns = []
    kmer_prime_counts = {}
    for sample in Dna:
        for i in range(len(sample) - k + 1):
            kmer = sample[i:i+k]
            ball = hamming_ball(kmer, d, 'ATCG')
            for kmer_prime in ball:
                if len(kmer_prime) < k:
                    continue
                presence_checks = []
                for check_sample in Dna:
                    count = ApproximatePatternCount(check_sample, kmer_prime, d)
                    if count > 0 :
                        presence_checks.append(1)
                if len(presence_checks) == len(Dna):
                    patterns.append(kmer_prime)
    patterns = sorted(set(patterns))
    return patterns

def get_potential_kmers(k):
    dna = 'A', 'C', 'G', 'T'
    return [''.join(i) for i in product(dna, repeat = k)]


def distance_score(dna, kmer):
    score = 0
    k = len(kmer)
    for sample in dna:
        sample_score = sys.maxsize
        for i in range(len(sample)-k+1):
            d = HammingDistance(kmer, sample[i:i+k])
            if (d < sample_score):
                sample_score = d
        score += sample_score
    return score

def MedianString(dna, k):
    potential_medians = get_potential_kmers(k)
    distance = sys.maxsize
    median = None
    for kmer in potential_medians:
       score = distance_score(dna, kmer)
       if (distance > score):
           distance = score
           median = kmer
    return median

if __name__=="__main__":
    print(MedianString(
    [
        'AAATTGACGCAT',
        'GACGACCACGTT',
        'CGTCAGCGCCTG',
        'GCTGAGCACCGG',
        'AGTACGGGACAG'
    ], 
    3))
''' 
   print(MotifEnumeration(
        [
            'ACGT',
            'ACGT',
            'ACGT'
        ],
        3, 0
    ))
'''

    