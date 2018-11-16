from week2 import hamming_ball, ApproximatePatternCount

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


if __name__=="__main__":
    print(MotifEnumeration(
        [
            'ACGT',
            'ACGT',
            'ACGT'

        ],
        3, 0
    ))