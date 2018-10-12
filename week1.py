import re
from collections import Counter


def PatternCount(Text, Pattern):
    '''
    Count the number of times a string appears as a 
    substring in a longer text. Overlaps allowed
    '''
    return len(re.findall('(?='+Pattern+')', Text))


def FrequentWords(Text, k, t=None):
    '''Find the most frequent k-mers in a string.'''

    words = [Text[i:i+k] for i in range(len(Text) - k + 1)]
    counts = Counter(words)
    
    # if t is not set, use the most frequent count
    if t is None:
        t = counts[max(counts, key=counts.get)]

    return [w for w, c in counts.items() if c == t]


def ReverseComplement(Pattern):
    '''
    Given a nucleotide p, we denote its complementary nucleotide as p*. 
    The reverse complement of a string Pattern = p1 … pn is the string 
    Pattern_rc = pn* … p1* formed by taking the complement of each nucleotide in Pattern, 
    then reversing the resulting string.
    '''
    complement_map = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    complement = [complement_map[x] for x in Pattern]
    return "".join(complement[::-1])


def PatternMatching(Pattern, Genome):
    '''Find all starting positions where Pattern appears as a substring of Genome.'''
    return [m.start(0) for m in re.finditer('(?=('+Pattern+'))', Genome)]


def ClumpFinding(genome, k, L, t):
    '''
    Find patterns forming clumps in a string.
    Input: 
        genome: string  
        k: int, k-mer size 
        L: int, window size for clump counting
        t: int, target clump count 
    Output: All distinct k-mers forming (L, t)-clumps in Genome.
    '''
    clumps = []
    for i in range(len(genome) - L + 1):
        clumps.extend(FrequentWords(genome[i:i+L], k, t))
    return set(clumps)
