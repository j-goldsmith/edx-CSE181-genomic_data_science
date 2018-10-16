import numpy as np

def skew(genome):
    skew = 0
    skew_history = [skew]
    for i in range(len(genome)):
        if genome[i] == 'C':
            skew += -1
        elif genome[i] == 'G':
            skew += 1
        skew_history.append(skew)
    return skew_history


def minimum_skew(genome):
    skew_history = np.array(skew(genome))
    return np.where(skew_history == skew_history.min())[0].tolist()


def HammingDistance(p, q):
    p = np.array(list(p))
    q = np.array(list(q))
    return (p != q).sum()


def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    pattern_length = len(Pattern)

    for i in range(len(Text) - pattern_length + 1):
        if HammingDistance(Pattern, Text[i:i+pattern_length]) <= d:
            positions.append(i)
    return positions


def ApproximatePatternCount(Text, Pattern, d):
    return len(ApproximatePatternMatching(Text,Pattern,d))

def FrequentWordsWithMismatches(Text, k, d):
    '''Code Vomit'''
    words = list(set([Text[i:i+k] for i in range(len(Text) - k + 1)]))
    jimmied_words = []
    chars = ['A','T','C','G']

    for word in words:
        for i in range(len(word)):
            for c in chars:
                w_list = list(word)
                w_list[i] = c
                jimmied_words.append("".join(w_list))

    jimmied_words = list(set(jimmied_words))
    approximate_counts = []
    for word in jimmied_words:
        approximate_counts.append(ApproximatePatternCount(Text, word, d))

    approximate_counts = np.array(approximate_counts)
    max_indices = np.where(approximate_counts == approximate_counts.max())[0].tolist()

    return [jimmied_words[i] for i in max_indices]

FrequentWordsWithMismatches('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1)
FrequentWordsWithMismatches('AAAAAAAAAA', 2, 1)