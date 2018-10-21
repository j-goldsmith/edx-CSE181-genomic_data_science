import numpy as np
from week1 import ReverseComplement
from itertools import chain, combinations, product

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


def max_skew(genome):
    skew_history = np.array(skew(genome))
    return np.where(skew_history == skew_history.max())[0].tolist()


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

def hamming_circle(s, n, alphabet):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.

    >>> sorted(hamming_circle('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_circle('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_circle('aaa', 2, 'ab'))
    ['abb', 'bab', 'bba']

    """
    for positions in combinations(range(len(s)), n):
        for replacements in product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin)

def hamming_ball(s, n, alphabet):
    """Generate strings over alphabet whose Hamming distance from s is
    less than or equal to n.

    >>> sorted(hamming_ball('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_ball('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'abc', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_ball('aaa', 2, 'ab'))
    ['aaa', 'aab', 'aba', 'abb', 'baa', 'bab', 'bba']

    """
    return chain.from_iterable(hamming_circle(s, i, alphabet)
                               for i in range(n + 1))
def FrequentWordsWithMismatchesAndReverseComplements(text, k, d):
    rc_text = ReverseComplement(text)
    words = [text[i:i+k] for i in range(len(text) - k + 1)]
    rc_words = [rc_text[i:i+k] for i in range(len(rc_text) - k + 1)]
    words.extend(rc_words)
    words = list(set(words))

    jimmied_words = []
    chars = ['A','T','C','G']
    for word in words:
        jimmied_words = hamming_ball(word,d,'ATCG')

    jimmied_words = list(set(jimmied_words))
    final_words = []
    for word in jimmied_words:
        final_words.append(word)
        final_words.append(ReverseComplement(word))
    final_words = list(set(final_words))

    approximate_counts = []
    for word in final_words:

        score = ApproximatePatternCount(text, word, d) + ApproximatePatternCount(text,ReverseComplement(word),d)
        approximate_counts.append(score)      

    approximate_counts = np.array(approximate_counts)
    max_indices = np.where(approximate_counts == approximate_counts.max())[0].tolist()
    return [final_words[i] for i in max_indices]