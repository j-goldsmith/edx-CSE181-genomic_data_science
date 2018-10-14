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
    return (p!=q).sum()

def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    pattern_length = len(Pattern)

    for i in range(len(Text) - pattern_length + 1):
        if HammingDistance(Pattern, Text[i:i+pattern_length]) <= d:
            positions.append(i)
    return positions
