def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i : i + len(Pattern)] == Pattern:
            count += 1
    return count


def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0 : n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i : i+n//2], symbol)
    return array


def SymbolArray(Genome, symbol):
    array = []
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0 : n//2]
    for i in range(n):
        count = PatternCount(ExtendedGenome[i : i + n//2], symbol)
        array.append(count)
    return array


def FastSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0 : n//2]
    array[0] = PatternCount(Genome[0 : n//2], symbol)
    for i in range(1, n):
        array[i] = array[i - 1]
        if ExtendedGenome[i - 1] == symbol:
            array[i] -= 1
        elif ExtendedGenome[i + (n//2) - 1] == symbol:
            array[i] += 1
    return array


def SkewArray(Genome):
    skew = {}
    skew[0] = 0
    n = len(Genome)
    for i in range(n):
        skew[i + 1] = skew[i]
        if Genome[i] == 'C':
            skew[i + 1] -= 1
        elif Genome[i] == 'G':
            skew[i + 1] += 1
    return skew


def MinimumSkew(Genome):
    position = []
    skew = SkewArray(Genome)
    mini = min(skew.values())
    for key, value in skew.items():
        if value == mini:
            position.append(key)
    return position


def HammingDistance(p, q):
    d = 0
    n = len(p)
    for i in range(n):
        if p[i] != q[i]:
            d += 1
    return d


def ApproximatePatternMatching(Text, Genome, d):
    position = []
    nG = len(Genome)
    nT = len(Text)
    for i in range(nG - nT + 1):
        if HammingDistance(Text, Genome[i : i + nT]) <= d:
            position.append(i)
    return position


def ApproximatePatternCount(Text, Genome, d):
    count = 0
    nG = len(Genome)
    nT = len(Text)
    for i in range(nG - nT + 1):
        if HammingDistance(Text, Genome[i : i + nT]) <= d:
            count += 1
    return count