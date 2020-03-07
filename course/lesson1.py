def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i: i + len(Pattern)] == Pattern:
            count += 1
    return count


def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n - k + 1):
        Pattern = Text[i:i + k]
        if Pattern not in freq.keys():
            freq[Pattern] = 1
        else:
            freq[Pattern] += 1
    return freq


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


def Reverse(Pattern):
    reverse = ''
    n = len(Pattern)
    for i in range(n):
        reverse += (Pattern[n - i - 1])
    return reverse


def Complementary(Pattern):
    complementary = ''
    n = len(Pattern)
    for i in range(n):
        if Pattern[i] == "T":
            complementary += "A"
        elif Pattern[i] == "A":
            complementary += "T"
        elif Pattern[i] == "G":
            complementary += "C"
        elif Pattern[i] == "C":
            complementary += "G"
    return complementary


def ReverseComplementary(Pattern):
    Pattern = Reverse(Pattern)
    Pattern = Complementary(Pattern)
    return Pattern


def PatternMatching(Pattern, Genome):
    positions = []
    n = len(Genome)
    for i in range(n - len(Pattern) + 1):
        if Genome[i:i + len(Pattern)] == Pattern:
            positions.append(i)
    return positions


def PatternMatching(Pattern, Genome):
    positions = []
    n = len(Genome)
    for i in range(n - len(Pattern) + 1):
        if Genome[i:i + len(Pattern)] == Pattern or ReverseComplementary(Pattern):
            positions.append(i)
    return positions

