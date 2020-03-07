def Count(Motifs):
    count = {}
    t = len(Motifs)
    j = len(Motifs[0])
    for s in 'ACGT':
        count[s] = []
        for l in range(j):
            count[s].append(0)
    for a in range(t):
        for b in range(j):
            count[Motifs[a][b]][b] += 1
    return count


def Profile(Motifs):
    t = len(Motifs)
    j = len(Motifs[0])
    profile = {}
    for s in 'ACGT':
        profile[s] = []
        for l in range(j):
            profile[s].append(0)
    for a in range(t):
        for b in range(j):
            profile[Motifs[a][b]][b] += 1 / t
    return profile


def Consensus(Motifs):
    t = len(Motifs)
    j = len(Motifs[0])
    profile = {}
    for s in 'ATGC':
        profile[s] = []
        for l in range(j):
            profile[s].append(0)
    for a in range(t):
        for b in range(j):
            profile[Motifs[a][b]][b] += 1 / t
    consensus = ''
    for b in range(j):
        m = 0
        for s in 'ACGT':
            if profile[s][b] > m:
                m = profile[s][b]
                FrequentSymbol = s
        consensus += FrequentSymbol
    return consensus


def HammingDistance(p, q):
    d = 0
    n = len(p)
    for i in range(n):
        if p[i] != q[i]:
            d += 1
    return d


def Score(Motifs):
    consensus = Consensus(Motifs)
    score = 0
    for i in Motifs:
        score += HammingDistance(i, consensus)
    return score


def Pr(Text, Profile):
    p = 1
    n = len(Text)
    for i in range(n):
        p = p * Profile[Text[i]][i]
    return  p


def ProfileMostProbableKmer(text, k, profile):
    nT = len(text)
    scorelist = []
    for i in range(nT - k + 1):
        scorelist.append(Pr(text[i : i + k], profile))
    j = scorelist.index(max(scorelist))
    return text[j : j + k]


def GreedyMotifSearch(Dna, k, t):
    n = len(Dna[0])
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0 : k])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i : i + k])
        for j in range(1, t):
            P = Profile(Motifs[0 : j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
