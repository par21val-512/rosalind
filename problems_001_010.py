# Rosalind Bioinformatics Stronghold problems

def dna(s: str):
    res = {}
    for i in range(len(s)):
        res[s[i]] = res.get(s[i], 0) + 1
    return sorted(res.items())


def rna(s: str):
    res = ''
    for i in range(len(s)):
        if s[i] == 'T':
            res += 'U'
        else:
            res += s[i]
    return res

def revc(s: str):
    rev = s[::-1].upper()
    sol = ''
    complement_dict = {'A' : 'T',
                       'C' : 'G',
                       'G' : 'C',
                       'T' : 'A'}
    for i in range(len(rev)):
        if rev[i] == '\n':
            pass
        sol += complement_dict[rev[i]]
    return sol


def iprb(k: int, m: int, n: int):
    total = (k + m + n) * ((k + m + n) - 1)  # total number of permutations, double count later !
    # homozygous dominant - AA ; heterozygous - Aa ; homozygous recessive - aa
    kk = k * (k-1)  # AA x AA => 100% AA
    km = 2 * k * m  # AA x Aa => 50% AA + 50% Aa
    kn = 2 * k * n  # AA x aa => 100% Aa
    mm = m * (m-1) * 0.75  # Aa x Aa => 25% AA + 50% Aa
    mn = 2 * m * n * 0.5  # Aa x aa => 50% Aa
    dom_pairs = kk + km + kn + mm + mn
    return dom_pairs / total


def fib(n: int, k: int):
    dp = [0] * (n+1)
    dp[0] = 0
    dp[1] = 1
    for i in range(2, n+1, 1):
        dp[i] = dp[i-1] + k * dp[i-2]
    return dp[n]


def gc(s: str):
    # converts FATSA formatted file to dictionary
    dic = {}
    strings = s.split('>')[1:]
    for s in strings:
        splt = s.split('\n')
        dic[''.join(splt[1:])] = splt[0]

    res = {}
    for s in list(dict.keys()):
        count = 0
        for i in range(len(s)):
            if s[i] == 'C' or s[i] == 'G':
                count += 1
        content = 100 * count / len(s)
        res[content] = dic[s]
    return res[max(res)], max(res)


def prot(s: str):
    codon_table = {
        "UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
        "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
        "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
        "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
        "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
        "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
        "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
        "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
        "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
        "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
        "UAA": "Stop", "CAA": "Q", "AAA": "K", "GAA": "E",
        "UAG": "Stop", "CAG": "Q", "AAG": "K", "GAG": "E",
        "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
        "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
        "UGA": "Stop", "CGA": "R", "AGA": "R", "GGA": "G",
        "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
    }

    res = ''
    for i in range(0, len(s), 3):
        amino_acid = codon_table[s[i:i+3]]
        if amino_acid == "Stop":
            break
        res += amino_acid
    return res


def subs(f: str):
    splt = f.split('\n')
    s, t = splt[0], splt[1]
    res = []
    for i in range(len(s) - len(t)):
        if s[i:i+len(t)] == t:
            res.append(i+1)
    return res


def hamm(f: str):
    splt = f.split('\n')
    s, t = splt[0], splt[1]
    dist = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            dist += 1
    return dist


def iev(a1: int, a2: int, a3: int, a4: int, a5: int, a6: int):
    e1 = a1 * 2
    e2 = a2 * 2
    e3 = a3 * 2
    e4 = a4 * 0.75 * 2
    e5 = a5 * 0.5 * 2
    e6 = a6 * 0
    return e1 + e2 + e3 + e4 + e5 + e6


if __name__ == '__main__':
    file_path = 'datasets/rosalind_prot.txt'
    # read in .txt file when input is one string
    with open(file_path, 'r') as file:
        input = file.read().strip()
    print(prot(input))