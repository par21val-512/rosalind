import math


def mmch(s: str):
    counts = {}
    for i in range(len(s)):
        counts[s[i]] = counts.get(s[i], 0) + 1
    au_matchings = math.perm(max(counts['A'], counts['U']), min(counts['A'], counts['U']))
    gc_matchings = math.perm(max(counts['G'], counts['C']), min(counts['G'], counts['C']))
    return au_matchings * gc_matchings

def inod(n: int):
    return n-2

def sset(n: int):
    return (2 ** n) % 1000000



def kmp(s: str):
    dna = ''.join(s.split('\n')[1:]) # preprocess string to just get dna bit
    m = -1
    res = [m]
    for i in range(len(dna)):
        while m >= 0 and dna[m] != dna[i]:
            m = res[m]
        m += 1
        res.append(m)
    res.pop(0)
    return res



if __name__ == '__main__':
    file_path = 'datasets/rosalind_kmp.txt'
    # read in .txt file when input is one string
    with open(file_path, 'r') as file:
        input = file.read().strip()
    print(*kmp('''>Rosalind_87
CAGCATGGTATCACAGCAGAG'''))