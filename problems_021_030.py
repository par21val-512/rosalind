from problems_011_020 import codon_table
from problems_001_010 import revc
import numpy as np


def revp(s: str):
    dna = ''.join(s.split('\n')[1:])
    res = []
    for i in range(len(dna)):
        for j in range(4, 13, 1):
            if i+j>len(dna):
                continue
            substr = dna[i:i+j]
            if substr == revc(substr):
                res.append((i+1, j))
    return res


def pper(n: int, k: int):
    res = 1
    for i in range(n, 0, -1):
        res = res * i
    for j in range(n-k, 0, -1):
        res = res // j
    return res % 1000000


def sseq(l: str):  # two pointers approach
    s = ''.join(l.split('>')[1].split('\n')[1:])
    t = ''.join(l.split('>')[2].split('\n')[1:])
    print(s)
    print(t)
    res = []
    i = j = 0
    while i < len(s) and j < len(t):
        if s[i] == t[j]:
            res.append(i+1)
            j += 1
        i += 1
    return res


def prob(s: str):
    dna, A = s.split('\n')[0], list(map(float, ''.join(s.split('\n')[1:]).split(' ')))

    def content_prob(c: str, x: float):
        if c == 'G' or c == 'C':
            return np.log10(x / 2)
        else:  # c == 'A' or c == 'T'
            return np.log10((1-x) / 2)
    res = []
    for i in range(len(A)):
        log_prob = 0
        for j in range(len(dna)):
            log_prob += content_prob(dna[j], A[i])
        res.append(log_prob)
    return res

if __name__ == '__main__':
    file_path = 'datasets/rosalind_prob.txt'
    # read in .txt file when input is one string
    with open(file_path, 'r') as file:
        input = file.read().strip()
    print(*prob(input))