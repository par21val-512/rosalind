from problems_011_020 import codon_table
from problems_001_010 import revc
import numpy as np
import math


def lscm(s: str):
    strs = []
    for dna in s.split('>')[1:]:
        strs.append(''.join(dna.split('\n')[1:]))
    strs = sorted(strs, key=len)

    # Step 1: Get shortest string in set
    short = strs[0]
    others = strs[1:]

    # Step 2: get all substrings of the shortest string, and sort in descending length order
    substrings = set()
    for j in range(len(short)+1):
        for i in range(j):
            substrings.add(short[i:j])
    substrings = sorted(substrings, key=len, reverse=True)

    # Step 3: return first substring that is a substring of all others
    for sub in substrings:
        if all(sub in string for string in others):
            return sub

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


def pmch(s: str):
    dna = ''.join(s.split('\n')[1:])
    au_pairs = 0
    gc_pairs = 0
    for i in range(len(dna)):
        if dna[i] == 'A':
            au_pairs += 1
        if dna[i] == 'C':
            gc_pairs += 1
    return math.factorial(au_pairs) * math.factorial(gc_pairs)

def pper(n: int, k: int):
    res = 1
    for i in range(n, 0, -1):
        res = res * i
    for j in range(n-k, 0, -1):
        res = res // j
    return res % 1000000

def tree(s: str):
    n, edges = int(s.split('\n')[0]), s.split('\n')[1:]
    return (n-1) - len(edges)

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

def tran(s: str):
    s1 = ''.join(s.split('>')[1].split('\n')[1:])
    s2 = ''.join(s.split('>')[2].split('\n')[1:])
    transitions = 0
    transversions = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            continue
        if sorted([s1[i], s2[i]]) == ['A', 'G'] or sorted([s1[i], s2[i]]) == ['C', 'T']:
            transitions += 1
        else:
            transversions += 1
    return (transitions / transversions)


def lexf(symbols: list, n: int):  # input symbols as list from string
    symbols = sorted(symbols)
    if n == 1:
        return symbols
    res = []
    for symbol in symbols:
        others = lexf(symbols, n-1)
        for s in others:
            res.append(symbol + s)
    return res

def sign(n : int):
    nums = [i for i in range(1, n+1, 1)]  # generate list of [1, ..., n] first

    def neg(x: int):
        return -1 * x

    def sign_helper(nums : list):
        if len(nums) == 1:
            return [nums[:], list(map(neg, nums[:]))]
        res = []
        for _ in range(2 * len(nums)):
            n = nums.pop(0)
            perms = sign_helper(nums)
            for p in perms:
                p.append(n)
            res.extend(perms)
            nums.append(-n)
        return res
    return sign_helper(nums)

if __name__ == '__main__':
    file_path = 'datasets/rosalind_lcsm.txt'
    # read in .txt file when input is one string
    with open(file_path, 'r') as file:
        input = file.read().strip()
    print(lscm(input))
