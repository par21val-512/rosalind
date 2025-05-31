import numpy as np
from scipy.stats import binom
import requests
import re
from problems_001_010 import revc

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

def fibd(n : int, m : int):
    inf = [0] * (n+1) # number of infant rabbits at week i
    mat = [0] * (n+1) # number of mature rabbits at week i
    inf[0], mat[0] = 0, 0
    inf[1], mat[1] = 1, 0
    for i in range(2, n+1, 1):
        if i >= m:
            mat[i] = mat[i-1] + inf[i-1] - inf[i-m]
        else:
            mat[i] = mat[i-1] + inf[i-1]
        inf[i] = mat[i-1]
    return mat[n] + inf[n]

def mrna(s : str):
    amino_acid_counts = {'F': 2, 'L': 6, 'I': 3, 'V': 4,
                         'M': 1, 'S': 6, 'P': 4, 'T': 4,
                         'A': 4, 'Y': 2, 'H': 2, 'N': 2,
                         'D': 2, 'Stop': 3, 'Q': 2, 'K': 2,
                         'E': 2, 'C': 2, 'R': 6, 'G': 4, 'W': 1}
    res = 1
    for i in range(len(s)):
        res *= amino_acid_counts[s[i]]
        if res >= 1000000:
            res = res % 1000000
    return (3 * res) % 1000000

def lia(k : int, N : int):
    num_ppl = 2 ** k # X ~ Bin(num_ppl, 1/4)
    return 1 - binom.cdf(N-1, num_ppl, 1/4) # P(X >= N) = P(X > N-1) = 1 - P(X <= N-1)


def prtm(s : str):
    protein_masses = {
        'A': 71.03711,
        'C': 103.00919,
        'D': 115.02694,
        'E': 129.04259,
        'F': 147.06841,
        'G': 57.02146,
        'H': 137.05891,
        'I': 113.08406,
        'K': 128.09496,
        'L': 113.08406,
        'M': 131.04049,
        'N': 114.04293,
        'P': 97.05276,
        'Q': 128.05858,
        'R': 156.10111,
        'S': 87.03203,
        'T': 101.04768,
        'V': 99.06841,
        'W': 186.07931,
        'Y': 163.06333
    }

    total = 0
    for i in range(len(s)):
        total += protein_masses[s[i]]
    return total

def grph(s : str):
    # converts FATSA formatted file to dictionary
    dict = {}
    strings = s.split('>')[1:]
    for s in strings:
        splt = s.split('\n')
        dict[splt[0]] = ''.join(splt[1:])


    res = []
    for k1 in dict:
        for k2 in dict:
            if (k1 != k2) and (dict[k1].endswith(dict[k2][:3])):
                res.append(k1 + ' ' + k2)
    return res

def mprt(s : str):
    ids = s.split('\n')
    prot_strings = []
    for id in ids:
        url = 'http://www.uniprot.org/uniprot/%s.fasta' % id.split('_')[0]
        content = requests.get(url).text
        prot_strings.append(''.join(content.split('\n')[1:]))

    motif = r'(?=N[^P][ST][^P])'
    res = {}
    for i in range(len(prot_strings)):
        res[ids[i]] = [m.start()+1 for m in re.finditer(motif, prot_strings[i])]
    return res


def cons(s : str):
    # converts FATSA formatted file to dictionary
    all_codes = []
    strings = s.split('>')[1:]
    for s in strings:
        splt = s.split('\n')
        all_codes.append(''.join(splt[1:]))

    default = [0] * len(all_codes[0])
    matrix = {"A" : default[:],
              "C" : default[:],
              "G" : default[:],
              "T" : default[:], }
    for code in all_codes:
        for i in range(len(code)):
            matrix[code[i]][i] += 1

    res = ''
    for i in range(len(all_codes[0])):
        max_k, max_v = None, 0
        for k in matrix.keys():
            if matrix[k][i] > max_v:
                max_k, max_v = k, matrix[k][i]
        res += max_k
    return matrix, res

def orf(s : str):

    def orf_helper(s : str):
        # get indices of all start codon positions
        s = s.replace('T', 'U')
        indices = []
        for i in range(len(s)-2):
            protein = codon_table[s[i:i+3]]
            if protein and protein == 'M':
                indices.append(i)
        # get string at each position until stop codon reached
        print(indices)
        proteins = []
        for idx in indices:
            protein = ''
            for j in range(idx, len(s), 3):
                if s[j:j+3] in codon_table.keys():
                    if codon_table[s[j:j+3]] == 'Stop':
                        proteins.append(protein)
                        break
                    protein += codon_table[s[j:j+3]]
        return proteins

    s = ''.join(s.split('\n')[1:])
    return sorted(set(orf_helper(s) + orf_helper(revc(s))))

def splc(s : str):
    main_string = ''.join(s.split('>')[1].split()[1:])
    introns = []
    for item in s.split('>')[2:]:
        introns.append(''.join(item.split()[1]))
    print(main_string)
    print(introns)
    for i in introns:
        main_string = main_string.replace(i, '')
    main_string = main_string.replace('T', 'U')
    res = ''
    for i in range(0, len(main_string), 3):
        codon = main_string[i:i+3]
        if codon_table[codon] == 'Stop':
                break
        res += codon_table[codon]
    return res

def perm(n : int):
    nums = [i for i in range(1, n+1, 1)]  # generate list of [1, ..., n] first

    def perm_helper(nums : list):
        if len(nums) == 1:
            return [nums[:]]
        res = []
        for _ in range(len(nums)):
            n = nums.pop(0)
            perms = perm_helper(nums)
            for p in perms:
                p.append(n)
            res.extend(perms)
            nums.append(n)
        return res
    return perm_helper(nums)

if __name__ == '__main__':
    file_path = 'datasets/rosalind_orf.txt'
    # read in .txt file when input is one string
    with open(file_path, 'r') as file:
        input = file.read().strip()
    out = perm(1)
    print(len(out))
    for item in out:
        print(*item)