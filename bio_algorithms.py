
import numpy as np

# ----------Read Fasta File----------


infile = open('HAPPENN_dataset.fasta')

sequence = ""
for line in infile:
    if not line.startswith(">"):
        sequence += line.strip()

infile.close()


# ----------GC Content----------

def GC_Content(seq):
    l=len(seq)
    num_G=seq.count("G")
    num_C=seq.count("C")
    total=num_C+num_G
    return total/l

# ----------Reverse----------

def Reverse(seq):
    s=list(seq)
    s=reversed(s)
    s="".join(s)
    return s

# ----------Complement----------

def Complement(seq):
    dic = {
        "A":"T", "T":"A",
        "C":"G", "G":"C",
        "N":"N", "Y":"Y", "R":"R", "S":"S", "W":"W", "K":"K", "M":"M"
    }
    return "".join(dic.get(base, base) for base in seq)

# ----------Reverse-Complement----------

def Reverse_Complement(seq):
    seq=Reverse(seq)
    seq=Complement(seq)
    return seq

# ----------DNA to Amino Acid Translation----------

CODON_TABLE = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}


def translate_dna(sequence):
    protein = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        protein.append(CODON_TABLE.get(codon, 'X'))

    return ''.join(protein)

# ----------Boyre Moore----------

# ----------Bad Character----------


# def bad_character(seq,sub_seq):
#     table=np.zeros([4,len(sub_seq)])
#     row=["A","C","G","T"]
#     for i in range (4):
#         num=-1
#         for j in range (len(sub_seq)):
#             if row[i]==sub_seq[j]:
#                 table[i,j]=-1
#                 num=-1
#             else:
#                 num+=1
#                 table[i,j]=num
#     x=-1
#     i=0
#     while(i<len(seq)-len(sub_seq)+1):
#         if sub_seq==seq[i:i+len(sub_seq)]:
#             x=i
#             break

#         else:
#             for j in range(len(sub_seq)-1,-1,-1):
#                 if seq[i+j] != sub_seq[j]:
#                     k=row.index(seq[i+j])
#                     i+=table[k,j]
#                     break
#         i=int(i+1)
#     return x

def bad_character(P):
    """
    Build bad character table for Boyer-Moore.
    Returns dictionary mapping base -> last index in pattern.
    """
    bc = {}
    for i, c in enumerate(P):
        bc[c] = i
    return bc

# ----------Good Suffix----------

def good_suffix_table(P):
    m = len(P)
    if m == 0:
        return []
    shift = [m] * m
    for i in range(m - 1):
        j = i
        k = m - 1
        while j >= 0 and P[j] == P[k]:
            j -= 1
            k -= 1
        shift[m - 1 - i] = m - 1 - k
    return shift
# ----------Boyer-Moore Search----------

def boyer_moore_search(text, pattern):
    if len(pattern) == 0 or len(text) == 0:
        return []

    bc = bad_character(pattern)
    gs = good_suffix_table(pattern)
    m = len(pattern)
    n = len(text)
    i = 0
    matches = []

    while i <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[i + j]:
            j -= 1
        if j < 0:
            matches.append(i)
            i += gs[0] if gs else 1
        else:
            bc_shift = j - bc.get(text[i + j], -1)
            gs_shift = gs[j + 1] if j + 1 < len(gs) else 1
            i += max(bc_shift, gs_shift)
    return matches

# ----------Indexed Matching----------

def index_sorted(seq, k):
    index = []
    for i in range(len(seq) - k + 1):
        index.append((seq[i:i+k], i))
    index.sort()
    return index

import bisect
from itertools import permutations
def query_indexed(text, pattern, index):
    k = len(index[0][0])
    keys = [r[0] for r in index]
    start = bisect.bisect_left(keys, pattern[:k])
    end = bisect.bisect_right(keys, pattern[:k])
    hits = index[start:end]


    offsets = []
    for _, pos in hits:
        if text[pos:pos+len(pattern)] == pattern:
            offsets.append(pos)
    return offsets

# ----------Suffix Array + Inverse Suffix Array----------

def Suffix_array_construction(T):
    result = []
    dec = {
    '$' : 0,
    'A' : 1,
    'C' : 2,
    'G' : 3,
    'T' : 4
    }
    table=[]
    i=2**0
    n=0
    while True:
        l=[]
        dec2={}
        if i>1:
            for j in range(len(T)):
                if not(table[n-1][j:j+i] in l):
                    l.append(table[n-1][j:j+i])
            l.sort()
            l
            for j in range(len(l)):
                dec2[tuple(l[j])]=j
        row=[]
        for j in range (len(T)):
            if i==1:
                row.append(dec[T[j]])
            else:
                row.append(dec2[tuple(table[n-1][j:j+i])])
        table.append(row)
        flag=0
        for j in range(len(row)):
            c=0
            c=row.count(j)
            if c>1:
                flag=1
                break
        result.append(f'iteration {n}: {row}')
        if flag==0:
            break
        n+=1
        i=2**n
    return result

def inverse_suffix_array(sa):
    inv = [0] * len(sa)
    for i, pos in enumerate(sa):
        inv[pos] = i
    return inv

def suffix_array(s):
        return sorted(range(len(s)), key=lambda i: s[i:])

def inverse_suffix_array(sa):
    inv = [0] * len(sa)
    for i, pos in enumerate(sa):
        inv[pos] = i
    return inv

# ----------Overlap & Overlap Graph----------

def overlap(a, b, min_length=1):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def overlap_graph(reads, k):
    graph = {}
    for a, b in permutations(reads, 2):
        olen = overlap(a, b, k)
        if olen > 0:
            graph.setdefault(a, []).append((b, olen))
    return graph

# ----------Approximate Matching----------

# ----------Hamming Distance----------


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strings must be equal length")

    distance = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            distance += 1
    return distance

# ----------Edit Distance----------

def edit_distance(x,y):
    d=[]
    for i in range(len(x)+1):
        d.append([0]*(len(y)+1))
    for i in range(len(y)+1):
        d[i][0]=i
    for i in range(len(y)+1):
        d[0][i]=i
    for i in range(1,len(x)+1):
        for j in range (1,len(y)+1):
            delta=1 if x[i-1] !=y[j-1] else 0
            d[i][j]= min(d[i-1][j-1]+delta, d[i][j-1]+1, d[i-1][j]+1)

    return d[-1][-1]


# ----------Hamming Distance----------

def hamming_distance(seq1, seq2):
    """Calculate Hamming distance between two sequences"""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have equal length")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


# ----------Read FASTA File----------

def read_fasta(filepath):
    """Read FASTA file and return list of (header, sequence) tuples"""
    sequences = []
    current_header = None
    current_seq = ""
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, current_seq))
                current_header = line[1:]
                current_seq = ""
            else:
                current_seq += line
        
        if current_header is not None:
            sequences.append((current_header, current_seq))
    
    return sequences