from itertools import permutations


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]


def approxMatch(p, t):
    # Create distance matrix
    D = []
    for i in range(len(p)+1):
        D.append([0]*(len(t)+1))
    # Initialize first row and column of matrix
    for i in range(len(p)+1):
        D[i][0] = i
    for i in range(len(t)+1):
        D[0][i] = 0
    # Fill in the rest of the matrix
    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if p[i-1] == t[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # best approx. match will be the min value in bottom row
    # for row in D:
    #    print(row)
    return min(D[-1])


p, t = 'GCGTATGC', 'TATTGGCTATACGGTT'
# print(approxMatch(p, t)) # 2

t = readGenome('chr1.GRCh38.excerpt.fasta')
p = 'GCTGATCGATCGTACG'
print(approxMatch(p, t))  # 3

p = 'GATTTACCAGATTGAG'
print(approxMatch(p, t))  # 2


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


def naive_overlap_map(reads, k):
    olaps = {}
    for a, b in permutations(reads, 2):
        # print(a, b)
        olen = overlap(a, b, min_length=k)
        if olen > 0:
            olaps[(a, b)] = olen
    return olaps

# reads = ['ACGGATC', 'GATCAAGT', 'TTCACGGA']
# print(naive_overlap_map(reads, 3))


def overlap_map(fastq_file, k):
    with open(fastq_file, 'r') as f:
        reads = f.readlines()
        reads = [r.strip() for r in reads[1::4]]
    # reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
    # print("reads:", reads)

    # For every k-mer in a read, we add the read to the set object corresponding
    # to that k-mer
    d = {}  # kmer:set(reads)
    for r in reads:
        i = 0
        for i in range(len(r) - k + 1):
            kmer = r[i:i+k]
            if kmer not in d:
                d[kmer] = set()
            d[kmer].add(r)
    # print("d:", d)

    # Now, for each read 'a' we find all overlaps involving a suffix of 'a'.
    # To do this, we take a's length-k suffix, find all reads containing that
    # k-mer (obtained from the corresponding set) and call
    # overlap(a, b, min_length=k) for each.
    # The most important point is that we do not call overlap(a, b, min_length=k)
    # if 'b' does not contain the length-k suffix of 'a'.

    olaps = {}
    for a in reads:
        suff = a[-k:]
        for b in d[suff]:
            if a == b:
                continue
            olen = overlap(a, b, min_length=k)
            if olen >= k:
                olaps[(a, b)] = olen
    return olaps


# print(overlap_map('test.fastq', 4))
# print(overlap_map('test.fastq', 5))

infile = 'ERR266411_1.for_asm.fastq'
k = 30
olaps = (overlap_map(infile, k))
print("len of olaps:", len(olaps))
print(len(set(k[0] for k in olaps)))
