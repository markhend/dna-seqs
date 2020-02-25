import itertools


def fastq_to_reads(file):
    with open(file, 'r') as f:
        reads = f.readlines()
        reads = [r.strip() for r in reads[1::4]]
        return reads


cache = {}

def overlap(a, b, min_length=3):
    ''' return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least length 'min_length'
        characters long. Return 0 if no such overlap exists. '''
    if (a, b, min_length) in cache:
        return cache[(a, b, min_length)]
    start = 0  # start all the way to the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's pref in a
        if start == -1:  # no more occurrences to the right
            cache[(a, b, min_length)] = 0
            return 0
        # found occurrence, check for full suff/pref match
        if b.startswith(a[start:]):
            ret = len(a) - start
            cache[(a, b, min_length)] = ret
            return ret
        start += 1  # move past previous match


def overlap_map(reads, k):
    # reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']

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

    olaps = {}
    for a in reads:
        suff = a[-k:]  # looking for suffix of 'a'
        for b in d[suff]:  # only check other reads that contain the suff of 'a'
            if a == b:
                continue
            olen = overlap(a, b, min_length=k)
            if olen >= k:
                olaps[(a, b)] = olen
    return olaps


def scs_list(ss):
    ''' return shortest common superstring of given
        strings, which must be the same length '''
    shortest_sup = [None]
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss) - 1):
            # overlap adjacent strings in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of b to sup
            sup += ssperm[i+1][olen:]
        if shortest_sup[0] is None or len(sup) < len(shortest_sup[0]):
            shortest_sup = [sup]
        else:
            if len(sup) == len(shortest_sup[0]):
                shortest_sup.append(sup)
    return shortest_sup


def pick_maximal_overlap(reads, k):
    ''' Return a pair of reads from the list with a
        maximal suffix/prefix overlap >= k.  Returns
        overlap length 0 if there are no such overlaps.'''
    reada, readb = None, None
    best_olen = 0
    for a, b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen


def greedy_scs(reads, k):
    ''' Greedy shortest-common-superstring merge.
        Repeat until no edges (overlaps of length >= k)
        remain. '''
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)


def greedy_scs_from_olaps(reads, olaps):
    '''Create a scs from a dict of (a, b):olen'''
    olaps = [(k) + (v,) for k, v in olaps.items()]  # convert dict to list of tuples
    olaps = sorted(olaps, key=lambda x: x[2], reverse=True)  # reverse sort list by olen
    print(olaps[-2:])
    for read_a, read_b, olen in olaps:
        if read_a not in reads or read_b not in reads:
            continue
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
    return ''.join(reads)


reads = fastq_to_reads('ads1_week4_reads.fq')
print("reads len:", len(reads))
print("unique reads:", len(set(reads)))
# olaps = overlap_map(reads, 8)
olaps = overlap_map(reads, 30)
olaps = [(k) + (v,) for k, v in olaps.items()]
olaps = sorted(olaps, key=lambda x: x[2])
print(len(olaps))

# scs = greedy_scs(reads, 30)
# print("scs len:", len(scs), "num 'A':", scs.count('A'), "num 'T':", scs.count('T'))
# len:15894, A:4633, T:3723


# scs = greedy_scs_from_olaps(reads, olaps)
# print("scs len:", len(scs), "num 'A':", scs.count('A'), "num 'T':", scs.count('T'))
