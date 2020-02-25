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


reads = fastq_to_reads('ads1_week4_reads.fq')
print("reads len:", len(reads))
print("unique reads:", len(set(reads)))

# scs = greedy_scs(reads, 30)
# print("scs len:", len(scs), "num 'A':", scs.count('A'), "num 'T':", scs.count('T'))
# len:15894, A:4633, T:3723
