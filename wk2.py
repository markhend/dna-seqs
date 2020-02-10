def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mm = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters and allow up to 2 mismatches
                mm += 1
                if mm > 2:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

# implement versions of the naive exact matching and 
# Boyer-Moore algorithms that additionally count and return 
# (a) the number of character comparisons performed and (b) 
# the number of alignments tried. Roughly speaking, these measure 
# how much work the two different algorithms are doing.

# Implement naive_with_counts by extending naive function
# from naive_with_counts import naive_with_counts

def naive_with_counts(p, t):
    occurrences = []
    num_aligns, num_chars = 0, 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        num_aligns += 1
        match = True
        for j in range(len(p)):  # loop over characters
            num_chars += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, num_aligns, num_chars


# Example 1
# p = 'word'
# t = 'there would have been a time for such a word'
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)
# ([40], 41, 46)

# Example 2
# p = 'needle'
# t = 'needle need noodle needle'
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)
# ([0, 19], 20, 35)


# Implement boyer_moore_with_counts by extending boyer_moore function
# from bm_with_counts import boyer_moore_with_counts
from bm_preproc import BoyerMoore

def boyer_moore_with_counts(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    num_aligns, num_chars = 0, 0
    while i < len(t) - len(p) + 1:
        num_aligns += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            num_chars += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, num_aligns, num_chars

# Example 1
# p = 'word'
# t = 'there would have been a time for such a word'
# lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
# p_bm = BoyerMoore(p, lowercase_alphabet)
# occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
# print(occurrences, num_alignments, num_character_comparisons)
# ([40], 12, 15)

# Example 2
# p = 'needle'
# t = 'needle need noodle needle'
# p_bm = BoyerMoore(p, lowercase_alphabet)
# occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
# print(occurrences, num_alignments, num_character_comparisons)
# ([0, 19], 5, 18)
print()


# 1. How many alignments does the naive exact matching algorithm try when matching the string 
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the 
# excerpt of human chromosome 1? (Don't consider reverse complements.)
# 2. How many character comparisons does the naive exact matching algorithm try when matching the string 
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the 
# excerpt of human chromosome 1? (Don't consider reverse complements.)

t = readGenome('chr1.GRCh38.excerpt.fasta')
p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print("occurrences, num_alignments, num_character_comparisons")
print("naive:", occurrences, num_alignments, num_character_comparisons)

# 3. How many alignments does Boyer-Moore try when matching the string 
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the 
# excerpt of human chromosome 1? (Don't consider reverse complements.) 

p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG" 
p_bm = BoyerMoore(p)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print("bm:", occurrences, num_alignments, num_character_comparisons)
print()

# 4. Index-assisted approximate matching. In practicals, we built a Python class called Index 
# implementing an ordered-list version of the k-mer index. The Index class is copied below.
# We also implemented the pigeonhole principle using Boyer-Moore as our exact matching algorithm.
# Implement the pigeonhole principle using Index to find exact matches for the partitions. 
# Assume P always has length 24, and that we are looking for approximate matches with up to 2 mismatches 
# (substitutions). We will use an 8-mer index.

import bisect

class Index(object):
    """ Holds a substring index for a text T """

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

# Write a function that, given a length-24 pattern P and given an Index object built on 8-mers, 
# finds all approximate occurrences of P within T with up to 2 mismatches. Insertions and deletions 
# are not allowed. Don't consider any reverse complements.

# def approximate_match(p, t, n):
def approximate_matches(p, t, index):
    # segment_length = int(round(len(p) / (n+1)))
    segment_length = 8
    all_matches = set()
    # for i in range(n+1):
    num_index_hits = 0
    for i in range(3):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        print("i:", i, " start:", start, " end:", end)
        # p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        # matches = boyer_moore(p[start:end], p_bm, t)
        matches = index.query(p[start:end])
        print("matches:", matches, "len:", len(matches), '\n')
        num_index_hits += len(matches)
        # Extend matching segments to see if whole p matches
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > 2:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > 2:
                        break
            if mismatches <= 2:
                all_matches.add(m - start)
    print("index hits", num_index_hits, '\n')
    return sorted(list(all_matches))

# How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, which is derived from a human Alu sequence, 
# occur with up to 2 substitutions in the excerpt of human chromosome 1? (Don't consider reverse 
# complements here.)
# Hint 1: Multiple index hits might direct you to the same match multiple times, but be careful not to count 
# a match more than once.
# Hint 2: You can check your work by comparing the output of your new function to that of the 
# naive_2mm function implemented in the previous module.  

index = Index(t, 8)
p = "GGCGCGGTGGCTCACGCCTGTAAT"
all_matches = approximate_matches(p, t, index)
print("4. all_matches:", all_matches, "len ", len(all_matches))
print()
print("naive_2mm:", naive_2mm(p, t))
print()

# 5. Using the instructions given in Question 4, how many total index hits are there when searching 
# for occurrences of GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of human 
# chromosome 1? (Don't consider reverse complements.) Hint: You should be able to use the boyer_moore 
# function (or the slower naive function) to double-check your answer. 


# 6. Let's examine whether there is a benefit to using an index built using subsequences of T rather 
# than substrings, as we discussed in the "Variations on k-mer indexes" video. We'll consider subsequences 
# involving every N characters. For example, if we split ATATAT into two substring partitions, we would 
# get partitions ATA (the first half) and TAT (second half). But if we split ATATAT into two subsequences 
# by taking every other character, we would get AAA (first, third and fifth characters) and TTT 
# (second, fourth and sixth).
# Another way to visualize this is using numbers to show how each character of P is allocated to a partition. 
# Splitting a length-6 pattern into two substrings could be represented as 111222, and splitting into two 
# subsequences of every other character could be represented as 121212
# The following class SubseqIndex is a more general implementation of Index that additionally handles 
# subsequences. It only considers subsequences that take every Nth character:

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

# Write a function that, given a length-24 pattern P and given a SubseqIndex object built with 
# k = 8 and ival = 3, finds all approximate occurrences of P within T with up to 2 mismatches.
# When using this function, how many total index hits are there when searching for GGCGCGGTGGCTCACGCCTGTAAT 
# with up to 2 substitutions in the excerpt of human chromosome 1? 
# (Again, don't consider reverse complements.)

def num_char_mm(s, t):
    # return the number of mismatched chars between 2 equal length strings
    if len(s) != len(t):
        print("strings not same length")
        return
    mm = 0
    for a, b in zip(s, t):
        if a != b:
            mm +=1
    return mm

def query_subseq(p, t, subseq_index):
    # return occurrences of lenth 24 p in t and num_index_hits, k = 8
    k, ival = 8, 3
    span = 1 + ival * (k - 1)
    total_index_hits = 0
    matches = []
    p_subs = []
    for i in range(3):
        p_subs.append(p[i:i+span:ival])
    print("p:",p)
    print("subseqs:", p_subs)
    for i in range(3):
        hits = subseq_index.query(p[i:])
        total_index_hits += len(hits)

        print("i:", i)
        print("hits", hits)
        print("total_index_hits:", total_index_hits)
        print()

        for h in hits:
            mismatches = 0
            if i == 0:
                for j in [1, 2]:
                    mismatches +=  num_char_mm(p_subs[i+j], t[h+j:h+j+span:ival])
            elif i == 1:
                for j in [-1, 1]:
                    mismatches +=  num_char_mm(p_subs[i+j], t[h+j:h+j+span:ival])
            elif i == 2:
                for j in [-2, -1]:
                    mismatches +=  num_char_mm(p_subs[i+j], t[h+j:h+j+span:ival])
            if mismatches <= 2 and h-i not in matches:
                matches.append(h-i)
        print("matches", sorted(matches), "len:", len(matches))
        print()
    return


p = "GGCGCGGTGGCTCACGCCTGTAAT"
subseq_index = SubseqIndex(t, 8, 3)
print("5.")
query_subseq(p, t, subseq_index)

t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_index = SubseqIndex(t, 8, 3)
query_subseq(p, t, subseq_index)

t = open('1110.txt.utf-8').read()
p = 'English measure backward'
subseq_index = SubseqIndex(t, 8, 3)
query_subseq(p, t, subseq_index)
