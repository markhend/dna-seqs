import itertools

def overlap(a, b, min_length=3):
    '''return length of longest suffix of 'a' matching
    a prefix of 'b' that is at least length 'min_length'
    characters long. Return 0 if no such overlap exists.
    '''
    start = 0  # start all the way to the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's pref in a
        if start == -1:  # no more occurrences to the right
            return 0
        # found occurrence, check for full suff/pref match
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1  # move past previous match


def scs_list(ss):
    '''return shortest common superstring of given
    strings, which must be the same length
    '''
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


ans1 = scs_list(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])
print("scs_list:", ans1, "len:", len(ans1), "str len:", len(ans1[0]))


