import matplotlib.pyplot as plt


def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def naive_with_rc(p, t):
    occ1 = naive(p, t)
    rc_p = reverseComplement(p)
    occ2 = naive(rc_p, t)
    return sorted(set(occ1 + occ2))


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip()  # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


# genome = readGenome('lambda_virus.fa')

# print(genome)
# print(len(naive('AGGT', genome)))
# print(len(naive('ACCT', genome)))
# print(len(naive_with_rc('AGGT', genome)))
# print(len(naive('TTAA', genome)))
# print(len(naive_with_rc('TTAA', genome)))
# print(naive_with_rc('ACTAAGT', genome)[0])
# print(naive_with_rc('AGTCGA', genome)[0])


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mm = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mm += 1
                if mm > 2:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

# print(naive_2mm('ACTTTA','ACTTACTTGATAAAGT'))


'''
phix_genome = readGenome('phix.fa')
occurrences = naive_2mm('GATTACA', phix_genome)
print('offset of leftmost occurrence: %d' % min(occurrences))
# offset of leftmost occurrence: 10
print('# occurrences: %d' % len(occurrences))
# occurrences: 79 
'''

# print(len(naive_2mm('TTCAAGCC', genome)))
# print(naive_2mm('AGGAGGTT', genome)[0])

seqs, quals = readFastq('ERR037900_1.first1000.fastq')
# print(quals[:10])


def phred33ToQ(qual):
    return ord(qual) - 33


def createHist(quals):
    # Create a histogram of quality scores
    hist = [0]*50
    for qual in quals:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist

# h = createHist(quals)
# print(h)
# Plot the histogram
# plt.plot(range(len(h)), h)
# plt.show()


def sum_quals(quals):
    # return a list of the total of the quality of each read
    # num_quals = len(quals)
    num_reads = len(quals[0])
    sum_quals = [0] * num_reads
    for qual in quals:
        # loop through all the quals and total the quality 
        # scores of each read position
        for i, q in enumerate(qual):
            sum_quals[i] += phred33ToQ(q)
    return sum_quals

sum_quals = sum_quals(quals)
min_sum = min(sum_quals)
print("min_sum is", min_sum)
print("which is at read number:", sum_quals.index(min_sum))


