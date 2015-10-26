import utils
import time
from collections import defaultdict

def edit_distance_alt(p, t):
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

    # Edit distance is the value in the bottom right corner of the matrix
    return min(D[-1])

def overlap(a, b, min_length=3):
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def overlap_graph(reads, k):
    # Make index
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i+k]].add(read)
    # Make graph
    graph = defaultdict(set)
    for r in reads:
        for o in index[r[-k:]]:
            if r != o:
                if overlap(r, o, k):
                    graph[r].add(o)
    edges = 0
    for read in graph:
        edges += len(graph[read])
    return(edges, len(graph))


def main():
    chr1 = utils.readGenome('chr1.GRCh38.excerpt.fasta')
    #Question 1
    start = time.clock()
    print edit_distance_alt('GCTGATCGATCGTACG', chr1)
    end = time.clock()
    print ">> %.2gs" % (end-start)
    #Question 2
    start = time.clock()
    print edit_distance_alt('GATTTACCAGATTGAG', chr1)
    end = time.clock()
    print ">> %.2gs" % (end-start)
    #Questions 3 and 4
    start = time.clock()
    seqs, _ = utils.readFastq('ERR266411_1.for_asm.fastq')
    edges, suffixes = overlap_graph(seqs, 30)
    print edges
    print suffixes
    end = time.clock()
    print ">> %.2gs" % (end-start)

if __name__ == '__main__':
    main()
