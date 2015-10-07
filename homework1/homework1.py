
# === Algorithms for DNA Sequencing ====
# On Coursera, provided by John Hopkins School of Medicine
# Programming Homework 1


# Methods developed during the lectures:
# call methods eg: basic_functions.readGenome(filename)
import basic_functions

# Custom methods 

#  1. naive_with_rc returns the same results when P equals its reverse complement

def naive_with_rc(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        for double_seq in (p, basic_functions.reverseComplement(p)): # check the original and complement
            match = True
            for j in range(len(double_seq)):  # loop over characters
                if t[i+j] != double_seq[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record and 
                break # break, no need to check the complement
    return occurrences

# 2. naive_2mm allows up to 2 mismatches per occurrence. 
# As per the question problem we do not consider the reverse complement here
# because We're looking for approximate matches for P itself, not its reverse complement. 

def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatches = 0 # counter
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.append(i)  # all chars matched; record
    return occurrences

# 3. The dataset has something wrong with it; one of the sequencing cycles is poor quality.
# With this function we report which sequencing cycle has the problem (return the minimum quality index). 

def lowest_quality_base(qualities):
    total = [0] * len(qualities[0])
    for q in qualities:
        for i, phred in enumerate(q):
            total[i] += phred33ToQ(phred)
    return total.index(min(total))


def phred33ToQ(quality):
    return ord(quality) - 33

# main() function:

def main():
    # read the .fa and .fastq files
    lambda_virus = basic_functions.readGenome('lambda_virus.fa')
    _, qualities = basic_functions.readFastq('ERR037900_1.first1000.fastq')
    # problem 1:
    print(len(naive_with_rc('AGGT', lambda_virus)))
    # problem 2:
    print(len(naive_with_rc('TTAA', lambda_virus)))
    # problem 3:
    print(naive_with_rc('ACTAAGT', lambda_virus)[0])
    # problem 4:
    print(naive_with_rc('AGTCGA', lambda_virus)[0])
    # problem 5:
    print(len(naive_2mm('TTCAAGCC', lambda_virus)))
    # problem 6:
    print(naive_2mm('AGGAGGTT', lambda_virus)[0])
    # problem 7:
    print(lowest_quality_base(qualities))

if __name__ == '__main__':
    main()
