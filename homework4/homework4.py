import utils

def main():
    # question 1: What is the length of the shortest common superstring of the following strings?
    shortest_string, num_strings = utils.scs(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])
    print len(shortest_string)
    # question 2: How many different shortest common superstrings are there for the input strings given in the previous question?
    print num_strings

    # question 3 and 4:

    unknown_virus_seq, _ = utils.readFastq('ads1_week4_reads.fq')
    for k in range(100, 1, -1):
        genome = utils.greedy_scs(unknown_virus_seq, k)
        if len(genome) == 15894:
            # q3: How many As are there in the full, assembled genome?
            print(genome.count('A'))
            # q4: How many As are there in the full, assembled genome?
            print(genome.count('T'))
            # g5: final genome that we can search for in the BLAST database
            print(genome)
            break

if __name__ == '__main__':
    main()
