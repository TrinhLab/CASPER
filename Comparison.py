"""This file is for analyzing sequences across organisms.  Use this to find repeated sequences in a group of organisms.
    To run, type the list of organisms you want to run in the list object directly below.
    WARNING: If you want comparisons to contain annotations, you will need to have GenBank files for all the organisms
    labelled in the directory."""

# ---------------------------------- CODE BELOW CAN BE IGNORED BY USER ---------------------------------------------- #

import os

from Algorithms import SeqTranslate
import itertools


class Compare_Orgs:

    def __init__(self, output_path, base_org_path, base_org, endo, other_genomes):
        # initialize SeqTranslate object
        self.ST = SeqTranslate()
        self.output_path = output_path
        # my_orgs contains just the
        self.organisms = other_genomes
        self.organisms.append(base_org)
        self.organisms = sorted(self.organisms)
        self.db_path = base_org_path[:base_org_path.find(base_org)]

        # Dictionary of dictionaries. Key1: generic total sequence Key2: org Value: position
        self.searchableseqs = {}

        # Container that stores all the sequences seen the combination of organisms defined by the key
        # An example key would be (sce, yli) for the shared sequences between S.cerevisiae and Y.lipolytica
        self.buckets = {}

        # Intitialize the self.buckets container to contain the tuple of every organism subset
        for i in range(2, len(self.organisms)):
            for subset in itertools.combinations(self.organisms, i):
                self.buckets[subset] = []
                print(subset)
        self.endo = endo

        # The object that is iterated over to decompress the output into readable form
        self.compressed_output = {}

        # Generates the sequence lists
        for org in self.organisms:
            print(org)
            self.make_lists(org)

        # Runs the comparison
        self.create_comparison()
        self.write_to_file()

    # Takes an organism and parses the target data into positions and repeated sequences containers
    def make_lists(self, org):
        name1 = self.db_path + org + self.endo + ".cspr"
        f = open(name1, 'r')
        curchrom = int()
        while True:
            position = f.readline()
            if position.find("CHROMOSOME") != -1:
                curchrom = position[position.find("#")+1:-1]
                print(curchrom)
            else:
                if position[0:-1] == "REPEATS":
                    break
                # adds to the positions container the unique position with organism, and chromosome as keys
                line = position[:-1].split(",")
                # change line into generic (no "+" or "-" change to generic .)
                totseq = self.ST.to_generic_compressed(line[1])
                self.add_to_sequence_matrix(totseq, org, curchrom, line[0])

        while True:
            seedseq = f.readline()[:-1]
            if seedseq.find("END_OF_FIL") != -1:
                break
            taillocs = f.readline().split('\t')[:-1]
            for item in taillocs:
                loctup = item.split(',')
                totseq = self.ST.to_generic_compressed([seedseq,loctup[2][1:]])
                self.add_to_sequence_matrix(totseq, org, loctup[0], loctup[1])
        f.close()

    # Takes in the variables of a sequence including what organism it is found on and adds it to the dict of dicts
    # named:self.searchableseqs
    def add_to_sequence_matrix(self, totseq, org, chrom, location):
        if totseq in self.searchableseqs.keys():
            # already seen this organism and sequence
            if org in self.searchableseqs[totseq].keys():
                self.searchableseqs[totseq][org].append((chrom, location))
            # already seen this sequence but not this sequence in the organism
            else:
                self.searchableseqs[totseq][org] = []
                self.searchableseqs[totseq][org].append((chrom, location))
        # new organism and new sequence
        else:
            self.searchableseqs[totseq] = {}
            self.searchableseqs[totseq][org] = []
            self.searchableseqs[totseq][org].append((chrom, location))

    def int_to_char(self, i):
        switcher = {
            0: 'A',
            1: 'T',
            2: 'C',
            3: 'G'
        }
        return switcher[i]

    def revcom(self, sequence):
        revseq = ""
        change = {'A':'T',
                  'T':'A',
                  'G':'C',
                  'C':'G',
                  'N':'N'}
        for nt in sequence:
            if nt in change:
                rnt = change[nt]
            else:
                rnt = nt
            revseq = rnt + revseq
        return revseq

    def create_comparison(self):
        # Put every sequence in the appropriate bucket
        tempdict = dict()
        for sequence in self.searchableseqs:
            # Look for the set of organisms containing this sequence
            if len(self.searchableseqs[sequence].keys()) > 1:
                # Make sure the tuple is in the right order
                orgs = self.searchableseqs[sequence].keys()
                orgs = tuple(sorted(orgs))
                # iterate through all the sequences contained in each organism
                for org in self.searchableseqs[sequence].keys():
                    tempdict[org] = []
                    for location in self.searchableseqs[sequence][org]:
                        tempdict[org].append(location)
                insert_tup = (sequence, tempdict)
                tempdict = {}
                # contains a list of tuples with the sequence then short dictionary of organisms containing sequence
                self.buckets[orgs].append(insert_tup)

    def write_to_file(self):
        filename = self.output_path + "compare_"
        for org in self.organisms:
            filename += org + "_"
        filename += self.endo + ".txt"
        f = open(filename, 'w')
        for key in self.buckets:
            f.write(str(key) + " " + str(len(self.buckets[key])) + "\n")
            for seq in self.buckets[key]:
                f.write(self.ST.decompress64(seq[0], True) + "\n")
                for suborg in seq[1]:
                    f.write(str(suborg) + ":")
                    for locs in seq[1][suborg]:
                        f.write(str(locs[0]) + "," + str(self.ST.decompress64(locs[1])) + "\t")
                        f.write("\n")
            f.write("\n")
        f.close()

