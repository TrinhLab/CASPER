"""This file is for generating data on multitargeting from a CASPER_Seq_Finder file.  It will simply decompress and
    print out in the command line in readable format the repeated sequences found in the organsim with the given
    endonuclease. Before doing anything, change the path object below to the appropriate path."""

# Please insert the appropriate path where you are storing your CASPER_Seq_Finder files

import operator
import os
from Algorithms import SeqTranslate


# ----------------------------------CODE BELOW CAN BE IGNORED BY USER------------------------------------------------ #


class Multitargeting:

    BAD_instances = {}
    sorted_instances = []

    def __init__(self, CASPER_Seq_Finder_file, output_path):
        self.file_name = CASPER_Seq_Finder_file
        self.get_instances()

    def get_instances(self):
            ST = SeqTranslate()
            os.chdir(path)
            f = open(self.file_name, 'r')
            while True:
                x = f.readline()
                if x == 'REPEATS\n':
                    print("reached repeat sequences")
                    break
            while True:
                t = f.readline()
                if t == 'END_OF_FILE':
                    print("reached end of repeat sequences")
                    break
                ukey = t[:-1]  # takes away the "\n" in the string
                key = ST.decompress64(ukey, True)
                key = ST.fill_As(key, 16)
                self.BAD_instances[key] = list()
                # Add sequences and locations to the list
                v = f.readline().split('\t')[:-1]
                for item in v:
                    loctup = item.split(',')
                    chrom = loctup[0]
                    location = ST.decompress64(loctup[1])
                    seq = ST.decompress64(loctup[2][1:],True)
                    seq = ST.fill_As(seq, 4)  # when A's get lost in the compression this fills them back in
                    mytup = (chrom, location, seq)
                    self.BAD_instances[key].append(mytup)
            f.close()
            print("currently sorting")
            for key in self.BAD_instances:
                size = len(self.BAD_instances[key])
                newtuple = (key, self.BAD_instances[key], size)  # sequence, location, size
                self.sorted_instances.append(newtuple)

    # Returns the container self.sorted_instances but removes all "single" repeats. Old Code to fix an off-by-1 error
    def return_all_seqs(self):
        myseqs = []
        for instance in self.sorted_instances:
            if instance[2] > 1:
                myseqs.append(instance)
        return myseqs

    def return_sorted(self):
        sorted_seqs = sorted(self.sorted_instances, key=operator.itemgetter(2), reverse=True)
        amounts = {}
        for instance in sorted_seqs:
            if instance[2] > 1:
                if instance[2] in amounts:
                    amounts[instance[2]] += 1
                else:
                    amounts[instance[2]] = 1
                print(str(instance[0]) + "," + str(instance[2]) + "," + str(instance[1]))
        for element in amounts:
            print("Number of seed sequences with " + str(element) + " appearances: " + str(amounts[element]))

    def return_positions(self):
        positions_mapped = []  # chromosme, beginning of range, end of range, and number of hits
        for instance in self.sorted_instances:
            if instance[2] > 1:
                for pos in instance[1]:
                    chrom = pos[0]
                    loc = int(pos[1])
                    # check to see if its already in the map
                    need_new = True
                    for position in positions_mapped:
                        if chrom == position[0]:
                            if position[1] < loc < position[2]:
                                position[3] += 1
                                position[4].append(instance[0])
                                need_new = False
                                print("position added")
                    if need_new:
                        newtuple = [chrom, loc-1000, loc+1000, 1, [" ", instance[0]]]
                        positions_mapped.append(newtuple)
        sorted_positions = sorted(positions_mapped, key=operator.itemgetter(3), reverse=True)
        for element in sorted_positions:
            print(str(element[0]) + "," + str(element[1]) + "," + str(element[2]) + "," + str(element[3]))
        for element in sorted_positions:
            sequences = ""
            for sequence in element[4]:
                sequences += sequence + ","
            print(sequences)
        return sorted_positions

    def int_to_char(self, i):
        switcher = {
            0: 'A',
            1: 'T',
            2: 'C',
            3: 'G'
        }
        return switcher[i]

#r = Multitargeting()
#r.return_sorted()
#r.return_positions()
