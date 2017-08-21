import Tkinter as tk
import operator
import re
from threading import Thread
import tkMessageBox
import subprocess as sub
from timeit import timeit
from bioservices import UniProt
from apis import Kegg
import os

class Random_Mutagenesis:

    path = "/Users/brianmendoza/Desktop/CrisprDB/"
    BAD_instances = {}
    sorted_instances = []

    def __init__(self, org, endo, isanti):
        self.file_name = org + endo + ".txt"
        self.anti = isanti
        self.get_instances()

    def get_instances(self):
            os.chdir(self.path)
            f = open(self.file_name, 'r')
            while True:
                x = f.readline()
                if x[0:3] == 'BAD':
                    print("reached bad sequences")
                    break
            isnewseq = True
            while True:
                if isnewseq:
                    t = f.readline()
                if t[0:3] == 'NAG':
                    print("reached end of bad sequences")
                    break
                ukey = t
                key = self.decompress(ukey)[0:-1]  # takes away the "\n" in the string
                self.BAD_instances[key] = []
                while True:
                    v = f.readline()
                    if "," in v:
                        value = tuple(v.split(","))
                        self.BAD_instances[key].append(value)
                    else:
                        t = v
                        isnewseq = False
                        break
            f.close()
            print("currently sorting")
            for key in self.BAD_instances:
                size = len(self.BAD_instances[key])
                newtuple = (key, self.BAD_instances[key], size)  # sequence, location, size
                self.sorted_instances.append(newtuple)

    # Returns the container self.sorted_instances but removes all "single" repeats
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
                if amounts.has_key(instance[2]):
                    amounts[instance[2]] += 1
                else:
                    amounts[instance[2]] = 1
                print str(instance[0]) + "," + str(instance[2])
        for element in amounts:
            print str(element) + "," + str(amounts[element])

                # print instance[2]
                # for pos in instance[1]:
                    # print pos[0] + "," + pos[1]

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
                                print "position added"
                    if need_new:
                        newtuple = [chrom, loc-1000, loc+1000, 1, [" ", instance[0]]]
                        positions_mapped.append(newtuple)
        sorted_positions = sorted(positions_mapped, key=operator.itemgetter(3), reverse=True)
        for element in sorted_positions:
            print str(element[0]) + "," + str(element[1]) + "," + str(element[2]) + "," + str(element[3])
        for element in sorted_positions:
            sequences = ""
            for sequence in element[4]:
                sequences += sequence + ","
            print sequences
        return sorted_positions



    def decompress(self, compressed_seq):
        matrixString = '!m#$%&()*+j-./23456789:;<=>?@BDEFHIJKLMNOPQRSkVWXYZ[]^_abcdefghi'
        seq = ''
        for c in compressed_seq:
            index = matrixString.find(c)
            if index == -1:
                triad = c
            else:
                z = index % 4
                y = (index/4) % 4
                x = index/16
                triad = self.int_to_char(x) + self.int_to_char(y) + self.int_to_char(z)
            seq += triad
        return seq

    def int_to_char(self, i):
        switcher = {
            0: 'A',
            1: 'T',
            2: 'C',
            3: 'G'
        }
        return switcher[i]

# r = Random_Mutagenesis()
# r.return_sorted()
# r.return_positions()
