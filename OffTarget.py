"""This file runs the off target analysis for CASPER. Use the CASPEROfflist.txt file to set up the information you want
    to run with this program.
    WARNING: Running this protocol on a large number of sequences is unwise and may take significant computing power/time."""

# Please set the location of your CASPEROfflist file!
casperofflist = " "  # This is the complete path of your CASPEROfflist.txt file
outfilepath = ""  # This is where your output file will go.  If left blank then it will go into the local directory


# -------------------------------------- USER CAN IGNORE CODE BELOW THIS LINE ---------------------------------------- #

import os, sys, math, datetime
from Bio import Seq
from Bio.Seq import Seq
from Algorithms import SeqTranslate
from Bio.Alphabet import IUPAC


class OffTargetAlgorithm:

    def __init__(self, targets, genomes):
        self.rSequences = []
        self.get_rseqs(targets)
        self.dna_seqs = list()
        self.fill_dna(genomes)
        self.matrixKeys = ["GT", "AC", "GG", "TG", "TT", "CA", "CT", "GA", "AA", "AG", "TC", "CC"]
        self.matrix = {}
        self.fill_matrix()

        outfile = open(outfilepath+"offtargetresults"+str(datetime.datetime.now()), 'w')
        for seq in self.rSequences:
            for dna_seq in self.dna_seqs:
                self.dSequence = Seq(dna_seq[1], IUPAC.unambiguous_dna).reverse_complement()
                hsu = self.get_hsu_score(seq)
                qual = self.get_qualt_score(seq)
                step = self.qualt_step_score(seq)
                output = (math.sqrt(hsu) + step)*pow(qual,6)
                if output > 0.3:
                    outfile.write(dna_seq + output)
        outfile.close()

    def get_rseqs(self, targets):
        for tar in targets:
            rseq = ""
            for nt in tar[0:-1]:
                rseq = nt + rseq
            self.rSequences.append(rseq)
        print(str(len(self.rSequences)))

    def fill_dna(self, genomes):
        ST = SeqTranslate()
        for genome in genomes:
            f = open(genome)
            fline = f.readline()
            while line[:-1] != "REPEATS":
                if line.find("CHROMOSOME") != -1:
                    mytup = fline[:-1].split(',')
                    loc = ST.decompress64(mytup[0])
                    seq = mytup[1][0] + ST.decompress64(mytup[1][1:], True)
                    mytup = (loc,seq)
                    self.dna_seqs.append(mytup)
            f.close()

    def fill_matrix(self):
        f = open('CASPERinfo', 'r')
        l = " "
        while True:
            l = f.readline()
            if l[0] == "H":
                break
        i = 0
        l = f.readline()
        while l[0] != '-':
            values = l.split("\t")
            self.matrix[self.matrixKeys[i]] = values
            i += 1
            l = f.readline()
        for element in self.matrix:
            self.matrix[element][18] = self.matrix[element][18][0:-1]

    def get_hsu_score(self, rSequence):
        score = 1.0
        for i in range(0,19):
            rnt = rSequence[i]
            dnt = self.dSequence[i]
            lookup = str(rnt) + str(dnt)
            if lookup in self.matrixKeys:
                hsu = self.matrix[lookup][18-i]
                score *= float(hsu)
        return score

    def get_qualt_score(self, rSequence):
        score = 3.5477
        for i in range(0, 19):
            lookup = rSequence[i] + self.dSequence[i]
            if lookup in self.matrixKeys:
                score -= 1.0/(i+1)
        return score/3.5477

    def qualt_step_score(self, rSequence):
        score = 1.0
        for i in range(0, 19):
            lookup = rSequence[i] + self.dSequence[i]
            if lookup in self.matrixKeys:
                if i < 6:
                    score -= 0.1
                elif i < 12:
                    score -= 0.05
                elif i < 20:
                    score -= 0.0125
        return score

    def separation_score(self, rSequence):
        misses = []
        delta = 0
        for i in range(0, 19):
            lookup = rSequence[i] + self.dSequence[i]
            if lookup in self.matrixKeys:
                misses.append(i)
        if len(misses) == 2:
            delta = (misses[1] - misses[0])/2.0
        if len(misses) == 3:
            delta = ((misses[1] - misses[0]) + (misses[2] - misses[1]))/3.0
        if len(misses) == 4:
            delta = ((misses[1] - misses[0]) + (misses[2] - misses[1]))/3.0
        retval = 1.0 - (delta/19.0)
        return retval

    def int_to_char(self, i):
        switcher = {
            0: 'A',
            1: 'T',
            2: 'C',
            3: 'G'
        }
        return switcher[i]

    def char_to_int(self, c):
        switcher = {
            'A': 0,
            'T': 1,
            'C': 2,
            'G': 3
        }
        return switcher[c]

f = open(casperofflist)
sequences = list()
targetgenomes = list()
line = f.readline()
while line.find('Add') != -1:
    line = f.readline()[:-1]
    sequences.append(line)
while line[0:3] != 'END':
    line = f.readline()[:-1]
    targetgenomes.append(line)

x = OffTargetAlgorithm(sequences, targetgenomes)
