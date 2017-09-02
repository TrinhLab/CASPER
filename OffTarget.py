"""This file runs the off target analysis for CASPER. Use the CASPEROfflist.txt file to set up the information you want
    to run with this program.
    WARNING: Running this protocol on a large number of sequences is unwise and may take significant computing power/time."""


# -------------------------------------- USER CAN IGNORE CODE BELOW THIS LINE ---------------------------------------- #

import os, sys, math, datetime
from Bio import Seq
from Bio.Seq import Seq
from Algorithms import SeqTranslate
from Bio.Alphabet import IUPAC


class OffTargetAlgorithm:

    def __init__(self, threshold, endo, base_org, csf_file, other_orgs, casperofflist, output_path):
        self.ST = SeqTranslate()
        self.rSequences = []
        self.get_rseqs(casperofflist)
        self.mypath = csf_file[:csf_file.find(base_org)]
        self.ref_genomes = [base_org]
        self.ref_genomes += other_orgs
        self.endo = endo
        self.threshold = threshold
        self.dSequence = str()  # global to class so that all scoring functions can use it

        # This is for autofilling the HsuMatrix
        self.matrixKeys = ["GT", "AC", "GG", "TG", "TT", "CA", "CT", "GA", "AA", "AG", "TC", "CC"]
        self.matrix = {}
        self.fill_matrix()

        # This is where the data is stored before it is written
        self.output_data = dict()
        for myseq in self.rSequences:
            self.output_data[myseq[0]] = list()

        # BEGIN RUNNING THROUGH SEQUENCES
        for sequence in self.rSequences:
            print(sequence)
            for genome in self.ref_genomes:
                f = open(self.mypath + genome + self.endo + ".cspr", 'r')
                while True:
                    line = f.readline()
                    if line.find("CHROMOSOME") != -1:
                        curchrom = line[line.find("#") + 1:-1]
                        print("Finished checking " + curchrom)
                    else:
                        if line[0:-1] == "REPEATS":
                            break
                        # Checks for a signifcant number of mismatches:
                        #locseq = line[:-1].split(",")
                        if self.critical_similarity(sequence[0], self.ST.decompress_csf_tuple(line)[1]):
                            # This is where the real fun begins: off target analysis
                            print('found a similarity')
                            seqscore = self.get_scores(sequence[1],self.ST.decompress_csf_tuple(line)[1])
                            if seqscore > self.threshold:
                                self.output_data[sequence[0]].append((str(curchrom),
                                                                      self.ST.decompress_csf_tuple(line[:-1]),
                                                                      int(seqscore*100), genome))

        # END SEQUENCES RUN
        # Output the data acquired:
        out = open(output_path + "off_results" + str(datetime.datetime.now().time()) + '.txt', 'w')
        out.write("Off-target sequences identified.  Scores are between O and 1.  A higher value indicates greater"
                  "probability of off-target activity at that location.\n")
        for sequence in self.output_data:
            out.write(sequence + "\n")
            for off_target in self.output_data[sequence]:
                outloc = off_target[0] + "," + str(off_target[1][0]) + "," + off_target[1][1]
                out.write(off_target[3] + "," + outloc + "\t" + str(off_target[2]/100) + '\n')
        out.close()

    def get_rseqs(self, offlist):
        targets = list()
        cofile = open(offlist, 'r')
        cofile.readline()
        while True:
            t = cofile.readline()[:-1]
            if t == 'EN':
                break
            targets.append(t)
        for tar in targets:
            compseed = self.ST.compress(tar[:16], 64)
            comptail = self.ST.compress(tar[16:], 64)
            compressed = compseed + "." + comptail
            rseq = ""
            for nt in tar[0:-1]:
                rseq = nt + rseq
            self.rSequences.append([tar, rseq])

    def get_scores(self,rseq, dseq):
        self.dSequence = Seq(dseq, IUPAC.unambiguous_dna).reverse_complement()
        hsu = self.get_hsu_score(rseq)
        qual = self.get_qualt_score(rseq)
        step = self.qualt_step_score(rseq)
        output = ((math.sqrt(hsu) + step) + pow(qual, 6))
        return output

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

    # If there is more than four mismatches it returns false, else it will return true
    def critical_similarity(self, cseq1, cseq2):
        mismatches = 0
        lim = min([len(cseq1), len(cseq2)])
        check = True
        for i in range(lim):  # Doesn't matter whether you use cseq1 or cseq2 they are the same length
            if cseq1[i] != cseq2[i]:
                mismatches += 1
                if mismatches == 5:
                    check = False
                    break
        return check

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
