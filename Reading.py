__author__ = 'brianmendoza'

import os
import operator
from testingfile import Random_Mutagenesis
from GenBankParse import GenBankFile


class AnnotationTable:

    def __init__(self, organism):
        self.org = organism

        # Output Containers
        self.repeats = []  # this should be filled by the table of positions on the excel document
        self.feature_targets = {}
        self.genes = {}

        self.get_positions()
        g = GenBankFile(self.org)
        # dictionary with chromosome as keys. Value is list with tuples in form: (srt, end, locus_tag, product)
        self.annotations = g.parseAnnotation()

        self.analyze_sequences()
        self.write_data_files()

    # use dictionary to but into chromosome buckets to make searching faster

    def get_positions(self):
        r = Random_Mutagenesis(self.org, "spCas9", False)
        self.repeats = r.return_all_seqs()  # entries exist in a list in the form (sequence, [locations], size)
        print "got instances"

    def translate_chromosome(self, chr):
        numbers = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16')
        letters = ('A','B','C','D','E','F','G','H')
        roman = ('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII', 'XIV',
                 'XV', 'XVI')
        types =[numbers,letters,roman]
        for list in types:
            if chr in list:
                ind = list.index(chr)
                return numbers[ind]
        return '1'

    def analyze_sequences(self):
        for sequence in self.repeats:
            realseq = sequence[0]
            self.feature_targets[realseq] = []  # establishes list for the sequence with the sequence as the key
            for position in sequence[1]:
                chromosome = position[0]
                for feature in self.annotations[chromosome]:
                    Start = int(feature[0])
                    End = int(feature[1])
                    if Start < int(position[1]) < End:
                        # Associates product with sequence. Extra indexing b/c parsing from GenBankFile returns list
                        self.feature_targets[sequence[0]].append(feature[3][0])
                        print "found a target"
                        # Adding to Gene Information Container
                        geneid = feature[2][0]
                        if geneid in self.genes:  # have we seen this gene before?
                            if realseq in self.genes[geneid][0]:  # have we seen this sequence on this gene before?
                                self.genes[geneid][0][realseq] += 1
                            else:
                                self.genes[geneid][0][realseq] = 1
                            self.genes[geneid][1] += 1  # adding to total times targeted
                        else:
                            self.genes[geneid] = [{realseq: 1}, 1]
        print "function ended."
        return self.feature_targets

    def write_data_files(self):
        newpath = "/Users/brianmendoza/Desktop/CrisprDB/multiData/" + self.org + "/"
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        os.chdir(newpath)

        # Writing Regular Data File
        filename = self.org + "data.txt"
        f = open(filename, "w")
        for element in self.feature_targets:
            seqidstring = element + "," + str(len(self.feature_targets[element]))
            for item in self.feature_targets[element]:
                myitem = "," + item
                seqidstring += myitem
            f.write(seqidstring + "\n")
        f.close()
        # Writing Gene Data File
        filename = self.org + "genedata.txt"
        f = open(filename, "w")
        gic = self.get_genes()
        for gene in gic:
            totstring = gene[0] + ";" + str(gene[1]) + ";" + str(gene[2])
            f.write(totstring + "\n")
        f.close()
        # Writing All Sequence Data File
        filename = self.org + "allseqdata.txt"
        f = open(filename, "w")
        for seq in self.repeats:
            inputstring = seq[0] + ";" + str(seq[2]) + ";"
            for location in seq[1]:
                l = ";" + str(location[0]) + ";" + str(location[1])
                inputstring += l
            f.write(inputstring + "\n")
        f.close()


    # Function returns the genes and the times each is targeted in descending order of times  each is targeted
    def get_genes(self):
        sorted_genes = []
        for gene in self.genes:
            tuple = (gene, self.genes[gene][1], len(self.genes[gene][0]))
            sorted_genes.append(tuple)
        actual_sorted = sorted(sorted_genes, key=operator.itemgetter(1), reverse=True)
        return actual_sorted

#orgs = ["cac", "cace", "cbe", "cbut", "ccb", "ckl", "csb", "csr", "cth", "ctx",
        #"dha", "med", "tro", "ttm", "tto"]
#for org in orgs:
#a = AnnotationTable("cce")

        # if S2 < S1 < E2 or S1 < S2 < E1 then they overlap



