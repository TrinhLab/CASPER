__author__ = 'brianmendoza'

import os

from GenBankParse import GenBankFile


class Compare_Orgs:

    def __init__(self, endo, orgs):
        # Dict of dicts. Key: Organism. Value: Dictionary of lists of unique Positions where key is chromosome
        self.positions = {}
        # Dict of lists. Key: Organism. Value: List of repeated Sequences
        self.sequences = {}

        # Dictionary of dictionaries. Key1: organism Key2: sequence Value: position
        self.searchableseqs = {}

        # Container for all the containers
        self.containers = {}

        # Intitialize containers with organism keys
        self.organisms = []
        for org in orgs:
            self.organisms.append(org)
            self.positions[org] = {}
            self.sequences[org] = []
            self.searchableseqs[org] = {}
        self.endo = endo

        # Generates the sequence lists
        for org in self.organisms:
            self.make_lists(org)
            self.find_sequences(org)

        # Runs the comparison
        self.compare()
        self.write_to_file()

    # Takes an organism and parses the target data into positions and repeated sequences containers
    def make_lists(self, org):
        name1 = "/Users/brianmendoza/Desktop/CrisprDB/" + org + self.endo + ".txt"
        f = open(name1, 'r')
        while True:
            position = f.readline()
            if position[0:3] == "BAD":
                break
            # adds to the positions container the unique position with organism, and chromosome as keys
            p = position.split(",")
            if p[0] in self.positions[org]:
                self.positions[org][p[0]].append(position)
            else:
                self.positions[org][p[0]] = []
                self.positions[org][p[0]].append(position)

        while True:
            compseq = f.readline()
            if compseq[0:3] == "NAG":
                break
            if not "," in compseq:
                self.sequences[org].append(self.decompress(compseq[0:-1]))  # -1 to get rid of the "\n"
        f.close()

    def decompress(self, compressed_seq):
        matrixString = '!"#$%&()*+j-./23456789:;<=>?@BDEFHIJKLMNOPQRSkVWXYZ[]^_abcdefghi'
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

    # This function takes in all the organisms compares them to each other generating a data set that includes
    def find_sequences(self, org):
        # Iterates through unique positions to add positions to seqsby container
        # Element corresponds to a list of positions in the element's chromosome
        for element in self.positions[org]:
            gparse = GenBankFile(org)
            chromstring = gparse.getChromSequence(int(element))
            # Searches the chromstring for the position and returns the sequence and position as a tuple
            for item in self.positions[org][element]:
                pos = item.split(",")
                gseq = ''
                begin = int(pos[1])
                if pos[2][0] == 'f':
                    gseq = chromstring[begin-21:begin]  # 21 is for seed length
                if pos[2][0] == 'r':
                    gseq = self.revcom(chromstring[begin:(int(pos[1])+20)])
                seqtuple = (pos[0], pos[1], pos[2][0])  # pos[3] is the score
                self.searchableseqs[org][gseq] = seqtuple
            # Iterates through the repeated sequences and assigns its position as an "error" string
            for item in self.sequences[org]:
                self.searchableseqs[org][item] = ("repeat. ", "look up for position")
            print "done."

    def compare(self):
        a, b = self.organisms[0], self.organisms[1]
        if len(self.organisms) == 3:
            c = self.organisms[2]
        for seq in self.searchableseqs[a]:
            if seq in self.searchableseqs[b]:
                print "found one"
                self.add_to_container(seq,a,b,'0')
                if len(self.organisms) == 3:  # if there are only two to be compared
                    if seq in self.searchableseqs[c]:
                        self.add_to_container(seq,a,b,c)
        if len(self.organisms) == 3:
            for seq in self.searchableseqs[c]:
                if seq in self.searchableseqs[a]:
                    self.add_to_container(seq,a,c,'0')
                elif seq in self.searchableseqs[b]:
                    self.add_to_container(seq,b,c,'0')

    def add_to_container(self, seq, a, b, c):
        if a+b not in self.containers:
            self.containers[a+b] = []
        if a+b+c not in self.containers:
            self.containers[a+b+c] = []
        apos = self.searchableseqs[a][seq]
        bpos = self.searchableseqs[b][seq]
        inputstring = seq + ";" + apos[0] + ";" + apos[1] + ";" + bpos[0] + ";" + bpos[1]
        if c != '0':
            cpos = self.searchableseqs[c][seq]
            inputstring += cpos[0] + ";" + cpos[1]
            self.containers[a+b+c].append(inputstring)
            return
        self.containers[a+b].append(inputstring)

    def write_to_file(self):
        os.chdir("/Users/brianmendoza/Desktop/CrisprDB/multiData/")
        filename = "compare_"
        for org in self.organisms:
            filename += org + "_"
        filename += ".txt"
        f = open(filename, 'w')
        for container in self.containers:
            original = 0
            repeats = 0
            f.write(str(container) + "\n")
            for item in self.containers[container]:
                if item.find("repeat") != -1:
                    repeats += 1
                else:
                    original += 1
                whatin = str(item)
                f.write(whatin + "\n")
            f.write("Statistics: " + str(original) + ";" + str(repeats) + "\n")
            f.write("NEXT" + "\n")
        f.close()

#clo = ["ctx", "cac", "cace", "cbe", "cbut", "ccb", "ckl", "csb", "csr", "ttm", "tto", "ebl", "ate", "tre", "cce"]
#for cac in clo:
pare = Compare_Orgs("SpCas9", ["ckl", "cac", "cce"])

