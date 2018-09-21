#make this the file with all the tedious KEGG and Uniprot parsing code

from bioservices import KEGG
import requests
from bs4 import BeautifulSoup
from Algorithms import SeqTranslate


class Kegg:

    k = KEGG()
    location = "http://www.genome.jp/dbget-bin/www_bget?"

    def gene_locator(self, gene_id):
        res = self.k.get(gene_id)
        d = self.k.parse(res)
        newstr = d['POSITION']
        cstop = newstr.find(':')
        if cstop == -1:
            chromosome = 1
        else:
            chrom = newstr[0:cstop]
            chromosome = self.translate_chromosome(chrom)
        sense = True
        if newstr.find('complement') != -1:  # it is on the opposite strand of DNA
            sense = False
            cstop = newstr.find('(')
        if newstr.find('join') != -1:
            cstop = newstr.find('(')
            srt = cstop + 1
            bothpos = newstr[srt:len(newstr)-1].split(",")
            totpos = []
            for i in range(0, len(bothpos)):
                spc = bothpos[i].find('..')
                spos = bothpos[i][0:spc]
                epos = bothpos[i][spc+2:len(bothpos[i])]
                totpos.append((spos, epos))
            startloc = totpos[0][0]
            endloc = totpos[len(bothpos)-1][1]
        else:
            srt = cstop + 1
            spc = newstr.find('..')
            startloc = newstr[srt:spc]
            if not sense:
                endloc = newstr[spc+2:len(newstr)-1]
            else:
                endloc = newstr[spc+2:len(newstr)]
        totloc = (chromosome, sense, int(startloc), int(endloc))
        return totloc

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
        return 1

    def revcom(self, sequence):
        revseq = ""
        change = {'A':'T',
                  'T':'A',
                  'G':'C',
                  'C':'G'}
        for nt in sequence:
            rnt = change[nt]
            revseq = rnt + revseq
        return revseq

    def added_nts(self, seqstart, seqend, vector, orgcode, chromosome):
        url = self.location + "FROM=" + seqstart + "&TO=" + seqend + "&VECTOR="\
              + vector + "&ORG=" + orgcode + "&CHR=" + chromosome
        source_code = requests.get(url)
        plain_text = source_code.text
        soup = BeautifulSoup(plain_text, "html.parser")
        exons = []
        x = soup.find('pre')
        for region in soup.findAll('font'):
            seq = str(region)
            st = seq.find('>') + 1
            z = seq.find('/font') - 2
            exons.append(seq[st:z])
        # get the sequence:
        sx = str(x)
        start = sx.find("\n", 10) + 1
        end = sx.find("/pre") - 2
        unfontseq = sx[start:end]
        trueseq = ""
        ingene = True
        i = 0
        for nt in unfontseq:
            if nt == "<":
                ingene = False
            if nt == ">":
                ingene = True
            elif ingene:
                trueseq += nt
            i += 1
        exon_position_tuples = []
        for exon in exons:
            spos = trueseq.find(exon)
            epos = spos + len(exon)
            pos = (spos, epos)
            exon_position_tuples.append(pos)
        return exon_position_tuples


class SeqFromCSPR:

    def __init__(self, filename):
        self.filename = filename
        self.targets = list()
        self.repeats = dict()

    def decode_targets(self):
        f = open(self.filename)
        # make sure to recognize chromosome number
        data = f.readline()[:-1]
        while data != "REPEATS":
            data = f.readline()[:-1]
            # parse location and sequence
            midpoint = data.find(',')
            location = data[:midpoint]
            sequence = data[midpoint+1:]
            # decompress the location and sequence information
            s = SeqTranslate()
            location = s.decompress64(location,toseq=False)
            sequence = s.decompress64(sequence,toseq=True)
            # add location to storage vector
            self.targets.append((location, sequence))

    def decode_repeats(self):
        f = open(self.filename)
        for line in f:
            if line[:-1] == "REPEATS":
                break
        r_line = f.readline()

    def print_targets(self):
        for item in self.targets:
            print(item)

class SeqFromFasta:

    def __init__(self):
        self.filename = ""
        self.genesequence = ""
        self.targets = []

    def getgenesequence(self):
        return self.genesequence

    def gettargets(self):
        return self.targets

    def complement(self, myseq):
        revseq = ""
        change = {'A': 'T',
                  'T': 'A',
                  'G': 'C',
                  'C': 'G'}
        for nt in myseq:
            rnt = change[nt]
            revseq = rnt + revseq
        return revseq

    def setfilename(self, filename):
        self.filename = filename

    def getsequenceandtargets(self, geneinfo, added_front, added_end, dbpath, endo):
        chrom = int(geneinfo[0])
        print(chrom)
        chromseq = ""
        if geneinfo[1]:
            spos = geneinfo[2] - added_front
            epos = geneinfo[3] + added_end
        else:
            spos = geneinfo[2] - added_end
            epos = geneinfo[3] + added_front
        f = open(self.filename)
        counter = 0
        for line in f:
            if line[0] == '>':
                counter += 1
                continue
            if counter == chrom:
                chromseq += line[0:-1]
            elif counter > chrom:
                break
        f.close()
        sequence = chromseq[spos-1:epos]
        if geneinfo[1]:
            self.genesequence = sequence
        else:
            self.genesequence = self.complement(sequence)

        # --- Target acquisition now ---- #

        filename = dbpath + '-' + endo + ".txt"
        f = open(filename)
        while True:
            line = f.readline()
            if line[0:-1] == "CHROMOSOME #" + str(chrom):
                while True:
                    pos = f.readline()
                    direct = pos[-2:-1]
                    pos = int(pos[0:-2])
                    if geneinfo[2] < pos:
                        targetpos = (pos, direct)
                        self.targets.append(targetpos)
                    if geneinfo[3] < pos:
                        break
                break
        f.close()


# https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:T00005
# https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:T01033