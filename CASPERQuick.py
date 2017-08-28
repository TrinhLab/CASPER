"""This file is a rapid way to investigate CASPER_Seq_Finder generated files for target locations. Run this python file
    AFTER you have created your setup file: Cquicksetup.txt."""

setup_path = "Cquicksetup.txt"  # Change this to the complete path of your setup file


# ------------------------------------------- CODE THE USER CAN IGNORE --------------------------------------------- #
import os
from Algorithms import SeqTranslate
from bioservices import KEGG
import requests
from bs4 import BeautifulSoup


class CasperQuick:

    def __init__(self):
        self.ST = SeqTranslate()
        self.allTargets = {}
        self.location = tuple()
        self.output = str()
        self.loadGenesandTargets()
        self.printoutresultstofile()

    def loadGenesandTargets(self):
        f = open(setup_path)
        info = dict()
        for line in f:
            if line.find('=') != -1:
                mytuple = line.split('=')
                info[mytuple[0]] = mytuple[1][:-1]
        f.close()
        casperseqfile = info['FILENAME']
        self.output = info['OUTPUT FILE PATH']
        region_keggs = info['REGION OR KEGG CODE'].split(';')
        for region_kegg in region_keggs:
            if region_kegg[0] == '(':
                region = region_kegg[1:-1].split(',')
                self.location = region
            else:
                k = Kegg()
                self.location = k.gene_locator(region_kegg)
            myfy = open(casperseqfile)
            while True:
                line = myfy.readline()
                if line.find('CHROMOSOME') != -1:
                    if line[line.find('CHROMSOME')+2:-1] == self.location[0]:
                        break
            curpos = int()
            while curpos < self.location[1]:
                line = myfy.readline()
                curpos = self.ST.decompress64(line.split(',')[0])
            while curpos < self.location[2]:
                line = myfy.readline().split(",")[:-1]
                curpos = self.ST.decompress64(line[0])
                seq = line[1][0] + self.ST.decompress64(line[1][1:-1], True)
                mytuple = (curpos,seq)
                self.allTargets[region_kegg].append(mytuple)
            myfy.close()

    def printoutresultstofile(self):
        f = open(self.output, 'w')
        for item, targets in self.allTargets:
            f.write(item)
            for target in targets:
                f.write(target)
        f.close()


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
        totloc = (chromosome, int(startloc), int(endloc))
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