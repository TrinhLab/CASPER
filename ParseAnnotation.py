"""File does what gbff parse does but makes sure that overlaps are taken care of so no sequences/genes are skipped."""

from SeqTranslate import SeqTranslate
from operator import itemgetter
import os,sys

gbff_file = "/Users/brianmendoza/Dropbox/JGI_CASPER/Kfedtschenkoi_gene_exons.gff3"
cspr_file = "/Users/brianmendoza/Drobpox/JGI_CASPER/kfdspCas9.cspr"
output_csv = "/Users/brianmendoza/Desktop/kfd_concise_output.csv"

grna_temp_storage = list()
# main storage vehicle
ScaffGeneDict = dict()

ScaffRNADict = dict()


def parse_gbff(file_path):
    f = open(file_path)
    for line in f:
        if line.startswith("##"):
            print(line)
        else:
            istics = line[:-1].split("\t")

            # grab only the CDS lines
            if istics[2] == 'CDS':
                keywds = istics[8].split(";")
                geneid = keywds[1][keywds[1].find("=")+1:-5]
                # reports the gene identification, positions, and direction
                # istics[3] is subtracted by 1 to make the indexing 0-base, istics[4] does not need any mod for inc.
                mistics = [geneid,int(istics[3])-1,int(istics[4]),istics[6]]
                # stores the data in a dict of lists for exon/gene and chromosome/scaffold lookup:
                if istics[0] not in ScaffGeneDict:
                    ScaffGeneDict[istics[0]] = [mistics]
                else:
                    ScaffGeneDict[istics[0]].append(mistics)
    f.close()
    # sort by the end location (for parsing index)
    for C_S in ScaffGeneDict:
        ScaffGeneDict[C_S] = sorted(ScaffGeneDict[C_S], key=itemgetter(1))

def parse_cspr(cspr_file):
    grna_temp_storage = list()
    S = SeqTranslate()
    f = open(cspr_file)
    for line in f:
        if line.startswith(">"):
            scaffold = line[1:line.find("(")-1]  # gets only the word scaffold and the number
            print(scaffold)
            if grna_temp_storage:
                if cur_cs in ScaffGeneDict:
                    assign_to_genes(cur_cs)
                grna_temp_storage.clear()
            cur_cs = line[1:line.find("(")-1]
            print("running " + cur_cs)
        elif line.startswith("REPEAT"):
            if cur_cs in ScaffGeneDict:
                assign_to_genes(cur_cs)
            grna_temp_storage.clear()
            break
        else:
            grna_temp_storage.append(S.decompress_csf_tuple(line[:-1]))
    f.close()
