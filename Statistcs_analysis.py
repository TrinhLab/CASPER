"""Get data of .cspr and .gbff files"""

from SeqTranslate import SeqTranslate


def quick_stats(get_repeats=False):
    S = SeqTranslate()
    file = open("/Users/brianmendoza/Desktop/bsuspCas9.cspr")

    Stats = {"Unique": 0, "Repeats": 0,}

    Repeats = list()

    repeats = False
    chrom = 0

    while True:
        line = file.readline()
        if line.startswith("END_OF"):
            break
        if line.startswith(">"):
            chrom += 1
        elif line.startswith("REPEAT"):
            repeats = True
        # Repeat section reading
        elif repeats:
            tag = line[:-1]
            hits = file.readline().split("\t")[:-1]
            Repeats.append((tag,hits))
            Stats["Repeats"] += 1
        # regular line reading
        else:
            Stats["Unique"] += 1
    file.close()

    sortedRepeats = sorted(Repeats, key=lambda sgrna: sgrna[1], reverse=True)
    print(sortedRepeats)
    if get_repeats:
        expandedRepeats = dict()  # contains the information
        for multi in sortedRepeats:
            baseseq = S.decompress64(multi[0], 16, True)
            mykey = multi[0] + "," + str(len(multi[1]))
            # Decompress each tail
            for tail in multi[1]:
                mytail = tail.split(",")
                chromosome = mytail[0]
                location = S.decompress64(mytail[1], False)
                # find the marker of the end:
                if mytail[3].find("+") != -1:
                    end = mytail[3].find("+")
                else:
                    end = mytail[3].find("-")
                totseq = S.decompress64(mytail[3][:end], 4, True) + baseseq
                keystr = chromosome + "," + location + "," + totseq
                expandedRepeats[mykey].append(keystr)
        # Output the repeats into a file:



def gbff_stats(myfile, codetype = "CDS"):
    gene_cds = dict()  # stores gene name and then the lengths of the total cds's
    f = open(myfile)
    for line in f:
        if line.startswith("##"):
            print(line)
        else:
            istics = line[:-1].split("\t")
            # grab only the exon lines
            if istics[2] == codetype:
                keywds = istics[8].split(";")
                gene_exonid = keywds[1][keywds[1].find("=")+1:-5]
                geneid = gene_exonid.split(".")[0]
                # reports the gene identification, positions, and direction
                mistics = [geneid,int(istics[3])-1,int(istics[4]),istics[6]]
                #print(gene_exonid,(mistics[2]-mistics[1])%3)
                if geneid in gene_cds:
                    gene_cds[geneid] += (mistics[2]-mistics[1])
                else:
                    gene_cds[geneid] = (mistics[2]-mistics[1])
    f.close()
    for gene in gene_cds:
        print(gene,gene_cds[gene]%3)


def get_subsequence(orgfile, scaffold, start, stop):
    f = open(orgfile)
    scaffold_seq = str()
    in_scaff = False
    for line in f:
        if line.startswith(">"):
            if in_scaff:
                break
            curscaff = line[1:-1]
            if scaffold == curscaff:
                in_scaff = True
                continue
        elif in_scaff:
            scaffold_seq += line
    f.close()
    # now get the subsequence of the scaffold belonging to the cds:
    print(scaffold_seq[start:stop])


quick_stats(True)
#gbff_stats("/Users/brianmendoza/Dropbox/JGI_CASPER/Kfedtschenkoi_gene_exons.gff3")
#get_subsequence("/Users/brianmendoza/Dropbox/JGI_CASPER/kfd.fna","Scaffold_1", 20713, 22117)