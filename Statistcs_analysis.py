"""Get data of .cspr and .gbff files"""


def quick_stats():
    file = open("/Users/brianmendoza/Dropbox/JGI_CASPER/kfdspCas9.cspr")

    Stats = {"Unique": 0, "Repeats": 0,}

    Repeats = list()

    repeats = False
    chrom = 0

    while True:
        line = file.readline()
        if line.startswith("END"):
            break
        if line.startswith(">"):
            chrom += 1
        elif line.startswith("REPEAT"):
            repeats = True
        # Repeat section reading
        elif repeats:
            tag = line
            tot_hits = len(file.readline().split("\t"))
            Repeats.append((tag,tot_hits))
            Stats["Repeats"] += 1
        # regular line reading
        else:
            Stats["Unique"] += 1
    file.close()

    print()
    #print(sorted(Repeats,key=lambda sgrna: sgrna[1],reverse=True))


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
    for line in f:
        if line.startswith(">"):
            curscaff = line[1:-1]
            if scaffold == curscaff:
                in_scaff = True
            


gbff_stats("/Users/brianmendoza/Dropbox/JGI_CASPER/Kfedtschenkoi_gene_exons.gff3")