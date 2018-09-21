"""This file will be merged with Statistics analysis into one class.  This file gives iSTOP stats."""

import sys

#  Arguments for running Stat2.py:
concise_file = sys.argv[1]
output_file = sys.argv[2]
is_gene = str(sys.argv[3]).lower()  # yes Yes or YES


def run_by_exon():
    gene_istops = dict()

    f = open(concise_file)
    # Info to keep track of the gene identification and exon:
    cur_gene = ""
    ginfo = list()
    exon_count = 1
    istops = list()
    for line in f:
        # skip scaffold lines:
        if line.startswith("SCAFF"):
            poo = 1
        elif line.startswith("GENE"):
            # consolidate the old info first:
            gene_exon = cur_gene + ",e" + str(exon_count)
            # Sort the istops list:
            gene_istops[gene_exon] = sorted(istops)
            istops = []

            # Find the new gene and its info:
            ginfo = line[line.find(":")+2:-1].split(",")  # +1 to deal with space

            # Keeps track of what gene you are in
            if ginfo[0] == cur_gene:
                exon_count += 1
            else:
                cur_gene = ginfo[0]
                exon_count = 1
        # Look into the line to see if an iSTOP occurs:
        elif line.find("YES") != -1:
            ihit = line.split(",")
            rel_loc = (int(ihit[0]) - int(ginfo[1]))/(int(ginfo[2]) - int(ginfo[1]))
            istops.append(rel_loc)
    f.close()

    # report in file:
    out = open(output_file, "w")
    for gene in gene_istops:
        if len(gene_istops[gene]) > 0:
            out.write(gene + "," + str(gene_istops[gene][0]) + "\n")
        else:
            out.write(gene + ",1\n")
    out.close()


def run_by_gene():
    gene_istops = dict()

    f = open(concise_file)
    # Info to keep track of the gene identification and exon:
    cur_gene = ""
    ginfo = list()
    past_gene_length = 0
    tot_gene_length = 0
    istops = list()
    for line in f:
        # skip scaffold lines:
        if line.startswith("SCAFF"):
            poo = 1
        elif line.startswith("GENE"):
            # Find the new gene and its info:
            ginfo = line[line.find(":") + 2:-1].split(",")  # +1 to deal with space
            # Keeps track of what gene you are in
            if ginfo[0] == cur_gene:
                # adapt the length
                past_gene_length += (int(ginfo[2]) - int(ginfo[1]))
                tot_gene_length += past_gene_length
            else:
                # Load the gene information into the output dict:
                rel_istops = list()
                for loc in istops:
                    rel_istops.append(loc/tot_gene_length)
                gene_istops[cur_gene] = sorted(rel_istops)
                # reset parameters for the new gene:
                istops = []
                cur_gene = ginfo[0]
                past_gene_length = 0
                tot_gene_length = int(ginfo[2])-int(ginfo[1])  # initializes gene length with first exon

        # Look into the line to see if an iSTOP occurs:
        elif line.find("YES") != -1:
            ihit = line.split(",")
            raw_loc = (int(ihit[0]) - int(ginfo[1]) + past_gene_length)
            istops.append(raw_loc)
    f.close()

    # report in file:
    out = open(output_file, "w")
    for gene in gene_istops:
        if len(gene_istops[gene]) > 0:
            out.write(gene + "," + str(gene_istops[gene][0]) + "\n")
        else:
            out.write(gene + ",1\n")
    out.close()


if is_gene == "yes":
    run_by_gene()
else:
    run_by_exon()


