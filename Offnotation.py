"""Right now this file parses the concise_data.csv generated for a whole organism and preps a file for
complete off-target analysis."""

import os

#  This function creates the off_query file
def off_prep():
    lines_for_off = str()

    f = open("/Users/brianmendoza/Dropbox/JGI_CASPER/kfd_concise_data.csv")
    for line in f:
        # Checks for scaffold or gene titles:
        if line.startswith("GENE") or line.startswith("SCAFF"):
            poo = 1
        else:
            lines_for_off += line
    f.close()

    x = open("/Users/brianmendoza/Dropbox/JGI_CASPER/KfdOFF_QUERY.txt",'w')
    x.write(lines_for_off)
    x.close()


#  This function appends the off-target data to the .csv "concise" file
def concise_data_append_off():
    off_target_dict = dict()  # stores all the off targets found by CasperOffTarget
    f = open("/Users/brianmendoza/Dropbox/CASPER/kfd_off_results.txt")
    for line in f:
        l = line[:-1].split(":")
        off_target_dict[l[0]] = float(l[1])
    f.close()
    output_string = str()
    condata = open("/Users/brianmendoza/Dropbox/JGI_CASPER/kfd_concise_data.csv")
    for line in condata:
        # skip scaffold and lines:
        if line.startswith("SCAFF") or line.startswith("GENE"):
            output_string += line
        elif line.split(",")[1] in off_target_dict:
            off_value = 0
            if off_target_dict[line.split(",")[1]] > 0.05:
                off_value = off_target_dict[line.split(",")[1]]
            newline = line[:-1] + "," + str(off_value) + "\n"
            output_string += newline
        else:
            output_string += line
    # now output the revised version:
    outfile = open("/Users/brianmendoza/Dropbox/kfd_concise_data_r.csv", "w")
    outfile.write(output_string)
    outfile.close()




#divide_and_conquer("/Users/brianmendoza/Dropbox/JGI_CASPER/KfdOFF_QUERY.txt","/Users/brianmendoza/Dropbox/kfd_batches/",1000000)
concise_data_append_off()
#off_prep()
