"""Currently the file takes a large genome and breaks it into separate files that can then be quickly interrogated."""


def genome_sub():
    genome = list()  #each item is a string that is a scaffold sequence
    f = open("/Users/brianmendoza/Dropbox/JGI_CASPER/kfd.fna")
    scaffold = str()
    for line in f:
        if line.startswith(">"):
            genome.append(scaffold)
            scaffold = ""
        else:
            scaffold += line
    genome.append(scaffold)
    f.close()

    # return the sequence in question:
    myscaff = 5
    start_end = "475543,476542"
    start = int(start_end.split(",")[0])
    end = int(start_end.split(",")[1])
    print(end-start)
    print(genome[myscaff][start:end])


def concise_info():
    f = open("/Users/brianmendoza/Dropbox/JGI_CASPER/kfd_concise_data.csv")
    count = 0
    for line in f:
        if line.startswith("GENE"):
            myl = line.split(",")
            if int(myl[2]) - int(myl[1]) < 100:
                count += 1
    print(count)
    f.close()

#concise_info()
genome_sub()