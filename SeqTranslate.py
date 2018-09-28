"""Algorithms.py File.  This file contains the following classes: SeqTranslate.
    SeqTranslate Class. Used for interpreting base64 representations of the target locations as well as their sequences.
    To interpret these run the class instance at the bottom of the file with the desired base64 representation into the
    decompress_tuple function."""


class SeqTranslate:

    def __init__(self):
        # Modification of MIME base64 coding so that +- can be used for strand direction
        self.base_array_64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789=/"


    # used to convert numbers in base4 back to nucleotides
    def int2nt(self, num):
        if num == 0:
            return 'A'
        elif num == 1:
            return 'T'
        elif num == 2:
            return 'C'
        elif num == 3:
            return 'G'
        else:
            return 'N'

    def nt2int(self,nt):
        if nt == 'A':
            return 0
        elif nt == 'T':
            return 1
        elif nt == 'C':
            return 2
        elif nt == 'G':
            return 3
        else:
            return 0

    def compress(self, uncompressed, base):
        compseq = 0
        # Checks to see if the uncompressed is a string or base10 and if it is a nt string,
        # This if statement takes the nucleotides and first turns it into a base10
        if type(uncompressed) == str:
            for i in range(len(uncompressed)):
                val = self.nt2int(uncompressed[i]) * pow(4, i)  # multiplying by power-4 converts to base10
                compseq += val
            uncompressed = compseq
        compreturn = str()
        while uncompressed >= base:
            rem = uncompressed%base
            uncompressed = int(uncompressed/base)
            compreturn = self.base_array_64[rem] + compreturn
        compreturn = self.base_array_64[uncompressed] + compreturn
        return compreturn

    def to_generic_compressed(self, seqobj):
        # Passed as a tuple from the repeats section of the .cspr file
        if type(seqobj) == list:
            gencomp = seqobj[0] + "." + seqobj[1][1:]
        else:
            split = seqobj.find("+")
            if split != -1:
                gencomp = seqobj[:split] + "." + seqobj[split+1:]
            else:
                split = seqobj.find("-")
                gencomp = seqobj[:split] + "." + seqobj[split+1:]
        return gencomp

    # Decompresses the base64 representation into base10.  If toseq is true it returns the sequence itself (nucleotides)
    def decompress64(self, base64seq, seqlength=0, toseq=False):
        base10seq = int()
        for i in range(len(base64seq)):
            power = len(base64seq) - (i+1)
            index = self.base_array_64.find(base64seq[i])
            if index != -1:
                base10seq += index*pow(64, power)
        if toseq:
            seq = str()
            number = base10seq
            while number>=4:
                rem = number % 4
                number = int(number/4)
                seq += self.int2nt(rem)
            seq += self.int2nt(number)
            while len(seq) < seqlength:
                seq = 'A' + seq
            return seq
        else:
            return base10seq

    def decompress_csf_tuple(self, locseq):
        mytuple = locseq.split(",")
        loc = self.decompress64(mytuple[0])
        seq = mytuple[1]
        scr = self.decompress64(mytuple[2])
        split = seq.find("+")
        if split != -1:
            dira = "+"
            sequence = seq[:split]
            pam = seq[split+1:]
        else:
            seq = seq.split("-")
            sequence = seq[0]
            pam = seq[1]
            dira = "-"
        sequence = self.decompress64(sequence, True)
        pam = self.decompress64(pam, True)
        # The for loops fixes the problem of A's not being added to the end because they are removed on compression
        for i in range(len(sequence),20):
            sequence += 'A'
        for j in range(len(pam),3):
            pam += 'A'
        return int(loc), sequence, pam, int(scr), dira
        #print("Location: " + str(myloc))
        #print("Sequence: " + myseq)

    def fill_As(self, seq, length):
        newseq = seq
        for i in range(len(seq),length):
            newseq += 'A'
        return newseq


S = SeqTranslate()
#print(S.decompress_csf_tuple("c,|LLmmwV/-8,t"))

print(S.decompress64("B/6o",toseq=False))
print(S.compress("CACCAAAAGCTACCTCCTGA",64))

