"""Algorithms.py File.  This file contains the following classes: SeqTranslate.
    SeqTranslate Class. Used for interpreting base64 representations of the target locations as well as their sequences.
    To interpret these run the class instance at the bottom of the file with the desired base64 representation into the
    decompress_tuple function."""


class SeqTranslate:

    def __init__(self):
        self.base_array_64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

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

    # Decompresses the base64 representation into base10.  If toseq is true it returns the sequence itself (nucleotides)
    def decompress64(self, base64seq, toseq=False):
        base10seq = int()
        for i in range(len(base64seq)):
            power = len(base64seq) - (i+1)
            index = self.base_array_64.find(base64seq[i])
            base10seq += index*pow(64, power)
        if toseq:
            seq = str()
            number = base10seq
            while number>=4:
                rem = number%4
                number = int(number/4)
                seq += self.int2nt(rem)
            seq += self.int2nt(number)
            return seq
        else:
            return base10seq

    def decompress_tuple(self,locseq):
        mytuple = locseq.split(",")
        loc = mytuple[0]
        seq = mytuple[1][1:-1]
        myloc = self.decompress64(loc)
        myseq = self.decompress64(seq,True)
        print("Location: " + str(myloc))
        print("Sequence: " + myseq)

S = SeqTranslate()
S.decompress_tuple("Q1k,+wDJcL7M")
