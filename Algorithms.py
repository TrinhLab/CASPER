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
                rem = number % 4
                number = int(number/4)
                seq += self.int2nt(rem)
            seq += self.int2nt(number)
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
            seed = seq[:split]
            tail = seq[split+1:]
        else:
            seq = seq.split("-")
            seed = seq[0]
            tail = seq[1]
        seed = self.decompress64(seed, True)
        tail = self.decompress64(tail, True)
        # The for loops fixes the problem of A's not being added to the end because they are removed on compression
        for i in range(len(seed),16):
            seed += 'A'
        for j in range(len(tail),4):
            tail += 'A'
        myseq = tail + seed
        return loc, myseq, scr
        #print("Location: " + str(myloc))
        #print("Sequence: " + myseq)

    def fill_As(self, seq, length):
        newseq = seq
        for i in range(len(seq),length):
            newseq += 'A'
        return newseq


#S = SeqTranslate()
#S.decompress_csf_tuple("RgK,BACMo0-w")
#print(S.decompress64("BALKN",True))

