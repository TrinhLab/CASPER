import sys
import math


class SeqTranslate:

    def __init__(self):
        self.base_array_64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

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

    def decompress64(self, base64seq, toseq=True):
        base10seq = int()
        for i in range(len(base64seq)):
            power = len(base64seq) - (i+1)
            index = self.base_array_64.find(base64seq[i])
            base10seq += index*pow(64, power)
        if toseq:
            seq = str()
            number = base10seq
            while number>4:
                rem = number%4
                number = int(number/4)
                seq += self.int2nt(rem)
            seq += self.int2nt(number)
            return seq
        else:
            return base10seq
