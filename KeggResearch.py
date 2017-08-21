__author__ = 'brianmendoza'

from bioservices import KEGG

from Reading import AnnotationTable


class KeggInfo:

    def __init__(self):
        self.k = KEGG()
        self.org = "lac"
        self.genelist = []
        self.genedict = {}
        self.a = AnnotationTable()
        self.targetdict = self.a.analyze_sequences()  # listed as sequence: list of genes
        for gene in self.a.genes:
            self.genelist.append(gene)
        #for gene in self.genelist:
            #self.get_info(gene)

    def get_info(self, gene):
        id = self.org + ":" + gene
        res = self.k.get(id)
        d = self.k.parse(res)
        ortho = "unknown"
        motif = "unknown"
        pfam = "unknown"
        definition = str(d['DEFINITION'])
        definition = definition[9:]
        if d.has_key('ORTHOLOGY'):
            ortho = str(d['ORTHOLOGY'])
        if d.has_key('MOTIF'):
            motif = d['MOTIF']
            if motif.has_key('Pfam'):
                pfam = str(motif['Pfam'])
            else:
                pfam = "unknown"
        # print gene + ";" + definition + ";" + pfam
        self.genedict[gene] = definition
        print gene + " info obtained"

    def make_file(self):
        f = open("/Users/brianmendoza/Desktop/CRISPRs/lac_multi_data.txt", 'w')
        for sequence in self.targetdict:
            sequenceLine = sequence + ";" + str(len(self.targetdict[sequence]))
            for gene in self.targetdict[sequence]:
                sequenceLine += ";" + gene[0:-2]  # self.genedict[gene]
            f.write(sequenceLine + "\n")
        f.close()
