__author__ = 'brianmendoza'

class ASummary:

    def __init__(self):
        self.orgcode = "eco"
        self.ids = []
        self.rawDefs = []
        self.catCount = {}

        self.get_categories()
        print self.ids
        self.add_categories()

    def get_categories(self):
        f = open("/Users/brianmendoza/Desktop/protein_tables/protein_families.txt", 'r')
        while True:
            line = f.readline()
            if line.find("OrgansimCode:" + self.orgcode):
                newfamily = f.readline()
                if newfamily == "\n":
                    break
                else:
                    self.ids.append(newfamily)
        f.close()

    def add_categories(self):
        for id in self.ids:
            self.catCount[id] = 0
        fname = "/Users/brianmendoza/Desktop/CRISPRs/" + self.orgcode + "_multi_data.txt"
        f = open(fname, 'r')
        for line in f:
            cat = line.split(";")
            cat = cat[2:]
            for element in cat:
                self.rawDefs.append(element)
        f.close()
        # go through all ids and see if one matches the definition
        for definition in self.rawDefs:
            for id in self.ids:
                if definition.find(id):
                    self.catCount[id] += 1

    def print_results(self):
        for item in self.ids:
            print item + ";" + str(self.catCount[item])


a = ASummary()
a.print_results()

        