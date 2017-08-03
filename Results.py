import sys
from PyQt5 import Qt, QtWidgets, uic
from Scoring import OnTargetScore

# =========================================================================================
# CLASS NAME: Results
# Inputs: Takes information from the main application window and displays the gRNA results
# Outputs: The display of the gRNA results search
# =========================================================================================


class Results(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super(Results, self).__init__(parent)
        uic.loadUi('resultsWindow.ui', self)
        self.setWindowTitle('Results')
        self.geneViewer.setReadOnly(True)
        # Scoring Class object #
        self.onscore = OnTargetScore()
        # Data containers #
        self.allGenes = []
        self.allTargets = {}
        self.allGeneSeqs = {}
        self.startpos = 0
        self.endpos = 0
        # Target Table settings #
        self.targetTable.setColumnCount(5)  # hardcoded because there will always be five columns
        self.targetTable.setShowGrid(False)
        self.targetTable.setHorizontalHeaderLabels("Sequence;PAM;Strand;Score;Off Targets".split(";"))
        self.targetTable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.targetTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.targetTable.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

    def loadGenesandTargets(self, genesequence, start, end, targets, genename):
        print(genesequence)
        print(targets)
        print(genename)
        self.startpos = start
        self.endpos = end
        self.allTargets[genename] = targets
        self.allGeneSeqs[genename] = genesequence
        self.comboBoxGene.addItem(genename)
        self.displayGeneData()
        self.comboBoxGene.currentIndexChanged.connect(self.displayGeneData)

    def displayGeneData(self):
        curgene = str(self.comboBoxGene.currentText())
        cg = self.allGeneSeqs[curgene]
        self.geneViewer.setPlainText(cg)
        #  --- Shifting numbers over based on start and end ---  #

        self.targetTable.setRowCount(len(self.allTargets[curgene]))
        print(self.allTargets[curgene])
        index = 0
        for item in self.allTargets[curgene]:
            st = item[0]-self.startpos
            print(st)
            seq = QtWidgets.QTableWidgetItem(cg[st-20:st])
            self.targetTable.setItem(index, 0, seq)
            pam = QtWidgets.QTableWidgetItem(cg[st:st+3])  # this is only for a 3nt PAM need to import pam info
            self.targetTable.setItem(index, 1, pam)
            strand = QtWidgets.QTableWidgetItem(item[1])
            self.targetTable.setItem(index, 2, strand)
            scr = self.onscore.returnScore(cg[st-6:st+30])
            print(scr)
            score = QtWidgets.QTableWidgetItem(str(scr))
            self.targetTable.setItem(index, 3, score)
            self.btn_sell = QtWidgets.QPushButton('Find Off Targets')
            self.btn_sell.clicked.connect(self.handleButtonClicked)
            self.targetTable.setCellWidget(index, 4, self.btn_sell)
            index += 1

        self.targetTable.resizeColumnsToContents()

    def handleButtonClicked(self):
        # button = QtGui.qApp.focusWidget()
        button = self.sender()
        index = self.targetTable.indexAt(button.pos())
        if index.isValid():
            print(index.row(), index.column())


# ----------------------------------------------------------------------------------------------------- #