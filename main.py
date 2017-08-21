import sys
import os
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from APIs import Kegg, SeqFromFasta
from Results import Results

# =========================================================================================
# CLASS NAME: CMainWindow
# Inputs: Takes in the path information from the startup window and also all input parameters
# that define the search for targets e.g. endonuclease, organism genome, gene target etc.
# Outputs: The results of the target search process by generating a new Results window
# =========================================================================================


class CMainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(CMainWindow, self).__init__()
        uic.loadUi('CASPER_main.ui', self)
        self.dbpath = ""
        self.data = {}
        self.orgcodes = {}
        #---Button Modifications---#
        #self.setWindowIcon(QtGui.QIcon("cas9image.png"))
        self.pushButton_FindTargets.clicked.connect(self.gather_settings)
        self.pushButton_ViewTargets.clicked.connect(self.view_results)
        self.pushButton_ViewTargets.setEnabled(False)

        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(100)
        self.progressBar.reset()


        #show functionalities on window
        self.view_my_results = Results()

        self.show()

    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def gather_settings(self):
        inputstring = str(self.geneEntryField.toPlainText())
        self.progressBar.setValue(10)
        if self.radioButton_Gene.isChecked():
            ginput = inputstring.split(',')
            self.run_results("gene", ginput)
        elif self.radioButton_Position.isChecked():
            pinput = inputstring.split(';')
            self.run_results("position", pinput)
        elif self.radioButton_Sequence.isChecked():
            sinput = inputstring
            self.run_results("sequence", sinput)


    # ---- IS ONLY CALLED FROM gather_settings!!!! ---- #
    def run_results(self, inputtype, inputstring):
        kegginfo = Kegg()
        org = self.orgcodes[str(self.orgChoice.currentText())]
        endo = str(self.endoChoice.currentText())
        ginfo = {}  # each entry ginfo[gene] = (chromosome number, t/f strand, start pos, end pos)
        progvalue = 15
        self.progressBar.setValue(progvalue)
        if inputtype == "gene":
            for gene in inputstring:
                g = org + ":" + gene
                ginfo[gene] = kegginfo.gene_locator(g)
                progvalue += 50/len(inputstring)
                self.progressBar.setValue(progvalue)
        if inputtype == "position":
            ginfo = inputstring[1:-1].split(",")
            self.progressBar.setValue(45)
        if inputtype == "sequence":
            self.progressBar.setValue(45)
        s = SeqFromFasta()
        filename = "/Users/brianmendoza/Desktop/GenBank_files/FASTAs/" + org + ".fna"
        s.setfilename(filename)
        progvalue = 75
        self.progressBar.setValue(progvalue)
        for gene in inputstring:
            s.getsequenceandtargets(ginfo[gene], 100, 100, self.dbpath+'/'+org, endo)
            progvalue += 25/len(inputstring)
            self.progressBar.setValue(progvalue)
            self.view_my_results.loadGenesandTargets(s.getgenesequence(), ginfo[gene][2]+100, ginfo[gene][3]-100,
                                                     s.gettargets(), gene)
        self.progressBar.setValue(100)
        self.pushButton_ViewTargets.setEnabled(True)

    def testexe(self):
        choice = QtWidgets.QMessageBox.question(self, "Extract!", "Are you sure you want to Quit?",
                                            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            print(self.orgChoice.currentText())
            sys.exit()
        else:
            pass

    def testcheckandradio(self):
         print(str(self.orgChoice.currentText()))

    # ----- CALLED IN STARTUP WINDOW ------ #
    def getData(self):
        mypath = os.getcwd()
        self.dbpath = mypath
        onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
        print(onlyfiles)
        orgsandendos = {}
        for file in onlyfiles[1:]:
            if file.find('.txt'):
                newname = file[0:-4]
                s = newname.split('-')
                species = str(s[0])
                endo = str(s[1])
                if species in orgsandendos:
                    orgsandendos[species].append(endo)
                else:
                    orgsandendos[species] = [endo]
        self.data = orgsandendos
        os.chdir('/Users/brianmendoza/PycharmProjects/CASPERapp/')
        f = open('CASPERinfo')
        accesstrue = False
        while True:
            line = f.readline()
            if line[0:4] == 'ORGA':
                while True:
                    orginfo = f.readline()
                    orginfo = orginfo[0:-1]
                    if orginfo[0] == '-':
                        break
                    stuff = orginfo.split(":")
                    self.orgcodes[stuff[1]] = stuff[0]
                break
        f.close()

        for item in self.orgcodes:
            self.orgChoice.addItem(item)
        self.endoChoice.addItems(self.data[self.orgcodes[str(self.orgChoice.currentText())]])
        self.orgChoice.currentIndexChanged.connect(self.changeEndos)

    def changeEndos(self):
        self.endoChoice.clear()
        self.endoChoice.addItems(self.data[self.orgcodes[str(self.orgChoice.currentText())]])

    @QtCore.pyqtSlot()
    def view_results(self):
        self.view_my_results.show()
        self.progressBar.setValue(0)
        self.pushButton_ViewTargets.setEnabled(False)

# ----------------------------------------------------------------------------------------------------- #

# =========================================================================================
# CLASS NAME: StartupWindow
# Inputs: Takes information from the main application window and displays the gRNA results
# Outputs: The display of the gRNA results search
# =========================================================================================


class StartupWindow(QtWidgets.QDialog):
    def __init__(self):
        super(StartupWindow, self).__init__()
        uic.loadUi('startupCASPER.ui', self)
        self.setWindowTitle('WELCOME TO CASPER!')

        #---Button Modifications---#
        self.setWindowIcon(QtGui.QIcon("cas9image.png"))
        pixmap = QtGui.QPixmap('mainart.jpg')
        self.labelforart.setPixmap(pixmap)
        self.pushButton_2.setDefault(True)

        self.gdirectory = os.path.expanduser("~")
        self.gdirectory = "/Users/brianmendoza/Desktop/CrisprDB/"
        self.lineEdit.setText(self.gdirectory)

        self.pushButton_3.clicked.connect(self.changeDir)
        self.pushButton_2.clicked.connect(self.show_window)
        self.pushButton.clicked.connect(self.errormsgmulti)

        self.show_main_window = CMainWindow()

        #show functionalities on window
        self.show()

    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def changeDir(self):
        filed = QtWidgets.QFileDialog()
        mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                       self.gdirectory, QtWidgets.QFileDialog.ShowDirsOnly)
        self.lineEdit.setText(mydir)
        cdir = self.lineEdit.text()
        os.chdir(cdir)

    def errormsgmulti(self):
        QtWidgets.QMessageBox.question(self, "Under Construction...", "Sorry this functionality is still"
                                            " under construction and will be available shortly!",
                                            QtWidgets.QMessageBox.Ok)

    @QtCore.pyqtSlot()
    def show_window(self):
        os.chdir(self.gdirectory)
        print(os.getcwd())
        self.show_main_window.getData()
        self.close()


if __name__ == '__main__':
    app = Qt.QApplication(sys.argv)
    startup = StartupWindow()
    sys.exit(app.exec_())
