import sys, os
import subprocess as sub
from threading import Thread
from queue import Queue, Empty
from PyQt5 import Qt, QtWidgets, uic
from bioservices import KEGG


def iter_except(function, exception):
    """Works like builtin 2-argument `iter()`, but stops on `exception`."""
    try:
        while True:
            yield function()
    except exception:
        return


class NewGenome(QtWidgets.QMainWindow):
    def __init__(self):
        super(NewGenome, self).__init__()
        uic.loadUi('NewGenome.ui', self)

        self.k = KEGG()

        #---Button Modifications---#
        self.setWindowIcon(Qt.QIcon("cas9image.png"))
        self.whatsthisButton.clicked.connect(self.whatsthisclicked)
        self.lineEdit_1.textEdited(self.updatekegglist)
        self.keggsearchresults.setReadOnly(True)


        #show functionalities on window
        self.show()

    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def findFasta(self):
        choice = QtWidgets.QMessageBox.question(self, "Extract!", "Are you sure you want to Quit?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            sys.exit()
        else:
            pass

    def testcheckandradio(self, state):
        if state == QtCore.Qt.Checked:
            pass

    def whatsthisclicked(self):
        QtWidgets.QMessageBox.information(self, "Organism Code", "The organism code is the manner in which CASPER will"
                                                                 "label it's data files and references for the organism"
                                                                 "you are importing here. It is HIGHLY RECOMMENDED that"
                                                                 "you use the 3-4 letter code used by KEGG as this will"
                                                                 "aid in automatic accession of annotations from the"
                                                                 "database.", QtWidgets.QMessageBox.Ok)

    def updatekegglist(self):
        self.keggsearchresults.clear()
        kegg_orglist = self.k.lookfor_organism(self.lineEdit_1.text())
        for item in kegg_orglist:
            item += "\n"
            self.keggsearchresults.insertPlainText(item)

    def run_job(self):
        newdatafile = self.lineEdit_3.text() + self.comboBoxEndo.currentText() + ".txt"
        newdatafile = str(os.curdir) + newdatafile
        f = open(newdatafile, 'w')
        # Write the organism name and subspecies/strain if applicable to the top of the file
        f.write(self.lineEdit_1.text() + self.lineEdit_2.text() + self.comboBoxEndo.currentText())

    def reset(self):
        self.lineEdit_1.clear()
        self.lineEdit_2.clear()
        self.lineEdit_2.clear()
        self.keggsearchresults.clear()



class DisplaySubprocessOutput:
    def __init__(self, window, outlabel, passinfo):  # location of .exe and infopath of the cfcode file with information
        self.root = window
        self.outlabel = outlabel

        # start dummy subprocess to generate some output
        self.process = sub.Popen(passinfo, stdout=sub.PIPE)

        # launch thread to read the subprocess output
        #   (put the subprocess output into the queue in a background thread,
        #    get output from the queue in the GUI thread.
        #    Output chain: process.readline -> queue -> label)
        q = Queue(maxsize=1024)  # limit output buffering (may stall subprocess)
        t = Thread(target=self.reader_thread, args=[q])
        t.daemon = True  # close pipe if GUI process exits
        t.start()

        # show subprocess' stdout in GUI
        self.label = tk.Label(self.root, text="  ", font=(None, 12))
        self.label.pack(ipadx=4, padx=4, ipady=4, pady=4, fill='both')
        self.clabel = tk.Label(self.root, text="  ", font=(None, 12))
        self.clabel.pack(ipadx=4, padx=4, ipady=4, pady=4, fill='both')
        self.update(q) # start update loop

    def reader_thread(self, q):
        """Read subprocess output and put it into the queue."""
        try:
            with self.process.stdout as pipe:
                for line in iter(pipe.readline, b''):
                    q.put(line)
        finally:
            q.put(None)

    def update(self, q):
        """Update GUI with items from the queue."""
        for line in iter_except(q.get_nowait, Empty):  # display all content
            if line is None:
                self.outlabel['text'] = "Program Complete."
                self.quit()
                return
            else:
                self.label['text'] = line # update GUI
                if "Chromosome" in line:
                    self.clabel['text'] = line
                if "Time" in line:
                    tk.Label(self.root, text=line, font=(None, 12)).pack(ipadx=2, fill="both")
                if "There" in line:
                    tk.Label(self.root, text=line, font=(None, 12)).pack(ipadx=2, fill="both")
                break  # display no more than one line per 40 milliseconds
        self.root.after(40, self.update, q)  # schedule next update

    def quit(self):
        self.process.kill()  # exit subprocess if GUI is closed (zombie!)