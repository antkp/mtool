
import pyqtgraph as pg
from pyqtgraph.dockarea import *
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets

try:
    # Include in try/except block if you're also targeting Mac/Linux
    from pyqtgraph.Qt.QtWinExtras import QtWin
    myappid = 'mycompany.myproduct.subproduct.version'
    QtWin.setCurrentProcessExplicitAppUserModelID(myappid)
except ImportError:
    pass

# ..or..
# import ctypes
# myappid = 'mycompany.myproduct.subproduct.version'
# ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

import sys

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.setWindowTitle("app Name")
        l = QtWidgets.QLabel("Text")
        l.setMargin(10)
        self.setCentralWidget(l)

        self.show()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('logo.ico'))
    w = MainWindow()
    app.exec()