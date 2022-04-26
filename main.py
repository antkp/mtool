import sys
from pyqtgraph.Qt import QtCore, QtGui , QtWidgets
from ui import UI
from data import Data
from control import Control

import ctypes
import platform

#ToDo:
# test the code !!!
# spectrogramm RIO line ?
# maybe use setup.py ! ?
# check number of values in "export" function
# add function "export to Zerodurmerge" ???

try:
    # Include in try/except block if you're also targeting Mac/Linux
    from pyqtgraph.Qt.QtWinExtras import QtWin
    myappid = 'mycompany.myproduct.subproduct.version'
    QtWin.setCurrentProcessExplicitAppUserModelID(myappid)
except ImportError:
    pass

def main():
    print('main')


    def make_dpi_aware():
        if int(platform.release()) >= 8:
            ctypes.windll.shcore.SetProcessDpiAwareness(True)
    # add this code before "app = QtWidgets.QApplication(sys.argv)"
    #make_dpi_aware()   # WINdows !?


    app = QtGui.QApplication([])
    win = QtGui.QMainWindow()
    ui = UI(win)
    app.setWindowIcon(QtGui.QIcon('log.png'))
    data = Data()
    control = Control(data, ui)
    win.show()
    control.update_ui_w2()

    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()


if __name__ == '__main__':
    main()