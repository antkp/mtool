import pyqtgraph as pg
import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.dockarea import *
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph.parametertree import Parameter, ParameterTree
import pickle
from pathlib import Path


# todo
#  Pickle Option only for config data not for Filter,Transform ... ???
#  standart-presets to load files
#  Autoscale within each graph


class ScalableGroup(pTypes.GroupParameter):

    def __init__(self, **opts):
        opts['type'] = 'action'
        opts['addText'] = 'add'
        opts['addList'] = ['newband']
        pTypes.GroupParameter.__init__(self, **opts)

    def addNew(self, typ):
        x  = (len(self.childs)/2+1)
        self.addChild(dict(name=" f_low_%d" % x, type='float', value= 1.1, removable=True, renamable=True))
        self.addChild(dict(name="f_high_%d" % x, type='float', value=10.1, removable=True, renamable=True))

class UI:
    def __init__(self, main_window):

        self.fftwindows = [['boxcar', 'blackman', 'hamming', 'hann', 'flattop'],  # windowName
                           [1.0, 2.3833, 1.8534, 2.002, 4.6433],  # amplitudeCorrection
                           [1.0, 1.97, 1.59, 1.63, 2.26]]  # energyCorrection

        self.average_arr = ['- none -', 'simple_moving_average', 'moving p2v', 'moving p2v no slope']
        filter_arr = ['- no filter -', 'gauss', 'butterworth', 'savitzky–golay', 'lowess', 'fft-filter']

        self.params = \
            [
            {'name': 'use last known config', 'type': 'bool', 'value': False},
            {'name': 'select file', 'type': 'action', 'tip': "select file"},
            {'name': 'filepath', 'type': 'str', 'value': '', 'visible': False, 'readonly': True},
            {'name': 'filename', 'type': 'str', 'value': '', 'visible': False, 'readonly': True},
            {'name': 'config data', 'type': 'group', 'expanded': True, 'visible': True, 'children': [
                {'name': 'head rows', 'type': 'int', 'value': 0, 'limits': (0, 20), 'tip': "ignore first rows of data head"},
                {'name': 'delimiter', 'type': 'list', 'values': {'tab':'\t', 'semicolon':';', 'comma':','},
                 'value': '\t', 'tip': 'delimiter between columns'},
                {'name': 'decimal separator', 'type': 'list', 'values': {'comma': ',', 'point':'.'},
                 'value': '.', 'tip': 'decimal separator'},
                {'name': 'x-col', 'type': 'list', 'values': [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
                 'value': 2, 'limits': (-1, 15), 'step': 1, 'tip': 'select row number for X values (first row = 0) \n -1 --> n_th value'},
                {'name': 'y-col', 'type': 'list', 'values': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                 'value': 0, 'limits': (0, 15), 'step': 1, 'tip': "select row number for Y values (first row = 0)"},
                {'name': 'use every n^th row', 'type': 'int', 'value': 1, 'limits': (1, 1000), 'step': 1},
                {'name': 'lengths compensation', 'type': 'bool', 'value': False, 'tip': "for measurements like x(x)"},
                {'name': 'enter', 'type': 'action', 'visible': True}]},
            {'name': 'transformation', 'type': 'group', 'expanded': True, 'visible': True, 'children': [
                {'name': 'plotcolor', 'type': 'color', 'value': '#ffffff', 'tip': 'change plot color'},
                {'name': 'scale x-axis', 'type': 'float', 'value': 1.00, 'tip': "scale coeff x-axis"},
                {'name': 'scale y-axis', 'type': 'float', 'value': 1.00, 'tip': "scale coeff y-axis"},
                {'name': 'x-offset', 'type': 'float', 'value': 0.00, 'tip': 'change offset in x'},
                {'name': 'y-offset', 'type': 'float', 'value': 0.00, 'tip': 'change offset in y'},
                {'name': 'negate', 'type': 'bool', 'value': False, 'tip': 'negate values'},
                {'name': 'reverse', 'type': 'bool', 'value': False, 'tip': 'reverse values --> first value will be the last one'},
                {'name': 'n^th derivative', 'type': 'int', 'value': 0, 'limits': (0, 2), 'step': 1},
                {'name': 'n^th integral', 'type': 'int', 'value': 0, 'limits': (0, 2), 'step': 1},
                {'name': 'linear regression', 'type': 'bool', 'value': False, 'tip': "remove linear regression curve from measurement"},
                {'name': 'first&last point on same level', 'type': 'bool', 'value': False},
                {'name': 'y-zero @ x-zero', 'type': 'bool', 'value': False, 'tip': 'set y offset XYZ'}]},
            {'name': 'filter', 'type': 'group', 'expanded': True, 'visible': True, 'children': [
                {'name': 'select filter', 'type': 'list', 'values':
                    filter_arr, 'value': filter_arr[0]},
                {'name': 'gauss', 'type': 'group', 'expanded': True, 'visible': False, 'tip': 'https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter1d.html','children': [
                    {'name': 'sigma', 'type': 'float', 'value': 3.0, 'limits': (0.0, 6.0), 'step': 0.1}]},
                {'name': 'butterworth', 'type': 'group', 'expanded': True, 'visible': False,
                 'tip': 'https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html',
                 'children': [
                     {'name': 'type', 'type': 'list', 'values':
                         ['lowpass', 'highpass', 'bandpass', 'bandstop'], 'value': 'lowpass'},
                     {'name': 'order', 'type': 'int', 'value': 4, 'limits': (2, 10), 'step': 2}, #warum hier immer nur zweier-Schritte
                     {'name': 'low_cut', 'type': 'int', 'value': 30, 'limits': (1, 10000), 'step': 1},
                     {'name': 'high_cut', 'type': 'int', 'value': 30, 'limits': (1, 10000.0), 'step': 1}]},
                {'name': 'savitzky–golay', 'type': 'group', 'expanded': True, 'visible': False, 'tip': 'https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html', 'children': [
                    {'name': 'size', 'type': 'int', 'value': 31, 'limits': (11, 501), 'step': 1},
                    {'name': 'order', 'type': 'int', 'value': 3, 'limits': (2, 10), 'step': 1}]},
                {'name': 'lowess', 'type': 'group', 'expanded': True, 'visible': False, 'tip': 'https://www.statsmodels.org/stable/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html', 'children': [
                    {'name': 'fraction', 'type': 'float', 'value': 0.025, 'limits': (0.0, 1.0), 'step': 0.01},
                    {'name': 'iteration', 'type': 'int', 'value': 0, 'limits': (0, 5), 'step': 1}]},
                ScalableGroup(name="fft-filter", visible= False, children=[
                ]),
            ]},


            {'name': 'moving section', 'type': 'group', 'expanded': True, 'visible': True, 'tip': 'größte Abweichung (ohne linearen Anteil) innerhalb eines Bereiches', 'children': [
                #{'name': 'show deviation', 'type': 'bool', 'value': False, 'tip': "plot "},
                {'name': 'average', 'type': 'list', 'values': self.average_arr,
                 'value': 'no avarage', 'tip': 'https://en.wikipedia.org/wiki/Moving_average'},
                #{'name': 'exclude linear proportion', 'type': 'bool', 'value': False, 'tip': "include linear proportion within section"},
                {'name': 'section', 'type': 'float', 'value': 100.0}]},
            {'name': 'extend axis', 'type': 'group', 'expanded': True, 'visible': True, 'children': [
                {'name': 'extend data', 'type': 'bool', 'value': False},
                {'name': 'extend', 'type': 'group', 'visible': False, 'children': [
                    {'name': '- fit length', 'type': 'float', 'value': 5.0, 'tip': 'section to fit into'},
                    {'name': '+ fit length', 'type': 'float', 'value': 5.0, 'tip': 'section to fit into'},
                    {'name': 'fit style', 'type': 'list', 'values':
                        ['flat', 'linear', 'quadratic', 'cubic', 'sin', 'exp'], 'value': 'flat'},
                    {'name': 'extend to -Pos', 'type': 'float', 'value': -5.0, 'tip': 'extend to this value'},
                    {'name': 'extend to +Pos', 'type': 'float', 'value': 100.0, 'tip': 'extend to this value'},
                    {'name': 'export for PPMAC', 'type': 'action', 'tip': 'export graph incl. extensions to .csv'}]},
                {'name': 'export to excel', 'type': 'action', 'tip': 'export graph incl. extensions to .csv'},]},
            {'name': 'fft', 'type': 'group', 'expanded': True, 'visible': False, 'children': [
                {'name': 'show fft', 'type': 'bool', 'value': False, 'tip': "This is a Tip"},
                {'name': 'fft config', 'type': 'group', 'expanded': True, 'visible': False, 'children': [
                    {'name': 'amplitude/phase', 'type': 'bool', 'value': False, 'tip': "This is a Tip"},
                    {'name': 'fft window', 'type': 'list', 'values': self.fftwindows[0], 'value': 'boxcar',
                     'tip': 'https://download.ni.com/evaluation/pxi/Understanding%20FFTs%20and%20Windowing.pdf'}]},
                {'name': 'show spectrogram', 'type': 'bool', 'value': False, 'tip': "This is a Tip"},
                {'name': 'spectrogram config', 'type': 'group', 'expanded': True, 'visible': False, 'children': [
                    {'name': 'segment length', 'type': 'int', 'value': 1024, 'limits': (128, 5000000), 'step': 1,
                     'tip': 'https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.spectrogram.html'},
                    {'name': 'fft window', 'type': 'list', 'values': self.fftwindows[0], 'value': 'hann',
                     'tip': 'https://download.ni.com/evaluation/pxi/Understanding%20FFTs%20and%20Windowing.pdf'},
                    {'name': 'overlap', 'type': 'int', 'value': 0, 'limits': (0, 5000), 'step': 1},
                    {'name': 'log_z_scale', 'type': 'bool', 'value': True, 'tip': "scale z-value by (20.*np.log10)"}]}]},
            {'name': 'polar plot', 'type': 'group', 'expanded': True, 'visible': True, 'children': [
                {'name': 'use polar plot', 'type': 'bool', 'value': False},
                {'name': 'polar config', 'type': 'group', 'expanded': True, 'visible': False, 'children': [
                    {'name': 'n', 'type': 'int', 'value': 3, 'tip': "number of revolution"},
                    {'name': 'scale coeff', 'type': 'float', 'value': 3.0, 'tip': "scale the polar plot"},
                    {'name': 'angle offset', 'type': 'float', 'value': 0.0, 'tip': "scale the polar plot"},
                    {'name': 'remove f_fund', 'type': 'bool', 'value': False},
                    {'name': 'show f_fund', 'type': 'bool', 'value': False}]},
]},

            {'name': 'autoscale graph', 'type': 'bool', 'value': True}]

        self.p = Parameter.create(name='params', type='group', children=self.params)

        self.treefile = Path('tree.prm')
        try:
            self.treefile.resolve(strict=True)
        except FileNotFoundError:
            with open('tree.prm', 'wb') as fp:
                pickle.dump(self.p.saveState(), fp)
        else:
            with open('tree.prm', 'rb') as fp:
                 self.p.restoreState(pickle.load(fp))

        self.noneParameter = [{'name': 'none', 'type': 'bool', 'value': False}]
        self.none = Parameter.create(name='none', type='group', children=self.noneParameter)
        self.area = DockArea()
        main_window.setCentralWidget(self.area)
        main_window.resize(2000, 1000)
        main_window.setWindowTitle('mtool_v1.46-ALPHA')
        self.filename = ''
        pg.setConfigOption('background', ('#141414'))
        pg.setConfigOption('foreground', (250, 250, 250))

        self.d1 = Dock("config", size=(300, 700))
        self.d1.hideTitleBar()
        self.d2 = Dock("raw data", size=(900, 100))
        self.d2.hideTitleBar()
        self.d3 = Dock("region", size=(900, 600))
        self.d3.hideTitleBar()
        self.d4 = Dock("fft", size=(900, 0))
        self.d4.hideTitleBar()
        self.d5 = Dock("spectrum", size=(900, 0))
        self.d5.hideTitleBar()

        self.area.addDock(self.d1, 'left')
        self.area.addDock(self.d2, 'right')
        self.area.addDock(self.d3, 'bottom', self.d2)
        self.area.addDock(self.d4, 'bottom', self.d3)
        self.area.addDock(self.d5, 'bottom', self.d4)

        self.t = ParameterTree()
        self.t.setParameters(self.p, showTop=False)

        self.w1 = pg.LayoutWidget()
        self.w1.addWidget(self.t)
        self.d1.addWidget(self.w1)

        self.w2 = pg.PlotWidget(title="-- no data --")
        self.w2.setAutoPan(y=True, x=True)

        self.region = pg.LinearRegionItem(values=(0, 1), brush=(70, 70, 70, 70), movable=True)

        self.wi2Text = ''
        self.wi2 = QtGui.QLabel(self.wi2Text)
        self.wi2.setFont(QtGui.QFont('Menlo', 9, QtGui.QFont.Bold))
        self.wi2.setAlignment(QtCore.Qt.AlignLeft)
        self.wi2.setAlignment(QtCore.Qt.AlignTop)
        self.d2.addWidget(self.w2, row=0, col=0)
        self.d2.addWidget(self.wi2, row=0, col=1)

        self.wi3Text = ''
        self.wi3 = QtGui.QLabel(self.wi3Text)
        self.wi3.setFont(QtGui.QFont('Menlo', 9, QtGui.QFont.Bold))
        self.wi3.setAlignment(QtCore.Qt.AlignLeft)
        self.wi3.setAlignment(QtCore.Qt.AlignTop)
        self.w3 = pg.PlotWidget(title="-- no data --")
        self.w3.setAutoPan(y=True, x=True)
        self.d3.addWidget(self.w3, row=0, col=0)
        self.d3.addWidget(self.wi3, row=0, col=1)
        self.w3vLine = pg.InfiniteLine(angle=90, pen={'color': (100, 100, 100), 'width': 1}, movable=False)
        self.w3hLine = pg.InfiniteLine(angle=0, pen={'color': (100, 100, 100), 'width': 1}, movable=False)

        self.wi4Text = ''
        self.wi4 = QtGui.QLabel(self.wi4Text)
        self.wi4.setFont(QtGui.QFont('Menlo', 9, QtGui.QFont.Bold))
        self.wi4.setAlignment(QtCore.Qt.AlignLeft)
        self.wi4.setAlignment(QtCore.Qt.AlignTop)
        self.w4 = pg.PlotWidget(title="-- no data --")
        self.w4.setAutoPan(y=True, x=True)
        self.d4.addWidget(self.w4, row=0, col=0)
        self.d4.addWidget(self.wi4, row=0, col=1)

        self.w4vLine = pg.InfiniteLine(angle=90, pen={'color': (100, 100, 100), 'width': 1}, movable=False)
        self.w4hLine = pg.InfiniteLine(angle=0, pen={'color': (100, 100, 100), 'width': 1}, movable=False)

        self.w5 = pg.ImageView(view=pg.PlotItem(title="-- no data --", labels={'bottom': ('time', ''), 'left': ('frequency', '')}, xMin=0, yMin=0, lockAspect=True))
        self.d5.addWidget(self.w5, row=0, col=0)
        self.w5vLine = pg.InfiniteLine(angle=90, pen={'color': (100, 100, 100), 'width': 1}, movable=False)
        self.w5hLine = pg.InfiniteLine(angle=0, pen={'color': (100, 100, 100), 'width': 1}, movable=False)

        # self.ra = pg.LineROI([0, 0], [-10, -10], width=1, pen=(1, 9))
        # self.w5.view.addItem(self.ra)
        self.w5.ui.histogram.show()
        self.w5.ui.roiBtn.hide()
        self.w5.ui.menuBtn.hide()

        font = QtGui.QFont()
        font.setPixelSize(12)
        self.plottext = pg.TextItem()
        self.plottext.setFont(font)
        self.plottext.setPos(0.0, 0.9)
