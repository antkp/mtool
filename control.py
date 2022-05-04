import os
import pyqtgraph as pg
import numpy as np
import pickle
from pathlib import Path


class Control:

    def __init__(self, data, ui):
        self.path = ''
        self.title = ''
        self.plotcolor = '#ffffff'
        self.data = data
        self.ui = ui
        self.Sxx = []
        self.spxscale = 1.0
        self.spyscale = 1.0

        self.region = self.ui.region.getRegion()

        self.x_scale = self.ui.p.child('transformation').child('scale x-axis').value()
        self.y_scale = self.ui.p.child('transformation').child('scale y-axis').value()
        self.x_offset = self.ui.p.child('transformation').child('x-offset').value()
        self.y_offset = self.ui.p.child('transformation').child('y-offset').value()
        self.neg = self.ui.p.child('transformation').child('negate').value()
        self.rev = self.ui.p.child('transformation').child('reverse').value()
        self.y_zero = self.ui.p.child('transformation').child('y-zero @ x-zero').value()
        self.derivative = self.ui.p.child('transformation').child('n^th derivative').value()
        self.integral = self.ui.p.child('transformation').child('n^th integral').value()
        self.lreg = self.ui.p.child('transformation').child('linear regression').value()
        self.flpz = self.ui.p.child('transformation').child('first&last point on same level').value()

        self.filter = self.ui.p.child('filter').child('select filter').value()
        self.sigma = self.ui.p.child('filter').child('gauss').child('sigma').value()
        self.size = self.ui.p.child('filter').child('savitzky–golay').child('size').value()
        self.order = self.ui.p.child('filter').child('savitzky–golay').child('order').value()
        self.fraction = self.ui.p.child('filter').child('lowess').child('fraction').value()
        self.iteration = self.ui.p.child('filter').child('lowess').child('iteration').value()

        self.transform_v_arr = ['n^th derivative', 'n^th integral', 'scale x-axis', 'scale y-axis', 'x-offset',
                              'y-offset', 'negate', 'reverse', 'y-zero @ x-zero', 'linear regression',
                              'first&last point on same level']
        self.transform_p_arr = [self.derivative, self.integral, self.x_scale, self.y_scale, self.x_offset,
                                self.y_offset, self.neg, self.rev, self.y_zero, self.lreg, self.flpz]

        self.ui.p.child('filepath').hide()
        self.ui.p.child('filename').hide()
        self.ui.p.child('config data').hide()
        self.ui.p.child('transformation').hide()
        self.ui.p.child('filter').hide()
        self.ui.p.child('moving section').hide()
        self.ui.p.child('fft').hide()
        self.ui.p.child('extend axis').hide()
        self.ui.p.child('polar plot').hide()
        self.ui.p.child('autoscale graph').hide()
        self.ui.p.child('select file').sigActivated.connect(self.on_selctbtn)
        self.data.sig_data_loaded.connect(self.update_ui_w2)

        self.ui.p.child('config data').child('head rows').sigValueChanged.connect(self.on_colselect)
        self.ui.p.child('config data').child('delimiter').sigValueChanged.connect(self.on_colselect)
        self.ui.p.child('config data').child('decimal separator').sigValueChanged.connect(self.on_colselect)
        self.ui.p.child('config data').child('x-col').sigValueChanged.connect(self.on_colselect)
        self.ui.p.child('config data').child('y-col').sigValueChanged.connect(self.on_colselect)
        self.ui.p.child('config data').child('use every n^th row').sigValueChanged.connect(self.on_colselect)

        self.ui.p.child('config data').child('lengths compensation').sigValueChanged.connect(self.on_colselect)
        self.ui.p.child('config data').child('enter').sigActivated.connect(self.on_enter)
        self.ui.region.sigRegionChanged.connect(self.on_region_change)

        self.ui.p.child('transformation').child('plotcolor').sigValueChanged.connect(self.on_plotcolor)
        self.ui.p.child('transformation').child('scale x-axis').sigValueChanged.connect(self.on_transformation)
        self.ui.p.child('transformation').child('scale y-axis').sigValueChanged.connect(self.on_transformation)
        self.ui.p.child('transformation').child('x-offset').sigValueChanged.connect(self.on_transformation)
        self.ui.p.child('transformation').child('y-offset').sigValueChanged.connect(self.on_transformation)
        self.ui.p.child('transformation').child('negate').sigValueChanged.connect(self.on_transformation)
        self.ui.p.child('transformation').child('reverse').sigValueChanged.connect(self.on_transformation)
        self.ui.p.child('transformation').child('first&last point on same level').sigValueChanged.connect(self.on_transformation)
        self.ui.p.child('transformation').child('linear regression').sigValueChanged.connect(self.on_transformation)
        self.ui.p.child('transformation').child('y-zero @ x-zero').sigValueChanged.connect(self.on_transformation)
        self.ui.p.child('transformation').child('n^th derivative').sigValueChanged.connect(self.on_transformation)
        self.ui.p.child('transformation').child('n^th integral').sigValueChanged.connect(self.on_transformation)

        self.ui.p.child('filter').child('select filter').sigValueChanged.connect(self.on_filter)
        self.ui.p.child('filter').child('gauss').child('sigma').sigValueChanged.connect(self.on_filter)
        self.ui.p.child('filter').child('butterworth').child('type').sigValueChanged.connect(self.on_filter)
        self.ui.p.child('filter').child('butterworth').child('order').sigValueChanged.connect(self.on_filter)
        self.ui.p.child('filter').child('butterworth').child('low_cut').sigValueChanged.connect(self.on_filter)
        self.ui.p.child('filter').child('butterworth').child('high_cut').sigValueChanged.connect(self.on_filter)

        self.ui.p.child('filter').child('savitzky–golay').child('size').sigValueChanged.connect(self.on_filter)
        self.ui.p.child('filter').child('savitzky–golay').child('order').sigValueChanged.connect(self.on_filter)
        self.ui.p.child('filter').child('lowess').child('fraction').sigValueChanged.connect(self.on_filter)
        self.ui.p.child('filter').child('lowess').child('iteration').sigValueChanged.connect(self.on_filter)

        self.ui.p.child('moving section').child('average').sigValueChanged.connect(self.on_moving_section)
        self.ui.p.child('moving section').child('section').sigValueChanged.connect(self.on_moving_section)

        self.ui.p.child('fft').child('show fft').sigValueChanged.connect(self.on_fft)
        self.ui.p.child('fft').child('fft config').child('fft window').sigValueChanged.connect(self.on_fft)
        self.ui.p.child('fft').child('fft config').child('amplitude/phase').sigValueChanged.connect(self.on_fft)
        self.ui.p.child('fft').child('show spectrogram').sigValueChanged.connect(self.on_spectrum)
        self.ui.p.child('fft').child('spectrogram config').child('segment length').sigValueChanged.connect(self.on_spectrum)
        self.ui.p.child('fft').child('spectrogram config').child('overlap').sigValueChanged.connect(self.on_spectrum)
        self.ui.p.child('fft').child('spectrogram config').child('fft window').sigValueChanged.connect(self.on_spectrum)
        self.ui.p.child('fft').child('spectrogram config').child('log_z_scale').sigValueChanged.connect(self.on_spectrum)
        self.ui.p.child('extend axis').child('extend data').sigValueChanged.connect(self.set_extend)
        self.ui.p.child('extend axis').child('extend').child('extend to +Pos').sigValueChanged.connect(self.on_extend)
        self.ui.p.child('extend axis').child('extend').child('extend to -Pos').sigValueChanged.connect(self.on_extend)
        self.ui.p.child('extend axis').child('extend').child('fit style').sigValueChanged.connect(self.on_extend)
        self.ui.p.child('extend axis').child('extend').child('- fit length').sigValueChanged.connect(self.on_extend)
        self.ui.p.child('extend axis').child('extend').child('+ fit length').sigValueChanged.connect(self.on_extend)
        self.ui.p.child('extend axis').child('extend').child('export for PPMAC').sigActivated.connect(self.on_export)
        self.ui.p.child('extend axis').child('export to excel').sigActivated.connect(self.on_export_excel)
        self.ui.p.child('polar plot').child('use polar plot').sigValueChanged.connect(self.on_polar)
        self.ui.p.child('polar plot').child('polar config').child('n').sigValueChanged.connect(self.on_polar)
        self.ui.p.child('polar plot').child('polar config').child('scale coeff').sigValueChanged.connect(self.on_polar)
        self.ui.p.child('polar plot').child('polar config').child('angle offset').sigValueChanged.connect(self.on_polar)
        self.ui.p.child('polar plot').child('polar config').child('remove f_fund').sigValueChanged.connect(self.on_polar)
        self.ui.p.child('polar plot').child('polar config').child('show f_fund').sigValueChanged.connect(self.on_polar)
        self.ui.p.sigValueChanged.connect(self.on_treechange)

    def on_selctbtn(self):
        print('on_selectbtn')
        self.ui.p.child('transformation').hide()
        self.ui.p.child('filter').hide()
        self.ui.p.child('moving section').hide()
        self.ui.p.child('fft').hide()
        self.ui.p.child('extend axis').hide()
        self.ui.p.child('polar plot').hide()
        self.ui.p.child('autoscale graph').hide()
        self.path = self.data.selectfile(self.ui.p.child('filepath').value())
        os.path.basename(self.path)
        self.ui.p.child('config data').show()
        self.ui.p.child('filepath').setValue(self.path)
        self.ui.p.child('filename').setValue(os.path.basename(self.path))
        self.ui.p.child('filename').show()
        self.plotcolor = self.ui.p.child('transformation').child('plotcolor').value()
        self.update_ui_w2()
        self.on_colselect()

    def on_colselect(self):
        print('on_colselect')
        self.ui.w3.plotItem.hideAxis('left')
        self.ui.w3.plotItem.hideAxis('bottom')
        self.ui.p.child('config data').show()
        self.data.head = self.ui.p.child('config data').child('head rows').value()
        self.data.delimiter = self.ui.p.child('config data').child('delimiter').value()
        self.data.separator = self.ui.p.child('config data').child('decimal separator').value()
        self.data.xcol = self.ui.p.child('config data').child('x-col').value()
        self.data.ycol = self.ui.p.child('config data').child('y-col').value()
        self.data.nthval = self.ui.p.child('config data').child('use every n^th row').value()
        if self.ui.p.child('config data').child('lengths compensation').value():
            self.ui.p.child('config data').child('head rows').setValue(1)
        try:
            len_comp = self.ui.p.child('config data').child('lengths compensation').value()
            self.data.load(len_comp)
            self.on_treechange()
        except:
            print('LOAD ERROR')
            self.ui.w2.clear()

    def on_treechange(self):
        print('---------------------on_treechange')
        self.ui.treefile = Path('tree.prm')
        with open('tree.prm', 'wb') as fp:
                pickle.dump(self.ui.p.saveState(), fp)

    def on_plotcolor(self):
        print('on_plotcolor')
        self.plotcolor = self.ui.p.child('transformation').child('plotcolor').value()
        self.update_ui_w2()
        self.ui.w2.addItem(self.ui.region, ignoreBounds=True)
        self.update_ui_w3()

    def update_ui_w2(self):
        print('update_ui_w2')
        self.ui.w2.clear()
        self.ui.w3.clear()
        self.ui.d2.setStretch(x=1000, y=100)
        self.ui.d3.setStretch(x=1000, y=600)
        self.ui.d4.setStretch(x=1000, y=0)
        self.ui.w3.addItem(self.ui.plottext)
        self.ui.plottext.setText(self.data.rows)
        try:
            self.ui.w2.plot(self.data.x_raw, self.data.y_raw, pen={'color': self.plotcolor, 'width': 1}, name= 'w2')
            self.ui.wi2.setText(self.data.wi2Text)
            self.ui.w2.autoRange()
        except:
            print('Plot ERROR')

    def on_enter(self):
        print('on-enter')
        self.update_ui_w2()
        self.ui.p.child('config data').hide()
        self.ui.p.child('transformation').show()
        self.ui.p.child('filter').show()
        self.ui.p.child('moving section').show()
        self.ui.p.child('fft').show()
        self.ui.p.child('extend axis').show()
        self.ui.p.child('polar plot').show()
        self.ui.p.child('autoscale graph').show()
        self.ui.wi2.setFixedWidth(180)
        self.ui.wi3.setFixedWidth(180)
        self.ui.wi4.setFixedWidth(180)
        self.ui.w2.addItem(self.ui.region, ignoreBounds=True)
        self.ui.region.setRegion([self.data.x_raw[0], self.data.x_raw[-1]])
        self.ui.w3.scene().sigMouseMoved.connect(self.mouseMovedW3)
        self.ui.w4.scene().sigMouseMoved.connect(self.mouseMovedW4)
        self.ui.w5.view.scene().sigMouseMoved.connect(self.mouseMovedW5)
        self.on_treechange()

    def on_region_change(self):
        print('on_region_change')
        with pg.BusyCursor():
            self.region = self.ui.region.getRegion()
            self.data.region_change(self.region)
            self.on_transformation(self.ui.none)

    def on_transformation(self, p):
        with pg.BusyCursor():
            print('on_transformation')
            print('p.value() =', p.value())

            self.transform_v_arr = ['n^th derivative', 'n^th integral', 'scale x-axis', 'scale y-axis', 'x-offset',
                                    'y-offset', 'negate', 'reverse', 'y-zero @ x-zero', 'linear regression',
                                    'first&last point on same level']
            self.transform_p_arr = [self.derivative, self.integral, self.x_scale, self.y_scale, self.x_offset,
                                    self.y_offset, self.neg, self.rev, self.y_zero, self.lreg, self.flpz]

            # ch = self.ui.p.child('transformation').children()
            # for i, key in enumerate(ch):
            #     if key.name() == p.name():
            #         print('p.name() = ', p.name())
            #         print('key.name() = ', key.name())
            #         break

            if p.name() == 'n^th derivative':
                self.derivative = p.value()
            if p.name() == 'n^th integral':
                self.integral = p.value()
            if p.name() == 'scale x-axis':
                self.x_scale = p.value()
            if p.name() == 'scale y-axis':
                self.y_scale = p.value()
            if p.name() == 'x-offset':
                self.x_offset = p.value()
            if p.name() == 'y-offset':
                self.y_offset = p.value()
            if p.name() == 'negate':
                self.neg = p.value()
            if p.name() == 'reverse':
                self.rev = p.value()
            if p.name() == 'y-zero @ x-zero':
                self.y_zero = p.value()
            if p.name() == 'linear regression':
                self.lreg = p.value()
            if p.name() == 'first&last point on same level':
                self.flpz = p.value()
            self.data.transformation(self.derivative, self.integral, self.x_scale, self.y_scale,  self.x_offset,
                                     self.y_offset, self.neg, self.rev, self.y_zero, self.lreg, self.flpz)
            self.ui.p.child('transformation').child('n^th derivative').sigValueChanged.connect(self.on_transformation)
            self.ui.p.child('transformation').child('n^th integral').sigValueChanged.connect(self.on_transformation)

            self.on_treechange()
            self.on_filter(self.ui.none)

    def on_filter(self, param):
        with pg.BusyCursor():
            print('on_filter')
            self.filter = self.ui.p.child('filter').child('select filter').value()
            if self.filter == '- no filter -':
                self.ui.p.child('filter').child('gauss').hide()
                self.ui.p.child('filter').child('butterworth').hide()
                self.ui.p.child('filter').child('savitzky–golay').hide()
                self.ui.p.child('filter').child('lowess').hide()

                self.data.gauss(0, False)
                config_filter = ''

            if self.filter == 'gauss':
                self.ui.p.child('filter').child('gauss').show()
                self.ui.p.child('filter').child('butterworth').hide()
                self.ui.p.child('filter').child('savitzky–golay').hide()
                self.ui.p.child('filter').child('lowess').hide()
                self.sigma = self.ui.p.child('filter').child('gauss').child('sigma').value()
                config_filter = ' - - - ' + str(self.data.gauss(self.sigma, True))

            if self.filter == 'butterworth':
                self.ui.p.child('filter').child('gauss').hide()
                self.ui.p.child('filter').child('butterworth').show()
                self.ui.p.child('filter').child('savitzky–golay').hide()
                self.ui.p.child('filter').child('lowess').hide()

                self.ui.p.child('filter').child('butterworth').child('high_cut').hide()
                self.ui.p.child('filter').child('butterworth').child('low_cut').hide()

                #self.sigma = self.ui.p.child('filter').child('gauss').child('sigma').value()
                butter_type = self.ui.p.child('filter').child('butterworth').child('type').value()
                butter_order = self.ui.p.child('filter').child('butterworth').child('order').value()
                butter_lowcut = self.ui.p.child('filter').child('butterworth').child('low_cut').value()
                butter_highcut = self.ui.p.child('filter').child('butterworth').child('high_cut').value()

                if butter_type == 'lowpass':
                    self.ui.p.child('filter').child('butterworth').child('high_cut').show()
                    self.ui.p.child('filter').child('butterworth').child('low_cut').hide()
                    config_filter = ' - - - ' + str(self.data.butterworth(butter_order, butter_highcut, butter_type))

                if butter_type == 'highpass':
                    self.ui.p.child('filter').child('butterworth').child('low_cut').show()
                    self.ui.p.child('filter').child('butterworth').child('high_cut').hide()
                    config_filter = ' - - - ' + str(self.data.butterworth(butter_order, butter_lowcut, butter_type))

                if butter_type == 'bandpass':
                    self.ui.p.child('filter').child('butterworth').child('low_cut').show()
                    self.ui.p.child('filter').child('butterworth').child('high_cut').show()
                    config_filter = ' - - - ' + str(self.data.butterworth(butter_order, [butter_lowcut, butter_highcut], butter_type))

                if butter_type == 'bandstop':
                    self.ui.p.child('filter').child('butterworth').child('low_cut').show()
                    self.ui.p.child('filter').child('butterworth').child('high_cut').show()
                    config_filter = ' - - - ' + str(
                    self.data.butterworth(butter_order, [butter_lowcut, butter_highcut], butter_type))

            if self.filter == 'savitzky–golay':
                self.ui.p.child('filter').child('gauss').hide()
                self.ui.p.child('filter').child('butterworth').hide()
                self.ui.p.child('filter').child('savitzky–golay').show()
                self.ui.p.child('filter').child('lowess').hide()
                self.size = self.ui.p.child('filter').child('savitzky–golay').child('size').value()
                self.order = self.ui.p.child('filter').child('savitzky–golay').child('order').value()
                config_filter = ' - - - ' + str(self.data.s_g(self.size, self.order))

            if self.filter == 'lowess':
                self.ui.p.child('filter').child('gauss').hide()
                self.ui.p.child('filter').child('butterworth').hide()
                self.ui.p.child('filter').child('savitzky–golay').hide()
                self.ui.p.child('filter').child('lowess').show()
                self.fraction = self.ui.p.child('filter').child('lowess').child('fraction').value()
                self.iteration = self.ui.p.child('filter').child('lowess').child('iteration').value()
                config_filter = ' - ' + str(self.data.lowess(self.fraction, self.iteration))
            self.title =  str(os.path.basename(self.path)) + ' - ' + config_filter
            self.ui.w2.setTitle(self.title)
            self.ui.w3.setTitle(self.title)
            self.ui.w4.setTitle(self.title)
            self.on_moving_section()

    def on_moving_section(self):
        p = self.ui.p.child('moving section').child('average').value()
        if p == '- none -':
            print('- none -')
            self.DWR = ''
        section = self.ui.p.child('moving section').child('section').value() / self.x_scale
        with pg.BusyCursor():
            if p == 'simple_moving_average':
                self.data.simple_moving_average(section)
            if p == 'moving p2v':
                self.data.moving_p2v(section)
            if p == 'moving p2v no slope':
                self.data.moving_p2v_no_slope(section)
        self.on_fft()


    def on_fft(self):
        print('on_fft')
        with pg.BusyCursor():
            show_fft = self.ui.p.child('fft').child('show fft').value()
            phase = self.ui.p.child('fft').child('fft config').child('amplitude/phase').value()
            self.ui.w4.setTitle(self.title)
            if self.ui.p.child('fft').child('show fft').value():
                self.ui.p.child('fft').child('show spectrogram').setValue(False)
                self.win = self.ui.p.child('fft').child('fft config').child('fft window').value()
                j = 0
                for i in self.ui.fftwindows[0]:
                    if i.__contains__(self.win):
                        win_coeff = self.ui.fftwindows[1][j]
                    else:
                        j = j + 1
                self.data.fft(self.win, win_coeff)
            self.update_ui_w4(True, show_fft, phase)
            self.on_spectrum()

    def on_spectrum(self):
        print('on_spectrum')
        with pg.BusyCursor():
            if self.ui.p.child('fft').child('show spectrogram').value():
                print('in')
                self.length = self.ui.p.child('fft').child('spectrogram config').child('segment length').value()
                self.overlap = self.ui.p.child('fft').child('spectrogram config').child('overlap').value()
                self.win = self.ui.p.child('fft').child('spectrogram config').child('fft window').value()
                self.log_z_scale = self.ui.p.child('fft').child('spectrogram config').child('log_z_scale').value()
                self.f, self.t, self.Sxx = self.data.spectrogram(self.length, self.overlap, self.win, self.log_z_scale)
                self.spxscale = self.t[-1] / self.Sxx.shape[0]
                self.spyscale = self.f[-1] / self.Sxx.shape[1]
            self.update_ui_w5()
            self.on_extend()

    def set_extend(self):
        print('set_extend')
        with pg.BusyCursor():
            self.ui.p.child('extend axis').child('extend').child('extend to -Pos').setValue(-10.0 + np.floor(self.data.current_x[0]))
            self.ui.p.child('extend axis').child('extend').child('extend to +Pos').setValue(10.0+np.ceil(self.data.current_x[-1]))
            self.on_extend()

    def on_extend(self):
        print('on_extend')
        self.ui.p.child('extend axis').child('extend').hide()
        self.extend_info = ''
        with pg.BusyCursor():
            if self.ui.p.child('extend axis').child('extend data').value():
                self.ui.p.child('extend axis').child('extend').show()
                plus = self.ui.p.child('extend axis').child('extend').child('extend to +Pos').value()
                minus = self.ui.p.child('extend axis').child('extend').child('extend to -Pos').value()
                fit_plus = self.ui.p.child('extend axis').child('extend').child('+ fit length').value()
                fit_minus = self.ui.p.child('extend axis').child('extend').child('- fit length').value()
                fit = self.ui.p.child('extend axis').child('extend').child('fit style').value()
                self.data.extend(minus, plus, fit_minus, fit_plus, fit)
            self.on_polar()

    def on_export(self):
        print('on_export')
        with pg.BusyCursor():
            stop = self.ui.p.child('extend axis').child('extend').child('extend to +Pos').value()
            offset = self.ui.p.child('extend axis').child('extend').child('extend to -Pos').value()
            length = stop - offset
            N = str(len(self.data.extend_array))
            self.data.export(length, N , offset)

    def on_export_excel(self):
        print('on_export_excel')
        with pg.BusyCursor():
            section = self.ui.p.child('moving section').child('section').value()
            if self.ui.p.child('moving section').child('average').value() != '- none -':
                av = True
            else:
                av = False
            self.data.excel_export(self.ui.p, av, section)

    def on_polar(self):
        print('on_polar')
        with pg.BusyCursor():
            if self.ui.p.child('polar plot').child('use polar plot').value():
                self.ui.p.child('polar plot').child('polar config').show()
                n = self.ui.p.child('polar plot').child('polar config').child('n').value()
                coeff = self.ui.p.child('polar plot').child('polar config').child('scale coeff').value()
                angle = self.ui.p.child('polar plot').child('polar config').child('angle offset').value()
                remove_f = self.ui.p.child('polar plot').child('polar config').child('remove f_fund').value()
                show_f = self.ui.p.child('polar plot').child('polar config').child('show f_fund').value()
                self.data.polar(n, coeff, angle, remove_f, show_f)
                fund_i =''
                if remove_f:
                    fund_i = 'f_fund removed'
                self.ui.w3.setTitle(self.title + ' - n='+str(n)+' - '+ fund_i)
            else:
                self.ui.p.child('polar plot').child('polar config').hide()
                self.ui.w3.setTitle(self.title)
            self.update_ui_w3()
            self.ui.wi3.setText(self.data.wi3Text)

    def update_ui_w3(self):
        print('update_ui_w3')
        with pg.BusyCursor():
            self.ui.w3.clear()
            self.ui.w3.addItem(self.ui.w3vLine, ignoreBounds=True)
            self.ui.w3.addItem(self.ui.w3hLine, ignoreBounds=True)
            self.ui.w3.plotItem.showAxis('left')
            self.ui.w3.plotItem.showAxis('bottom')
            if self.ui.p.child('polar plot').child('use polar plot').value():
                self.ui.w3.plotItem.hideAxis('left')
                self.ui.w3.plotItem.hideAxis('bottom')
                self.ui.w3.setAspectLocked(lock=True, ratio=1)
                self.ui.w3.addLine(x=0, pen={'color': (200, 200, 200), 'width': 1})
                self.ui.w3.addLine(y=0, pen={'color': (200, 200, 200), 'width': 1})
                self.ui.w3.plot(self.data.xx, self.data.yy, pen={'color': self.plotcolor, 'width': 1})
                self.ui.w3.addItem(self.data.circleMax)
                self.ui.w3.addItem(self.data.circleMin)
                if self.ui.p.child('polar plot').child('polar config').child('remove f_fund').value():
                    self.ui.w3.plot(self.data.radius_f_x, self.data.radius_f_y, pen={'color': (0, 0, 250), 'width': 3})
                    self.ui.w3.setTitle(self.title +'polar plot - fft-window = ' + self.win+ ' - f_fund removed')
                else:
                    self.ui.w3.plot(self.data.radius_x, self.data.radius_y, pen={'color': self.plotcolor, 'width': 1})
                    self.ui.w3.setTitle(self.title + 'polar plot - fft-window = ' + self.win)
                if self.ui.p.child('polar plot').child('polar config').child('show f_fund').value():
                    self.ui.w3.plot(self.data.fund_x, self.data.fund_y, pen={'color': (0, 250, 250), 'width': 1})
            else:
                self.ui.w3.plotItem.showAxis('left')
                self.ui.w3.plotItem.showAxis('bottom')
                self.ui.w3.setAspectLocked(lock=False, ratio=None)
                self.ui.w3.plot(self.data.current_x, self.data.current_y, pen={'color': self.plotcolor, 'width': 1}, name='region')
                if self.ui.p.child('moving section').child('average').value() != '- none -':
                    self.ui.w3.plot(self.data.x_transformed, self.data.y_devdist, pen={'color': (230, 50, 50), 'width': 3}, name='dev sec')
                self.ui.wi3.setText(self.data.wi3Text)
                if self.ui.p.child('extend axis').child('extend data').value():
                    self.ui.w3.plot(self.data.X_extend_1, self.data.Y_extend_1, pen={'color': (100, 100, 255), 'width': 3}, name='extended_1')
                    self.ui.w3.plot(self.data.X_extend_2, self.data.Y_extend_2, pen={'color': (100, 100, 255), 'width': 3}, name='extended_2')
                    self.ui.w3.plot(self.data.X_fit_1, self.data.Y_fit_1, pen={'color': (100, 255, 100), 'width': 3}, name='fit_1')
                    self.ui.w3.plot(self.data.X_fit_2, self.data.Y_fit_2, pen={'color': (100, 255, 100), 'width': 3},name='fit_2')
            if self.ui.p.child('autoscale graph').value():
                self.ui.w3.autoRange()
            self.on_treechange()
            self.data.info_w3()

    def mouseMovedW3(self, pos):
        print('mouseMovedW3')
        with pg.BusyCursor():
            vb = self.ui.w3.getPlotItem().vb
            self.mousePoint = vb.mapSceneToView(pos)
            self.ui.w3vLine.show()
            self.ui.w3hLine.show()
            self.ui.w3vLine.setPos(self.mousePoint.x())
            self.ui.w3hLine.setPos(self.mousePoint.y())
            if self.ui.p.child('polar plot').child('use polar plot').value():
                z = self.mousePoint.x() + self.mousePoint.y() * 1j
                radius = np.absolute(z)
                angle = np.angle(z, deg=True)
                self.data.wi3Pos = ' radius = ' + str(round(radius - self.data.offset, 6)) + '\n' + \
                                   ' angle = ' + str(round(angle, 6)) + '\n' + '\n'
            else:
                self.data.wi3Pos = ' X = ' + str(round(self.mousePoint.x(), 6)) + '\n' + \
                                   ' Y = ' + str(round(self.mousePoint.y(), 6)) + '\n' + '\n'
            self.data.info_w3()
            self.ui.wi3.setText(self.data.wi3Text)

    def update_ui_w4(self, mouse, show_fft, phase):
        with pg.BusyCursor():
            if mouse:
                self.ui.w4.clear()
                self.ui.w4.addItem(self.ui.w4vLine, ignoreBounds=True)
                self.ui.w4.addItem(self.ui.w4hLine, ignoreBounds=True)
                if show_fft:
                    self.ui.p.child('fft').child('fft config').show()
                    self.ui.d2.setStretch(x=1000, y=100)
                    self.ui.d3.setStretch(x=1000, y=100)
                    self.ui.d4.setStretch(x=1000, y=500)
                    self.ui.d5.setStretch(x=1000, y=0)
                    if phase:
                      self.ui.w4.plot(self.data.x_f, self.data.y_p, pen={'color': self.plotcolor, 'width': 1})
                      self.ui.w4.setTitle( self.title + ' - polar spectrum - fft-window=' + self.win )
                    else:
                      self.ui.w4.plot(self.data.x_f, self.data.y_f, pen={'color': self.plotcolor, 'width': 1})
                      self.ui.w4.setTitle( self.title + ' - amplitude spectrum - fft-window=' + self.win )

                else:
                    self.ui.p.child('fft').child('fft config').hide()
                    self.ui.d2.setStretch(x=1000, y=100)
                    self.ui.d3.setStretch(x=1000, y=600)
                    self.ui.d4.setStretch(x=1000, y=0)
                    self.ui.d5.setStretch(x=1000, y=0)
            self.data.info_w4()
            self.ui.wi4.setText(self.data.wi4Text)
            if self.ui.p.child('autoscale graph').value():
                self.ui.w4.autoRange()

    def mouseMovedW4(self, pos):
        print('mouseMovedW4')
        with pg.BusyCursor():
            vb = self.ui.w4.getPlotItem().vb
            self.mousePoint = vb.mapSceneToView(pos)
            self.ui.w4vLine.show()
            self.ui.w4hLine.show()
            self.ui.w4vLine.setPos(self.mousePoint.x())
            self.ui.w4hLine.setPos(self.mousePoint.y())
            self.data.w4_x_point = round(self.mousePoint.x(), 6)
            self.data.w4_y_point = round(self.mousePoint.y(), 6)
            self.update_ui_w4(False, False, False)

    def update_ui_w5(self):
        global Sxx
        print('update_ui_w5')
        with pg.BusyCursor():
            self.ui.w5.clear()
            if self.ui.p.child('fft').child('show spectrogram').value():
                self.ui.p.child('fft').child('show fft').setValue(False)
                self.ui.p.child('fft').child('fft config').hide()
                self.ui.p.child('fft').child('spectrogram config').show()
                self.ui.w5.addItem(self.ui.w5vLine, ignoreBounds=True)
                self.ui.w5.addItem(self.ui.w5hLine, ignoreBounds=True)
                self.ui.p.child('fft').child('spectrogram config').show()
                self.ui.p.child('fft').child('show fft').setValue(False)
                self.ui.d2.setStretch(x=1000, y=0)
                self.ui.d3.setStretch(x=1000, y=100)
                self.ui.d4.setStretch(x=1000, y=0)
                self.ui.d5.setStretch(x=1000, y=600)
                self.ui.p.child('fft').child('spectrogram config').show()
                self.ui.w5.setImage(self.Sxx, scale=[self.spxscale, self.spyscale], pos=[0, 0])
                self.ui.w5.view.invertY(False)
                self.ui.w5.setPredefinedGradient('greyclip')
                self.ui.w5.view.setTitle(str(self.title) + ' - spectrogram - fft-win=' + str(self.win) \
                                         + ' - length=' + str(self.length) + ' - overlap=' + str(self.overlap))
                self.ui.w5.view.setAspectLocked(lock=False)
            else:
                self.ui.p.child('fft').child('spectrogram config').hide()
                self.ui.d2.setStretch(x=1000, y=100)
                self.ui.d3.setStretch(x=1000, y=600)
                self.ui.d5.setStretch(x=1000, y=0)
            if self.ui.p.child('autoscale graph').value():
                self.ui.w5.view.autoRange()

    def mouseMovedW5(self, pos):
        print('mouseMovedW5')
        with pg.BusyCursor():
            vb = self.ui.w5.view.vb
            self.mousePoint = vb.mapSceneToView(pos)
            self.ui.w5vLine.show()
            self.ui.w5hLine.show()
            self.ui.w5vLine.setPos(self.mousePoint.x())
            self.ui.w5hLine.setPos(self.mousePoint.y())
            self.data.w5_x_point = round(self.mousePoint.x(), 3)
            self.data.w5_y_point = round(self.mousePoint.y(), 3)
            if (self.data.w5_x_point >= 0) and (self.data.w5_x_point <= self.t[-1]) and \
                        (self.data.w5_y_point >= 0) and (self.data.w5_y_point <= self.f[-1]):
                self.data.w5_z_point = round(self.Sxx.item(int(self.data.w5_x_point/self.spxscale), int(self.data.w5_y_point/self.spyscale)), 9)
                self.ui.w5.view.setTitle(str(self.title) + ' - spectrogram - fft-win=' + str(self.win) \
                                    + ' - length=' + str(self.length) + ' - overlap=' + str(self.overlap)+ '\n'+\
                                     ' X='+str()+' y='+str(self.data.w5_y_point)+' z='+str(self.data.w5_z_point))
            else:
                self.ui.w5.view.setTitle(str(self.title) + ' - spectrogram - fft-win=' + str(self.win) \
                                         + ' - length=' + str(self.length) + ' - overlap=' + str(self.overlap) + '\n'+\
                                         ' X=' + str(self.data.w5_x_point) + ' y=' + str(
                    self.data.w5_y_point) + ' z= nan')

