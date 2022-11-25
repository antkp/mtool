import os
import pyqtgraph as pg
import numpy as np
import scipy as sp
import itertools

from xlsxwriter import Workbook
from scipy.signal import get_window
from scipy import signal
from scipy.optimize import curve_fit
from pyqtgraph.Qt import QtCore
from scipy.ndimage import gaussian_filter1d
from scipy import integrate
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

# todo
#  autodetect delimiter & dezimal separator
#  "switch" paramtree preset doesen´t work
#  bei Anwendung von Filtern im "Hintergrund" die ungefilterten Daten anzeigen (auch bei Excel export)

class Data(QtCore.QObject):
    sig_data_loaded = QtCore.pyqtSignal()

    def __init__(self):
        super().__init__()
        self.filepath = ''
        self.av_section_info = ''
        self.av_name = ''
        self.RMS = ''
        self.head = 0
        self.delimiter = ','
        self.separator = '.'
        self.xcol = 2
        self.ycol = 1
        self.nthval = 1
        self.xscale = 1.0
        self.yscale = 1.0
        self.rows = ''
        self.n = 0
        self.T_raw = 0.0
        self.T =0.0
        self.p2v = 0.0
        self.w3_x_point = 0.0
        self.w3_y_point = 0.0
        self.w4_x_point = 0.0
        self.w4_y_point = 0.0
        self.fund_ampl = 0.0
        self.fund_phase = 0.0

        self.x_raw = []
        self.y_raw = []
        self.raw_region_x = []
        self.raw_region_y = []
        self.current_x = []
        self.current_y = []
        self.x_scaled = []
        self.y_scaled = []
        self.x_transformed = []
        self.y_transformed = []
        self.comparray = []
        self.x_filtered = []
        self.y_filtered = []
        self.xx = []
        self.yy = []
        self.x_f = []
        self.y_f = []
        self.y_p = []
        self.plus_fit = []
        self.minus_fit = []
        self.plus_extend = []
        self.minus_extend = []
        self.extend_array = []
        self.y_fund = []
        self.fft_band_arr =[]


        self.circleMax = pg.QtGui.QGraphicsEllipseItem()
        self.circleMid = pg.QtGui.QGraphicsEllipseItem()
        self.circleMin = pg.QtGui.QGraphicsEllipseItem()

        self.delta = 0
        self.reg_min = 0
        self.reg_max =1
        self.dx = 0.0
        self.minXi = 1
        self.maxXi = 2

        self.wi2Text = ''
        self.wi3Text = ''
        self.wi3Pos = ''
        self.wi4Text = ''
        self.extend_info = ''
        self.fund_info = ''
        self.val_dict = {
            "N": 0,
            "t": 0.0,
            "min": 0.0,
            "max": 0.0,
            "p2v": 0.0,
            "RMS": 0.0,
            'sec_info': 0.0}



    def selectfile(self, path):

        # todo evtl. "preload mit Vorschau nur wenige sec versuchen danach stoppen"
        print('selectfile')
        qfd = pg.FileDialog()
        filter = "Text files (*.txt *lsdf *.csv)"
        self.selectedFile = pg.FileDialog.getOpenFileName(qfd, "select file", path, filter)
        print('selected 1 = ', self.filepath)
        self.filepath = self.selectedFile[0]
        self.rows = ''
        with open(self.filepath, newline='', mode='r') as f:
            for i in range(20):
                self.rows = self.rows + f.readline()
        return self.filepath


    def load(self, len_comp):
        print('-load-')

        #num_lines = sum(1 for _ in open(self.filepath)) # number lines in file
        #print('num_lines = ', num_lines)


        #search file reversed for lines onliy containing charater listed in l
        skip_footer = 0
        l = ['0', ';']
        for line in reversed(list(open(self.filepath, 'r'))):
            result = all(elem in l for elem in list(set(line.rstrip())))
            if not result:
                break
            skip_footer += 1

        text = open(self.filepath, 'r').read()

        #print('---text = ', text)

        if self.separator == ',':
            text = text.replace(',', '.')
            print(len(text))

        #print('list(set(data))', list(set(text)))  # todo use also to autodetect separators !

        #todo .find wenn möglich zeilenweise --> dann '0;0;0' zu begin der Zeile suchen?
        #
        #p = text.find('0;0;0;0;0;')
        #text = text[0:p]

        text = text.replace('nan', '0')
        temp_file = open("temp", "w+")
        temp_file.write(text)

        if self.xcol != -1:
            print('xcol = ', self.xcol)
            if self.nthval != 1:
                print('nthval = ', self.nthval)
                with open('temp') as f_in:
                    data = np.genfromtxt(itertools.islice(f_in, 0, None, self.nthval), skip_footer= skip_footer,
                        skip_header=self.head,delimiter=self.delimiter, usecols=[self.xcol, self.ycol], autostrip=True)
            else:
                print('nthval ++ = ', self.nthval)
                data = np.genfromtxt('temp', skip_header=self.head, delimiter=self.delimiter, skip_footer= skip_footer,
                                usecols=[self.xcol, self.ycol], autostrip=True)


            print('type(data) = ', type(data))
            print('data = ',  data)

            self.x_raw = data[:, 0]
            self.y_raw = data[:, 1]
        else:
            print('xcol ## = ', self.xcol)
            if self.nthval != 1:
                print('nthval = ', self.nthval)
                with open('temp') as f_in:
                    self.y_raw = np.genfromtxt(itertools.islice(f_in, 0, None, self.nthval), skip_header=self.head, skip_footer= skip_footer,
                                delimiter=self.delimiter, usecols=[self.ycol], autostrip=True)
            else:
                print('nthval ++ = ', self.nthval)
                print('###___###')
                self.y_raw = np.genfromtxt('temp', skip_header=self.head, delimiter=self.delimiter, skip_footer= skip_footer,
                                usecols=[self.ycol], autostrip=True)

            self.x_raw = np.arange(start=0, stop=len(self.y_raw))

        temp_file.close()
        os.remove("temp")


        if len_comp:
            g = self.x_raw[0] - self.x_raw[-1]
            h = self.y_raw[0] - self.y_raw[-1]
            signg = np.sign(g)
            signh = np.sign(h)

            y_1 = self.y_raw - self.y_raw[0]
            x_1 = self.x_raw - self.x_raw[0]
            if signg == signh:
                self.y_raw = y_1 - x_1
            else:
                self.y_raw = y_1 + x_1

        self.current_x = self.x_raw * self.xscale
        self.current_y = self.y_raw * self.yscale
        self.T = abs(round(self.x_raw[11] - self.x_raw[10], 6))
        #print('self.T = ', self.T)
        self.wi2Text = '\n' + ' T = ' + str(self.T) + '\n' + \
                  ' fs = ' + str(round((1 / self.T), 8)) + '\n' + \
                  ' N = ' + str(len(self.x_raw)) + '\n' + \
                  ' t = ' + str(round(len(self.x_raw) * self.T, 8))
        #print('self.wi2Text = ', self.wi2Text)
        self.sig_data_loaded.emit()

    def region_change(self, region):
        print('region_change')
        a, b = region
        self.reg_min = (np.abs(self.x_raw - a)).argmin()
        self.reg_max = (np.abs(self.x_raw - b)).argmin()
        self.raw_region_x = self.x_raw[self.reg_min:self.reg_max]
        self.raw_region_y = self.y_raw[self.reg_min:self.reg_max]
        self.cut_x = self.raw_region_x - self.raw_region_x[0]
        self.cut_y = self.raw_region_y
        self.current_x = self.cut_x
        self.current_y = self.cut_y

    def linreg(self, x, y):
        res_a = stats.linregress(x, y)
        i = (res_a.intercept)
        s = (res_a.slope)
        reg = -(i + s * x) + y
        return reg

    def transformation(self, derivative, integral, x_scale, y_scale, x_offset, y_offset, neg, rev,  y_zero, lreg, flpz):
        print('transformation')
        self.x_transformed = (self.cut_x * x_scale) + x_offset
        self.y_transformed = (self.cut_y * y_scale) + y_offset

        dx = self.x_transformed[1] - self.x_transformed[0]
        print('dx = ', dx)
        if neg:
            self.y_transformed = -self.y_transformed
        if rev:
            self.y_transformed = self.y_transformed[::-1]
        if lreg:
            self.y_transformed = self.linreg(self.x_transformed, self.y_transformed)
        if flpz:
            self.delta = self.y_transformed[-1] - self.y_transformed[0]
            self.n = len(self.y_transformed)
            self.comparray = np.arange(0, self.delta, (self.delta / self.n))
            if self.n != len(self.comparray):
                self.y_transformed = self.y_transformed - self.comparray[:self.n]
            else:
                self.y_transformed = self.y_transformed - self.comparray
        if derivative != 0:
            self.x_transformed = self.x_transformed[:-derivative]+(derivative*dx)
            self.y_transformed = np.diff(self.raw_region_y, n=derivative)
        if integral == 1:
            self.y_transformed = integrate.cumtrapz(self.y_transformed, self.x_transformed, initial=0.0)
        if integral == 2:
            self.y_transformed = integrate.cumtrapz(self.y_transformed, self.x_transformed, initial=0.0)
            self.y_transformed = integrate.cumtrapz(self.y_transformed, self.x_transformed, initial=0.0)
        if y_zero:
            idx = np.searchsorted(self.x_transformed, 0, side="left")
            self.y_transformed = self.y_transformed - self.y_transformed[idx]
        self.current_y = self.y_transformed
        self.current_x = self.x_transformed
        self.y_filtered = self.y_transformed

    def gauss(self, sigma, act):
        f_control = 'active filter: gauss - sigma =' + str(sigma)
        print(f_control)
        if act:
            self.y_filtered = gaussian_filter1d(self.y_transformed, sigma)
            self.current_y = self.y_filtered
        else:
            self.y_filtered = self.y_transformed
            self.current_y = self.y_transformed
        return f_control

    def butterworth(self, order, Wn, type):
        ny = (1 / self.T)/2
        b, a = signal.butter(order, Wn /ny, type, analog=False)
        self.y_filtered = signal.filtfilt(b, a, self.y_transformed)
        self.current_y = self.y_filtered

    def s_g(self, size, order):
        f_control = 'active filter: s_g - size=' + str(size) + ' order=' + str(order)
        print(f_control)
        self.y_filtered = signal.savgol_filter(self.y_transformed, size, order)
        self.current_y = self.y_filtered
        return f_control

    def lowess(self, fraction, iteration):
        f_control = 'active filter: lowess - fraction=' + str(fraction) + ' iteration=' + str(iteration)
        print('in lowess ' + str(fraction) + ' ' + str(iteration))
        self.y_filtered = lowess(self.y_transformed, np.arange(len(self.y_transformed)), is_sorted=True, frac=fraction,
                                 it=iteration)
        self.y_filtered = self.y_filtered[:, 1]
        self.current_y = self.y_filtered
        return f_control

    def fft_filter(self):

        # todo update fft graph if fiter was modified

        x_f = np.fft.fftfreq(len(self.x_transformed), self.T * self.xscale)
        x_f = np.fft.fftshift(x_f)

        df = (1 / (self.T** self.xscale)) / len(self.x_transformed)

        fft_filter_koeff_array_x = np.arange(0, (1 / self.T) / (2 * self.xscale) + df, df)
        print('fft_filter_koeff_array_x = ', fft_filter_koeff_array_x)

        a = np.ones(round(len(self.x_transformed)/2))
        print('a_1 = ', a)
        print('len(a_1) = ', len(a))

        for i in self.fft_band_arr:
            print('i = ', i )
            print('(i[1]-i[0])/df = ', (i[1]-i[0])/df)
            band =round((i[1]-i[0])/df)

            print('df_fftfreq= ', df)
            print('x_f[2] -x_f[1]', x_f[2] -x_f[1])

            print('i[1] = ', i[1])
            print('i[0] = ', i[0])
            print('band = ', band)
            W = signal.windows.hann(band)

            b = 1-W
            b = b * 0
            start = round(i[0] / df)

            stop = start + band
            print('start = ', start)
            print('stop = ', stop)

            print('len a_2 = ', len(a))
            print('len b = ', len(b))

            c = np.append(a[:start], b)
            c = np.append(c, a[stop:])
            print('len(c) = ', len(c))
            a = c

        a_rev = a[::-1]
        a = np.append(a_rev, a)
        print('len(a_3) = ', len(a))

        sp = np.fft.fft(self.y_transformed)
        sp = np.fft.fftshift(sp)
        print('len(a) = ', len(a))
        print('len(sp) = ', len(sp))
        spa = sp * a
        spa = np.fft.ifftshift(spa)
        print('x_f = ', x_f)
        print('sp.real = ', spa.real)

        iff = np.fft.ifft(spa)
        self.current_y =  iff.real
        self.y_filtered = self.current_y

        return x_f, spa.real,  iff.real

    def simple_moving_average(self, section):
        print('simple_moving_average')
        self.av_name = 'simple_moving_average'
        self.y_devdist = np.arange(len(self.y_filtered)) * np.nan
        for i in range(len(self.y_filtered) - int(section)):
            y_arr = self.y_filtered[i:i + int(section)]
            self.y_devdist[i + int(section / 2)] = np.sum(y_arr)/len(y_arr)
        x = self.y_devdist
        x = x[~np.isnan(x)]
        self.av_section_info = self.av_name + '\n max= ' + str(round(max(x), 6)) + '\n min= ' +\
                               str(round(min(x), 6)) + '\n p2v= ' + str(round(max(x)-min(x), 6))

    def moving_p2v(self,  section):
        print('moving_p2v')
        self.av_name = 'moving_p2v'
        self.y_devdist = np.arange(len(self.y_filtered)) * np.nan
        for i in range(len(self.y_filtered) - int(section)):
            y_arr = self.y_filtered[i:i + int(section)]
            self.y_devdist[i + int(section / 2)] = np.ptp(y_arr)
        x = self.y_devdist
        x = x[~np.isnan(x)]
        self.av_section_info = self.av_name + '\n max= ' + str(round(max(x), 6))

    def moving_p2v_no_slope(self, section):
        print('moving_p2v_no_slope')
        self.av_name = 'moving_p2v_no_slope'
        self.y_devdist = np.arange(len(self.y_filtered)) * np.nan
        for i in range(len(self.y_filtered) - int(section)):
            y_arr = self.y_filtered[i:i + int(section)]
            x_arr = self.x_transformed[i:i + int(section)]
            y_reg = self.linreg(x_arr, y_arr)
            self.y_devdist[i + int(section / 2)] = np.ptp(y_reg)
        x = self.y_devdist
        x = x[~np.isnan(x)]
        self.av_section_info = self.av_name + '\n max= ' + str(round(max(x), 6))

    # todo fft log scale Y ?

    def fft(self, win, win_coeff):
        N = len(self.y_filtered)  # number of samples
        self.x_f = np.linspace(0.0, 1.0 / (2.0 * self.T * self.xscale), N // 2)
        if win == 'boxcar':
            amp = np.fft.fft(self.y_filtered)
            # angle_f = np.rad2deg(np.angle(self.y_filtered[0:N // 2]))
        else:
            W = sp.signal.windows.get_window(win, len(self.y_filtered))
            a = self.y_filtered * W * win_coeff
            amp = np.fft.fft(a)
        self.y_f = 2.0 / N * np.abs(amp[0:N // 2])
        self.y_p = np.rad2deg(np.angle(amp[0:N // 2]))
        ind_a = np.where(self.y_f == max(self.y_f))[0]
        self.fund_ampl = max(self.y_f)
        self.fund_phase = self.y_p[ind_a[0]]

    def spectrogram(self, length, overlap, win, log_z_scale):
        print('spectrogram')
        f, t, Sxx = signal.spectrogram(self.y_filtered, fs=(1 / self.T), nperseg=length, noverlap=overlap, window=win, scaling='spectrum')
        Sxx = np.sqrt(Sxx)  # from V**2 to V
        Sxx = np.transpose(Sxx)
        if log_z_scale:
            Sxx = 20 * np.log10(np.abs(Sxx))
        return f, t, Sxx

    def extend(self, minus, plus, fit_minus, fit_plus, fit):
        # todo extend offset evtl manuell einstellen ?
        print('extend')
        l_minus = int(fit_minus / self.T)
        l_plus = int(fit_plus / self.T)
        self.X_fit_1 = self.current_x[: l_minus]
        self.X_fit_2 = self.current_x[-l_plus:]
        self.Y_fit_1 = self.y_filtered[: l_minus]
        self.Y_fit_2 = self.y_filtered[-l_plus:]
        print('minus = ', minus, '   self.current_x[0] = ', self.current_x[0], ' self.T = ', self.T)
        self.X_extend_1 = np.arange(minus, self.current_x[0], self.T)
        self.X_extend_2 = np.arange(self.current_x[-1] + self.T, plus, self.T)

        if fit == 'flat':
            self.Y_extend_1 = np.full(len(self.X_extend_1), self.y_filtered[0])
            self.Y_extend_2 = np.full(len(self.X_extend_2), self.y_filtered[-1])
        if fit == 'linear':
            popt_1, pcov_1 = curve_fit(self.lin_func, self.X_fit_1, self.Y_fit_1)
            popt_2, pcov_2 = curve_fit(self.lin_func, self.X_fit_2, self.Y_fit_2)
            self.Y_extend_1 = self.lin_func(self.X_extend_1, *popt_1)
            self.Y_extend_2 = self.lin_func(self.X_extend_2, *popt_2)
        if fit == 'quadratic':
            popt_1, pcov_1 = curve_fit(self.qadr_func, self.X_fit_1, self.Y_fit_1)
            popt_2, pcov_2 = curve_fit(self.qadr_func, self.X_fit_2, self.Y_fit_2)
            self.Y_extend_1 = self.qadr_func(self.X_extend_1, *popt_1)
            self.Y_extend_2 = self.qadr_func(self.X_extend_2, *popt_2)
        if fit == 'cubic':
            popt_1, pcov_1 = curve_fit(self.cube_func, self.X_fit_1, self.Y_fit_1)
            popt_2, pcov_2 = curve_fit(self.cube_func, self.X_fit_2, self.Y_fit_2)
            self.Y_extend_1 = self.cube_func(self.X_extend_1, *popt_1)
            self.Y_extend_2 = self.cube_func(self.X_extend_2, *popt_2)
        if fit == 'sin':
            popt_1, pcov_1 = curve_fit(self.sin_func, self.X_fit_1, self.Y_fit_1)
            popt_2, pcov_2 = curve_fit(self.sin_func, self.X_fit_2, self.Y_fit_2)
            self.Y_extend_1 = self.sin_func(self.X_extend_1, *popt_1)
            self.Y_extend_2 = self.sin_func(self.X_extend_2, *popt_2)
        if fit == 'exp':
            popt_1, pcov_1 = curve_fit(self.exp_func, self.X_fit_1, self.Y_fit_1)
            popt_2, pcov_2 = curve_fit(self.exp_func, self.X_fit_2, self.Y_fit_2)
            self.Y_extend_1 = self.exp_func(self.X_extend_1, *popt_1)
            self.Y_extend_2 = self.exp_func(self.X_extend_2, *popt_2)

        self.extend_array = np.append(self.Y_extend_1, self.y_filtered)
        self.extend_array = np.append(self.extend_array, self.Y_extend_2)
        self.extend_info = '\n' + '\n' + ' T = ' + str(self.T) + '\n' + ' N_ext = ' + str(len(self.extend_array))

    def lin_func(self, x, a, b, c, d):
        return a*x+b

    def qadr_func(self, x, a, b, c, d):
        return a*x*x+b*x+c

    def cube_func(self, x, a, b, c, d):
        return a*x*x*x+b*x*x+c*x+d

    def sin_func(self, x, a, b, c, d):
        return a*np.sin(b*x+c)+d

    def exp_func(self, x, a, b, c, d):
        return a * np.exp(-b * x) + c

    def export(self, length, N, offset): # todo export to PPMAC
        print('export')
        dirname = os.path.dirname(self.filepath)
        exportfilename = os.path.basename(self.filepath)
        exportfilename = str(os.path.splitext(exportfilename)[0])
        exportfilename = str(exportfilename)+'_length='+str(length)+'_points='+str(N)+'_offset='+str(offset)+'.csv'
        expfile = os.path.join(dirname, exportfilename)

        export_filename = pg.FileDialog.getSaveFileName(None, "export file", expfile, "CSV files (*.csv)")
        file = open(export_filename[0], 'w')
        np.savetxt(file, self.extend_array.reshape(-1, 1), delimiter=',', fmt='%10.9f')
        file.close()

    def excel_export(self, p, av, section, fft):
        dirname = os.path.dirname(self.filepath)
        exportfile = os.path.basename(self.filepath)
        exportfile = str(os.path.splitext(exportfile)[0])
        exportfilename = str(exportfile)+'.xlsx'
        expfile = os.path.join(dirname, exportfilename)
        expfile = dirname+'/'+ exportfilename
        export_filename = pg.FileDialog.getSaveFileName(None, "export file as", expfile, "EXCEL files (*.xlsx)")

        workbook = Workbook(str(export_filename[0]))

        worksheet1 = workbook.add_worksheet(exportfile)  # Required for the chart data.
        worksheet1.write_column('A2', self.current_x)
        worksheet1.write('A1', 'base')
        worksheet1.write_column('B2', self.current_y)
        worksheet1.write('B1', exportfilename)

        i = 0
        for key, value in self.val_dict.items():
            worksheet1.write(i, 14, key)
            worksheet1.write(i, 15, value)
            i = i+1

        chart_1 = workbook.add_chart({'type': 'scatter', 'subtype': 'straight'})
        chart_1.add_series({
            'name': exportfilename,
            'line': {'width': 1, 'color': '#002060'},
            'categories': '='+exportfile+'!$A$2:$A$'+str(len(self.current_x)),
            'values': '='+exportfile+'!$B$2:$B$'+str(len(self.current_y))})

        chart_1.set_size({'width': 604.5, 'height': 312.41})
        chart_1.set_plotarea({'layout': {'x': 1, 'y': 0.1, 'width': 0.85, 'height': 0.75, }})
        chart_1.set_title({'name': exportfilename, 'font': {'size': 9}, 'layout': {'x': 0.1, 'y': 0.0}})
        chart_1.set_x_axis({'name': 'x-axis', 'label_position': 'low', 'name_layout': {'x': 0.45, 'y': 0.95, },
                      'major_gridlines': {'visible': True, 'line': {'width': 0.2, 'color':'#808080', }}, })
        chart_1.set_y_axis({'name': 'y-axis', 'label_position': 'low', 'name_layout': {'x': 0.0, 'y': 0.45, },
                      'major_gridlines': {'visible': True, 'line': {'width': 0.2, 'color':'#808080', }}, })
        #chart_1.set_legend({'position': 'top'})
        chart_1.set_legend({'none': True})
        worksheet1.insert_chart('D2', chart_1)

        if av:
            for row, value in enumerate(self.y_devdist):
                try:
                    worksheet1.write(row+1, 2,  value)
                except:
                    pass
            worksheet1.write('C1', self.av_name + '_' +exportfilename + '_period='+str(section))
            chart_2 = workbook.add_chart({'type': 'scatter', 'subtype': 'straight'})
            chart_2.add_series({
                'name': self.av_name + exportfilename,
                'line': {'width': 1, 'color': '#C00000'},
                'categories': '=' + exportfile + '!$A$2:$A$' + str(len(self.current_x)),
                'values': '=' + exportfile + '!$C$2:$C$' + str(len(self.current_y))})

            chart_2.set_size({'width': 604.5, 'height': 312.41})
            chart_2.set_plotarea({'layout': {'x': 1, 'y': 0.1, 'width': 0.85, 'height': 0.75, }})
            chart_2.set_title({'name': exportfilename + ' moving PV', 'font': {'size': 9}, 'layout': {'x': 0.1, 'y': 0.0}})
            chart_2.set_x_axis({'name': 'x-axis', 'label_position': 'low', 'name_layout': {'x': 0.45, 'y': 0.95, },
                              'major_gridlines': {'visible': True, 'line': {'width': 0.2, 'color': '#808080', }}, })
            chart_2.set_y_axis({'name': 'y-axis', 'label_position': 'low', 'name_layout': {'x': 0.0, 'y': 0.45, },
                              'major_gridlines': {'visible': True, 'line': {'width': 0.2, 'color': '#808080', }}, })
            #chart_2.set_legend({'position': 'top'})
            chart_2.set_legend({'none': True})
            worksheet1.insert_chart('D18', chart_2)

        if fft:

            worksheet2 = workbook.add_worksheet('FFT')
            worksheet2.write_column('A2', self.x_f)
            worksheet2.write('A1', 'frequency')
            worksheet2.write_column('B2', self.y_f)
            worksheet2.write('B1', 'FFT-'+exportfilename+' amplitude')
            worksheet2.write_column('C2', self.y_p)
            worksheet2.write('C1', 'FFT-'+exportfilename+' phase')


            chart_3 = workbook.add_chart({'type': 'scatter', 'subtype': 'straight'})
            chart_3.add_series({
                'name': 'FFT-'+exportfilename+' amplitude',
                'line': {'width': 1, 'color': '#002060'},
                'categories': 'FFT'+'!$A$2:$A$'+str(len(self.x_f)),
                'values': 'FFT'+'!$B$2:$B$'+str(len(self.y_f))})

            chart_3.set_size({'width': 604.5, 'height': 312.41})
            chart_3.set_plotarea({'layout': {'x': 1, 'y': 0.1, 'width': 0.85, 'height': 0.75, }})
            chart_3.set_title({'name': 'FFT-'+exportfilename+' amplitude', 'font': {'size': 9}, 'layout': {'x': 0.1, 'y': 0.0}})
            chart_3.set_x_axis({'name': 'frequency', 'label_position': 'low', 'name_layout': {'x': 0.45, 'y': 0.95, },
                          'major_gridlines': {'visible': True, 'line': {'width': 0.2, 'color':'#808080', }}, })
            chart_3.set_y_axis({'name': 'amplitude', 'label_position': 'low', 'name_layout': {'x': 0.0, 'y': 0.45, },
                          'major_gridlines': {'visible': True, 'line': {'width': 0.2, 'color':'#808080', }}, })
            #chart_3.set_legend({'position': 'top'})
            chart_3.set_legend({'none': True})
            worksheet2.insert_chart('D2', chart_3)

            chart_4 = workbook.add_chart({'type': 'scatter', 'subtype': 'straight'})
            chart_4.add_series({
                'name': 'FFT-'+exportfilename+' phase',
                'line': {'width': 1, 'color': '#002060'},
                'categories': 'FFT' + '!$A$2:$A$'+str(len(self.x_f)),
                'values': 'FFT' + '!$C$2:$C$'+str(len(self.y_p))})

            chart_4.set_size({'width': 604.5, 'height': 312.41})
            chart_4.set_plotarea({'layout': {'x': 1, 'y': 0.1, 'width': 0.85, 'height': 0.75, }})
            chart_4.set_title({'name': 'FFT-'+exportfilename+' phase', 'font': {'size': 9}, 'layout': {'x': 0.1, 'y': 0.0}})
            chart_4.set_x_axis({'name': 'frequency', 'label_position': 'low', 'name_layout': {'x': 0.45, 'y': 0.95, },
                          'major_gridlines': {'visible': True, 'line': {'width': 0.2, 'color':'#808080', }}, })
            chart_4.set_y_axis({'name': 'phase', 'label_position': 'low', 'name_layout': {'x': 0.0, 'y': 0.45, },
                          'major_gridlines': {'visible': True, 'line': {'width': 0.2, 'color':'#808080', }}, })
            #chart_4.set_legend({'position': 'top'})
            chart_4.set_legend({'none': True})
            worksheet2.insert_chart('D18', chart_4)

        worksheet3 = workbook.add_worksheet('config')
        self.i_row = 0
        self. i_column = 0
        self.treech(p, worksheet3)
        workbook.close()

    def treech(self, p , sheet):
        if p.hasChildren():
            self.i_row = self.i_row + 1
            print('row =', self.i_row, ' column', self.i_column, ' BRANCH name = ', p.name())
            sheet.write(self.i_row, self.i_column, p.name())
            self.i_column = self.i_column + 1
            for key in p.children():
                self.treech(p.child(key.name()), sheet)
            self.i_column = self.i_column - 1
        else:
            self.i_row = self.i_row + 1
            print('row =', self.i_row, 'column', self.i_column, ' name = ', p.name(), ' value =', p.value())
            sheet.write(self.i_row, self.i_column, p.name())
            sheet.write(self.i_row, self.i_column + 1, str(p.value()))

    def polar(self, n, coeff, angle, remove_f, show_f):
        self.fund_info = ''
        self.p2v = max(self.y_filtered) - min(self.y_filtered)
        self.n = n
        self.offset = min(self.y_filtered) + (self.p2v / 2) + coeff * self.p2v
        self.theta = np.linspace(np.deg2rad(angle), (2 * np.pi * self.n) + np.deg2rad(angle), len(self.y_filtered))
        fft_amp_arr = np.fft.fft(self.y_filtered)
        f_a = np.abs(fft_amp_arr)
        fund = np.where(f_a == max(f_a[1:]))[0]
        fund_arr = fft_amp_arr * 0
        fund_arr[fund] = fft_amp_arr[fund]
        fund_arr[-fund] = fft_amp_arr[-fund]
        self.fund_sin = np.fft.ifft(fund_arr)
        self.fund_sin = self.fund_sin.real
        self.radius_f = self.y_filtered - self.fund_sin + self.offset
        self.radius = self.y_filtered + self.offset
        self.radius_x = self.radius * np.cos(self.theta)
        self.radius_y = self.radius * np.sin(self.theta)
        self.radius_f_x = self.radius_f * np.cos(self.theta)
        self.radius_f_y = self.radius_f * np.sin(self.theta)
        self.radius_fund = self.fund_sin.real + self.offset
        self.fund_x = self.radius_fund * np.cos(self.theta)
        self.fund_y = self.radius_fund * np.sin(self.theta)

        if remove_f:
            self.current_y = self.y_filtered - self.fund_sin
            self.r = self.radius_f
        else:
            self.current_y = self.y_filtered
            self.r = self.radius

        if show_f:
            print('show_f on')
            self.fund_info = 'p2v_fund =' + str(round(max(self.fund_sin) - min(self.fund_sin), 6))
        self.circleMax = pg.QtGui.QGraphicsEllipseItem(-max(self.r), -max(self.r), max(self.r) * 2, max(self.r) * 2)
        self.circleMin = pg.QtGui.QGraphicsEllipseItem(-min(self.r), -min(self.r), min(self.r) * 2, min(self.r) * 2)
        self.circleMax.setPen(pg.mkPen(color=(200, 200, 0), width=2, style=QtCore.Qt.DotLine))
        self.circleMin.setPen(pg.mkPen(color=(200, 200, 0), width=2, style=QtCore.Qt.DotLine))

    def rms(self, array):
        return np.sqrt(np.mean(np.square(array)))

    def info_w3(self):
        print('info_w3')

        self.val_dict = {
            "N": len(self.current_y),
            "t": round(len(self.current_y) * self.T, 6),
            "min": round(min(self.current_y), 6),
            "max": round(max(self.current_y), 6),
            "p2v": round(self.p2v, 6),
            "RMS": round(self.rms(self.current_y), 6),
            'sec_info': self.av_section_info}

        self.p2v = round(max(self.current_y) - min(self.current_y), 8)
        self.T = round(self.x_transformed[2] - self.x_transformed[1], 8)
        self.wi3Text = \
                 '\n' + \
                 self.wi3Pos +\
                 ' N = '   + str(len(self.current_y)) + '\n' + \
                 ' t = '   + str(round(len(self.current_y) * self.T, 6)) + '\n' + '\n' + \
                 ' min = ' + str(round(min(self.current_y), 6)) + '\n' + \
                 ' max = ' + str(round(max(self.current_y), 6)) + '\n' + \
                 ' p2v = ' + str(round(self.p2v, 6)) + '\n' + \
                 ' RMS = ' + str(round(self.rms(self.current_y), 6)) + '\n' + '\n' + '\n' +\
                 self.av_section_info + '\n' +\
                 self.extend_info + '\n' + \
                 self.fund_info

    def info_w4(self):
        print('info_w4')
        self.wi4Text = '\n' + ' X = ' + str(self.w4_x_point) + '\n' + ' Y = ' + str(self.w4_y_point) + '\n'