import os
import pyqtgraph as pg
import numpy as np
import scipy as sp
from scipy.signal import get_window
from scipy import signal
from scipy.optimize import curve_fit
from pyqtgraph.Qt import QtCore
from scipy.ndimage import gaussian_filter1d
from scipy import integrate
from statsmodels.nonparametric.smoothers_lowess import lowess


class Data(QtCore.QObject):
    sig_data_loaded = QtCore.pyqtSignal()

    def __init__(self):
        super().__init__()
        self.filepath = ''
        self.DWR = ''
        self.head = 0
        self.delimiter = ','
        self.separator = '.'
        self.xcol = 2
        self.ycol = 1
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
        self.y_devdist = []
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

    def selectfile(self, path): #filedialog gehört eigentlich in die GUI ?
        print('selectfile')
        qfd = pg.FileDialog()
        filter = "Text files (*.txt *lsdf *.csv)"
        self.selectedFile = pg.FileDialog.getOpenFileName(qfd, "select file", path, filter)
        self.filepath = self.selectedFile[0]
        self.rows = ''
        with open(self.filepath, newline='', mode='r') as f:
            for i in range(20):
                self.rows = self.rows + f.readline()
        return self.filepath

    def load(self, len_comp):
        print('load')
        text = open(self.filepath, 'r').read()
        p = text.find('0;0;0;0;0;')
        rtext = text[::-1]
        #print(rtext)
        xx = rtext.rfind(';0')
        print('rtext[0:xx] =',rtext[0:xx])

        print('len.text = ', len(text))
        print('p = ', p)
        print('xx= ', xx)
        print('len(text) - xx= ', len(text)-xx)


        text = text[0:p]
        text = text.replace('nan', '0')
        if self.separator == ',':
            text = text.replace(',', '.')
        temp_file = open("temp", "w+")
        temp_file.write(text)
        data = np.genfromtxt('temp', skip_header=self.head, delimiter=self.delimiter, usecols=[self.xcol, self.ycol], autostrip=True)
        temp_file.close()
        os.remove("temp")
        self.x_raw = data[:, 0]
        self.y_raw = data[:, 1]
        print('self.x_raw[0:5] = ', self.x_raw[0:5])
        print('self.y_raw[0:5] = ', self.y_raw[0:5])
        if len_comp:
            g = self.x_raw[0] - self.x_raw[-1]
            h = self.y_raw[0] - self.y_raw[-1]
            signg = np.sign(g)
            signh = np.sign(h)

            print('g = ', g)
            print('h = ', h)
            print('signg = ', signg)
            print('signh = ', signh)
            y_1 = self.y_raw - self.y_raw[0]
            x_1 = self.x_raw - self.x_raw[0]
            if signg == signh:
                self.y_raw = y_1 - x_1
            else:
                self.y_raw = y_1 + x_1
            print('self.x_raw[0:5] = ', self.x_raw[0:5])
            print('self.y_raw[0:5] = ', self.y_raw[0:5])
        self.x_scaled = self.x_raw * self.xscale
        self.y_scaled = self.y_raw * self.yscale
        self.y_lincomp = self.y_raw * self.yscale
        self.current_x = self.x_raw * self.xscale
        self.current_y = self.y_raw * self.yscale
        self.T = abs(round(self.x_scaled[11] - self.x_scaled[10], 6))
        print('self.T = ', self.T)
        self.wi2Text = '\n' + ' T = ' + str(self.T) + '\n' + \
                  ' fs = ' + str(round((1 / self.T), 8)) + '\n' + \
                  ' N = ' + str(len(self.x_raw)) + '\n' + \
                  ' t = ' + str(round(len(self.x_raw) * self.T, 8))
        print('self.wi2Text = ', self.wi2Text)
        self.sig_data_loaded.emit()

    def region_change(self, region):
        print('region_change')
        a, b = region
        self.reg_min = (np.abs(self.x_scaled - a)).argmin()
        self.reg_max = (np.abs(self.x_scaled - b)).argmin()
        self.raw_region_x = self.x_scaled[self.reg_min:self.reg_max]
        self.raw_region_y = self.y_scaled[self.reg_min:self.reg_max]
        self.current_x = self.raw_region_x
        self.current_y = self.raw_region_y



    def transformation(self, derivative, integral, x_offset, y_offset, neg, rev,  y_zero, lcomp):
        print('transformation')
        self.x_transformed = self.raw_region_x + x_offset
        self.y_transformed = self.raw_region_y + y_offset
        dx = self.x_transformed[1] - self.x_transformed[0]
        print('dx = ', dx)
        if neg:
            self.y_transformed = -self.y_transformed
        if rev:
            self.y_transformed = self.y_transformed[::-1]
        if lcomp:
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
            print(len(self.raw_region_x))
            print(len(self.x_transformed))

        if integral == 1:
            self.y_transformed = integrate.cumtrapz(self.y_transformed, self.x_transformed, initial=0.0)
            print(len(self.raw_region_x))
            print(len(self.x_transformed))

        if integral == 2:
            self.y_transformed = integrate.cumtrapz(self.y_transformed, self.x_transformed, initial=0.0)
            self.y_transformed = integrate.cumtrapz(self.y_transformed, self.x_transformed, initial=0.0)
            print(len(self.raw_region_x))
            print(len(self.x_transformed))

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

    def s_g(self, size, order):
        f_control = 'active filter: s_g - size=' + str(size) + ' order=' + str(order)
        print(f_control)
        self.y_filtered = signal.savgol_filter(self.y_transformed, size, order)
        self.current_y = self.y_filtered
        return f_control

    def lowess(self, fraction, iteration):
        f_control = 'active filter: lowess - fraction=' + str(fraction) + ' iteration=' + str(iteration)
        print('in lowess ' + str(fraction) + ' ' + str(iteration))
        self.y_filtered = lowess(self.y_transformed, np.arange(len(self.y_transformed)), is_sorted=True, frac=fraction, it=iteration)
        self.y_filtered = self.y_filtered[:, 1]
        self.current_y = self.y_filtered
        return f_control

    def dev_dist(self, section, ch):
        print('dev_dist')
        if ch:
            self.y_devdist = np.arange(len(self.y_filtered)) * np.nan
            for i in range(len(self.y_filtered) - int(section)):
                arr = self.y_filtered[i:i + int(section)]
                self.y_devdist[i + int(section / 2)] = max(arr) - min(arr)
            x = self.y_devdist
            x = x[~np.isnan(x)]
            self.DWR = '\n' + 'DWR_max = ' + str(max(x))
        else:
            self.DWR = ''

    def fft(self, win, win_coeff):
        N = len(self.y_filtered)  # number of samples
        self.x_f = np.linspace(0.0, 1.0 / (2.0 * self.T * self.xscale), N // 2)
        if win == 'boxcar':
            amp = np.fft.fft(self.y_filtered)
            angle_f = np.rad2deg(np.angle(self.y_filtered[0:N // 2]))
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
        print('win = ', win)
        f, t, Sxx = signal.spectrogram(self.y_filtered, fs=(1 / self.T), nperseg=length, noverlap=overlap, window=win, scaling='spectrum')
        Sxx = np.sqrt(Sxx)  # from V**2 to V
        Sxx = np.transpose(Sxx)
        if log_z_scale:
            Sxx = 20 * np.log10(np.abs(Sxx))
        return f, t, Sxx

    def extend(self, minus, plus, fit_minus, fit_plus, fit):
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

    def export(self):
        print('export')
        exportfile = os.path.basename(self.filepath)
        exportfile = os.path.splitext(exportfile)[0]
        export_filename = pg.FileDialog.getSaveFileName(None, "export file", os.path.basename(exportfile), "CSV files (*.csv)")
        file = open(export_filename[0], 'w')
        np.savetxt(file, self.extend_array.reshape(-1, 1), delimiter=',', fmt='%10.9f')
        file.close()

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

    def info_w3(self):
        print('info_w3')
        self.p2v = round(max(self.current_y) - min(self.current_y), 8)
        self.T = round(self.x_transformed[2] - self.x_transformed[1], 8)
        self.wi3Text = '\n' + ' min = ' + str(round(min(self.current_y), 6)) + '\n' + \
                 ' max = ' + str(round(max(self.current_y), 6)) + '\n' + \
                 ' p2v = ' + str(round(self.p2v, 6)) + '\n' + \
                 ' N = ' + str(len(self.current_y)) + '\n' + \
                 ' t = ' + str(round(len(self.current_y) * self.T, 6)) + '\n' + '\n' + \
                 self.wi3Pos + self.DWR + '\n' + '\n' + \
                 self.extend_info + '\n' + '\n' + \
                 self.fund_info

    def info_w4(self):
        print('info_w4')
        self.wi4Text = '\n' + ' X = ' + str(self.w4_x_point) + '\n' + ' Y = ' + str(self.w4_y_point) + '\n'
