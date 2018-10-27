#!/usr/bin/env python
#-*-coding: utf-8-*-
"""

@author: xuyihe

History:
    Yihe Xu Mar 2 2018: Initial version is copyed from
            check_and_modify_disp_qt4.py
"""

#DEBUG = True
DEBUG = False
import warnings
warnings.filterwarnings('ignore')

import os
import sys
from ConfigParser import ConfigParser

from GVImage import GrpVelImg
from GVImage_pick import (track_peak, track_peak_forward,
                               track_peak_backward)
from Disp import Disp

#import EGFpy as ftpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from PyQt4 import QtCore, QtGui

if DEBUG:
    import pdb

from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas)
from matplotlib.figure import Figure

class DispCheckModifyWidget(FigureCanvas):

    selectChanged = QtCore.pyqtSignal(object, bool)
    stateChanged = QtCore.pyqtSignal(object, object)

    def __init__(self, parent=None, width=5, height=4,
                 dpi=100):

        self.fig = Figure(figsize=(width, height), dpi=dpi,
                         facecolor='none', frameon=True,
                         linewidth=5)
        self.ax = self.fig.add_axes([0.1, 0.1, 0.8, 0.8], facecolor='none')

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding,
                                   #QtGui.QSizePolicy.Minimum,
                                   #QtGui.QSizePolicy.Minimum,
                                  )
        FigureCanvas.updateGeometry(self)

    def set_EGF_file(self, EGF_file):
        self.EGFfile = EGF_file

    def set_Disp_file(self, Disp_file):
        self.dispfile = Disp_file

    def set_sta_pair(self, sta_pair):
        self.sta_pair = sta_pair

    def get_EGF_file(self):
        return self.EGFfile

    def get_Disp_file(self):
        return self.dispfile

    def get_sta_pair(self):
        return self.sta_pair

    def set_gv_para(self, gvmin, gvmax, gvdelta):
        self.gvmin = gvmin
        self.gvmax = gvmax
        self.gvdelta = gvdelta

    def set_T_para(self, Tmin, Tmax, Tdelta):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Tdelta = Tdelta

    def set_cfile_format(self, cfile_format):
        self.cfile_format = cfile_format

    def set_accept_state(self, state):
        assert(state in ['P', 'A', 'R'])
        self.accept = state

    def get_gv_para(self):
        return self.gvmin, self.gvmax, self.gvdelta

    def get_T_para(self):
        return self.Tmin, self.Tmax, self.Tdelta

    def get_cfile_format(self):
        return self.cfile_format

    def get_accept_state(self):
        return self.accept

    def save_disp(self, folder, filename=None):
        if filename is None:
            filename = self.disp.filename
        self.disp.save(os.path.join(folder, filename))

    def plot_init(self):
        self.disp = Disp(self.dispfile)

        self.gvi = GrpVelImg(self.EGFfile,
                      cfile_format=self.cfile_format,
                      gvmin   = self.gvmin,
                      gvmax   = self.gvmax,
                      gvdelta = self.gvdelta,
                      Tmin    = self.Tmin,
                      Tmax    = self.Tmax,
                      Tdelta  = self.Tdelta)

        mpl.rcParams['xtick.labelsize'] = 8
        mpl.rcParams['ytick.labelsize'] = 8

        self.gvimage = self.ax.imshow(self.gvi.img, origin='lower',
                   extent=self.gvi.extent,
                   cmap='jet', aspect=self.gvi.aspect)

        self.mask_origin = self.disp.masked_disp.mask.copy()
        self.displine_origin, = self.ax.plot(self.disp.T, self.disp.masked_disp,
                                        '.', color='gray', lw=2)
        self.displine, = self.ax.plot(self.disp.T, self.disp.masked_disp, '.w', lw=2)

        for vg, c in zip([1.0, 2.0, 3.0], ['y--', 'g--', 'b--']):
            x = self.gvi.dist / (2.0*vg)
            self.ax.plot([x,x],[0,1], c, transform=self.ax.get_xaxis_transform())

        if DEBUG:
            print(self.disp.masked_disp)

        #self.accept = 'P'
        self.accept_dict = dict(P='Pending', A='Accept', R='Reject')
        self.accept_color_dict = dict(P='none', A='green', R='red')
        self.accept_text = self.ax.text(1,1, self.accept_dict[self.accept],
                                        ha='right', va='top',
                                    transform=self.ax.transAxes)
        self.fig.set_facecolor(self.accept_color_dict[self.accept])

        self.sta_pair_text = self.ax.text(0.5, 1, self.sta_pair, ha='center',
                                          va = 'top',
                                          transform=self.ax.transAxes)

        self.ax.set_xlim(*self.gvi.extent[:2])
        self.ax.set_ylim(*self.gvi.extent[2:])

        self.ax.set_xlabel('Period (s)', fontsize=10)
        self.ax.set_ylabel('Velocity (km/s)', fontsize=10)

        self.bline = None
        self.eline = None
        self.b = 2*self.gvi.T[0] - self.gvi.T[1]
        self.e = 2*self.gvi.T[-1] - self.gvi.T[-2]

        self.modify_state = False
        self.shift_state = False
        self.alt_state = False

        self.selected = False

        self.track_disp_line = None

    def connect(self):
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('figure_enter_event', self.mouse_enter)
        self.fig.canvas.mpl_connect('figure_leave_event', self.mouse_leave)

    def onclick(self, event):
        if DEBUG:
            if event.xdata is None:
                print("Out of axes")
            else:
                if event.button == 1:
                    print("      Click on (%.2f,%.2f)" % (event.xdata, event.ydata))
                elif event.button == 3:
                    print("Right Click on (%.2f,%.2f)" % (event.xdata, event.ydata))
        if event.button == 1:
            if self.modify_state:
                x, y = event.xdata, event.ydata
                xi = np.argmin(np.abs(x - self.disp.T))
                self.disp.masked_disp.data[xi] = y
            elif self.shift_state:
                T0, gv0 = event.xdata, event.ydata
                gvmin, gvmax, gvdelta = self.get_gv_para()
                Tmin, Tmax, Tdelta = self.get_T_para()
                i0 = int((T0-Tmin)/Tdelta)
                j0 = int((gv0-gvmin)/gvdelta)
                self.tracked_disp = track_peak(i0, j0, self.gvi.img)
                self.tracked_disp = self.tracked_disp * gvdelta + gvmin
                self.disp.masked_disp.data = self.tracked_disp.copy()
                #if self.track_disp_line is None:
                #    self.track_disp_line, = \
                #        self.ax.plot(self.gvi.T, self.tracked_disp, 'w+')
                #else:
                #    self.track_disp_line.set_ydata(self.tracked_disp)
            elif self.alt_state:
                pass

            else:
                self.b = event.xdata
                if self.b is None: # out of axes
                    self.b = 2*self.gvi.T[0]-self.gvi.T[1]
                if self.bline is None:
                    self.bline, = self.ax.plot([self.b, self.b], [0, 1],
                                         'w', lw=2,
                                         transform=self.ax.get_xaxis_transform())
                else:
                    self.bline.set_data([self.b,self.b], [0,1])

        elif event.button == 3:
            if self.modify_state:
                if not self.shift_state:
                    T0, gv0 = event.xdata, event.ydata
                    gvmin, gvmax, gvdelta = self.get_gv_para()
                    Tmin, Tmax, Tdelta = self.get_T_para()
                    i0 = int((T0-Tmin)/Tdelta)
                    j0 = int((gv0-gvmin)/gvdelta)
                    tracked_disp_f = track_peak_forward(i0, j0, self.gvi.img)
                    tracked_disp_f = tracked_disp_f * gvdelta + gvmin
                    self.disp.masked_disp.data[i0:] = tracked_disp_f[i0:]
                else:
                    T0, gv0 = event.xdata, event.ydata
                    gvmin, gvmax, gvdelta = self.get_gv_para()
                    Tmin, Tmax, Tdelta = self.get_T_para()
                    i0 = int((T0-Tmin)/Tdelta)
                    j0 = int((gv0-gvmin)/gvdelta)
                    tracked_disp_b = track_peak_backward(i0, j0, self.gvi.img)
                    tracked_disp_b = tracked_disp_b * gvdelta + gvmin
                    self.disp.masked_disp.data[:i0+1] = tracked_disp_b[:i0+1]

            elif self.shift_state:
                T0, gv0 = event.xdata, event.ydata
                gvmin, gvmax, gvdelta = self.get_gv_para()
                Tmin, Tmax, Tdelta = self.get_T_para()
                i0 = int((T0-Tmin)/Tdelta)
                j0 = int((gv0-gvmin)/gvdelta)
                self.tracked_disp = track_peak(i0, j0, self.gvi.img)
                self.tracked_disp = self.tracked_disp * gvdelta + gvmin
                self.disp.masked_disp.data[:] = self.tracked_disp[:]

            elif self.alt_state:
                pass
            else:
                self.e = event.xdata
                if self.e is None: # out of axes
                    self.e = 2*self.gvi.T[-1]-self.gvi.T[-2]
                if self.eline is None:
                    self.eline, = self.ax.plot([self.e, self.e], [0, 1],
                                         'y', lw=2,
                                         transform=self.ax.get_xaxis_transform())
                else:
                    self.eline.set_data([self.e,self.e], [0,1])

        self.newmask = (self.disp.T < self.b) | (self.disp.T > self.e) \
                | self.mask_origin

        if DEBUG:
            print self.disp.T
            print self.b, self.e
            print self.mask_origin
            print self.newmask

        self.disp.masked_disp.mask = self.newmask
        self.displine.set_data(self.disp.T, self.disp.masked_disp)

        self.accept = 'A'
        self.stateChanged.emit(self.sta_pair, 'A')
        self.fig.set_facecolor('green')
        self.accept_text.set_text('Accept')

        self.fig.canvas.draw()

    def onmotion(self, event):
        if DEBUG:
            print('Mouse position (matplotlib): {}, {}'.format(\
                  event.xdata, event.ydata))

    def mouse_enter(self, event):
        self.selected = True
        self.selectChanged.emit(self.sta_pair, self.selected)
        print('mouse_enter (emit): {} {}'.format(self.sta_pair, self.selected))
        self.fig.set_edgecolor('red')
        self.fig.canvas.draw()
        if DEBUG:
            print('Mouse Enter:')

    def mouse_leave(self, event):
        self.selected = False
        self.selectChanged.emit(self.sta_pair, self.selected)
        print('mouse_enter (emit): {} {}'.format(self.sta_pair, self.selected))
        self.fig.set_edgecolor('white')
        self.fig.canvas.draw()
        if DEBUG:
            print('Mouse Leave:')

    def onkeypress(self, event):
        if DEBUG:
            print('Key pressed: %s' % event.key())
        if event.key() == QtCore.Qt.Key_R:
            if DEBUG:
                print('Key pressed: R')
            self.accept = 'R'
            self.stateChanged.emit(self.sta_pair, 'R')
            self.fig.set_facecolor('red')
            self.accept_text.set_text('Reject')
            self.fig.canvas.draw()
        elif event.key() == QtCore.Qt.Key_A:
            if DEBUG:
                print('Key pressed: A')
            self.accept = 'A'
            self.stateChanged.emit(self.sta_pair, 'A')
            self.fig.set_facecolor('green')
            self.accept_text.set_text('Accept')
            self.fig.canvas.draw()
        elif event.key() == QtCore.Qt.Key_P:
            self.accept = 'P'
            self.stateChanged.emit(self.sta_pair, 'P')
            self.fig.set_facecolor('none')
            self.accept_text.set_text('Pending')
            self.fig.canvas.draw()
        elif event.key() == QtCore.Qt.Key_Control:
            if DEBUG:
                print('Key pressed: Control')
            self.modify_state = True
        elif event.key() == QtCore.Qt.Key_Shift:
            self.shift_state = True
            if DEBUG:
                print('Key pressed: Shift')
                print('Shift_state: {}'.format(self.shift_state))
        elif event.key() == QtCore.Qt.Key_Alt:
            self.alt_state = True

        elif event.key() == QtCore.Qt.Key_F1:
            self.pop_help_window()
        #elif event.key() == QtCore.Qt.Key_F:
        #    self.forward_tracking = True
        #elif event.key() == QtCore.Qt.Key_B:
        #    self.backward_tracking = True

    def onkeyrelease(self, event):
        if DEBUG:
            print('Key released: %s' % event.key())
        if event.key() == QtCore.Qt.Key_Control:
            if DEBUG:
                print('Key released: Control')
            self.modify_state = False
        elif event.key() == QtCore.Qt.Key_Shift:
            self.shift_state = False
            if DEBUG:
                print('Key released: Shift')
                print('Shift_state: {}'.format(self.shift_state))
        elif event.key() == QtCore.Qt.Key_Alt:
            self.alt_state = False
        #    self.forward_tracking = False
        #elif event.key() == QtCore.Qt.Key_B:
        #    self.backward_tracking = False


    def pop_help_window(self):
        Help_message = """Keys shortcut:
Left click: L Limit T   Right click: L Limit T
Shift:
    click      : pick (automatic)
Control:
    click      : pick (manually)
    right click: pick (auto, leftwards)
    Shift +
    right click: pick (auto, rightwards)
        """
        reply = QtGui.QMessageBox.information(self,
            "Help",
            Help_message,
            QtGui.QMessageBox.Ok)


class MainWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle('Check and Modify Program (Version 0.1)')

        self.main_widget = QtGui.QWidget(self)

        l = QtGui.QVBoxLayout(self.main_widget)

        self.main_cnm_figure = DispCheckModifyWidget(self.main_widget,
                width=5, height=4, dpi=100)
        self.main_cnm_figure.set_gv_para(0.5, 4, 0.002)
        self.main_cnm_figure.set_T_para(0.5, 5, 0.1)
        self.main_cnm_figure.set_EGF_file('../demo_data/CF.BA.040-BA.116.sac')
        self.main_cnm_figure.set_cfile_format('SAC_CF')
        self.main_cnm_figure.set_Disp_file('../demo_data/GDisp.BA.040-BA.116.dat')
        self.main_cnm_figure.set_sta_pair('040-116')
        self.main_cnm_figure.set_accept_state('P')
        self.main_cnm_figure.plot_init()
        self.main_cnm_figure.connect()

        l.addWidget(self.main_cnm_figure)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        self.statusBar().showMessage('Ready!', 2000)

    def closeEvent(self, event):
        self.close()

    def mouseMoveEvent(self, event):
        if DEBUG:
            print('Mouse pos: {}'.format(event.pos()))

    def keyPressEvent(self, event):
        if DEBUG:
            print('Key pressed: {}'.format(event.key()))
            print('event.key(): {}'.format(type(event.key())))
            #print('GVImage selected: {}'.format(\
                #self.main_cnm_figure.selected))

        if self.main_cnm_figure.selected:
            self.main_cnm_figure.onkeypress(event)

    def keyReleaseEvent(self, event):
        if DEBUG:
            print('Key pressed: {}'.format(event.key()))
        if self.main_cnm_figure.selected:
            self.main_cnm_figure.onkeyrelease(event)

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    win = MainWindow()
    if DEBUG:
        print isinstance(win, QtGui.QWidget)
        print dir(win)
    win.show()
    sys.exit(app.exec_())

