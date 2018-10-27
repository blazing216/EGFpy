#!/usr/bin/env python
#-*-coding:utf-8-*-

#DEBUG = True
DEBUG = False
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

from PyQt4 import QtCore, QtGui

from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas)
from matplotlib.figure import Figure

from InterStaDist import InterStaDist

class StationMapWidget(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4,
                 dpi=100):
        self.fig = Figure(figsize=(width,height), dpi=dpi,
                          facecolor='none')
        self.ax = self.fig.add_axes([0.1, 0.1, 0.8, 0.8],
                                    facecolor='none')

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


    def set_InterStaDist(self, inter_sta_dist):
        self.inter_sta_dist = inter_sta_dist

    def plot_init(self):
        self.ax.plot(*self.inter_sta_dist.get_lonlat(),
                     marker='v', mfc='none',
                     linestyle = 'none')
        self.ax.set_aspect(1.0)
        self.sta_group1_h = None
        self.sta_group2_h = None
        self.sta1_sta2_h = None

    def connect(self):
        pass

    def set_sta_radius(self, sta1, sta2, radius1, radius2):
        self.sta1 = sta1
        self.sta2 = sta2
        self.radius1 = radius1
        self.radius2 = radius2

    def update_sta_selected(self):
        assert(self.sta1 is not None and \
               self.sta2 is not None and \
               self.radius1 is not None and \
               self.radius2 is not None)
        lon1, lat1 = self.inter_sta_dist.get_lonlat(self.sta1)
        lon2, lat2 = self.inter_sta_dist.get_lonlat(self.sta2)

        sta_group1 = self.inter_sta_dist.get_adjacent_sta(self.sta1,
                self.radius1)
        sta_group2 = self.inter_sta_dist.get_adjacent_sta(self.sta2,
                self.radius2)
        lons1, lats1 = self.inter_sta_dist.get_lonlat(sta_group1)
        lons2, lats2 = self.inter_sta_dist.get_lonlat(sta_group2)

        if self.sta_group1_h is None:
            self.sta_group1_h, = self.ax.plot(lons1, lats1, 'r+')
        else:
            self.sta_group1_h.set_data(lons1, lats1)

        if self.sta_group2_h is None:
            self.sta_group2_h, = self.ax.plot(lons2, lats2, 'g+')
        else:
            self.sta_group2_h.set_data(lons2, lats2)

        if self.sta1_sta2_h is None:
            self.sta1_sta2_h, = self.ax.plot([lon1, lon2],
                                            [lat1, lat2],
                                            'k')
        else:
            self.sta1_sta2_h.set_data([lon1, lon2],
                                      [lat1, lat2])


        self.fig.canvas.draw()

        return sta_group1, sta_group2

    def plot_sta_radius(self, *args):
        self.set_sta_radius(*args)
        #lon1, lat1 = self.inter_sta_dist.get_lonlat(sta1)
        #lon2, lat2 = self.inter_sta_dist.get_lonlat(sta2)
        sta_group1, sta_group2 = self.update_sta_selected()
        return sta_group1, sta_group2

    #def changeStaPath(self, sta_pair, selected):


#class StationMapWidget(QtGui.QWidget):
#    def __init__(self, parent=None):
#        QtGui.QWidget.__init__(self, parent)
#        #self.setWindowTitle('Station Map Program (Version 0.1)')
#        #l = QtGui.QVBoxLayout(self)
#
#        self.main_widget = QtGui.QWidget(self)
#        l = QtGui.QVBoxLayout(self.main_widget)
#
#        self.main_cnm_figure = StationMapCanvas(self.main_widget,
#                width=5, height=4, dpi=100)
#        a = InterStaDist.from_file('staloc_BC', 4, 3, 1, 0)
#        self.main_cnm_figure.set_InterStaDist(a)
#        self.main_cnm_figure.plot_init()
#        self.main_cnm_figure.connect()
#
#        self.initUI()
#
#        l.addWidget(self.main_cnm_figure)
#
#        #self.main_widget.setFocus()
#        #self.setCentralWidget(self.main_widget)
#
#        #self.statusBar().showMessage('Ready!', 2000)
#
#    def initUI(self):
#        self.main_cnm_figure.plot_sta_radius('025', '657', 7, 7)
#
#    #def closeEvent(self, event):
#    #    self.close()

class MainWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle('Station Map Program (Version 0.1)')

        self.main_widget = QtGui.QWidget(self)
        #self.main_widget = StationMapWidget(self)
        l = QtGui.QVBoxLayout(self.main_widget)

        self.main_cnm_figure = StationMapWidget(self.main_widget,
                width=5, height=4, dpi=100)
        a = InterStaDist.from_file('staloc_BC', 4, 3, 1, 0)
        self.main_cnm_figure.set_InterStaDist(a)
        self.main_cnm_figure.plot_init()

        self.main_cnm_figure.connect()

        self.initUI()

        #l = QtGui.QVBoxLayout()

        l.addWidget(self.main_cnm_figure)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        #self.setLayout(l)

        self.statusBar().showMessage('Ready!', 2000)

    def initUI(self):
        self.main_cnm_figure.plot_sta_radius('025', '657', 7, 7)

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
    win.show()
    sys.exit(app.exec_())

