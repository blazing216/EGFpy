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
import math
from ConfigParser import ConfigParser

import EGFpy as ftpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from PyQt4 import QtCore, QtGui

if DEBUG:
    import pdb

def safe_add(j, dj, n):
    assert(j>=0 and j<n)
    newj = j + dj
    if newj >= n: newj = n-1
    if newj < 0: newj = 0
    return newj

def search_peak(i0, j0, ImageData, step=3):
    n = ImageData.shape[0]
    if DEBUG:
        print n
    j = j0
    j1 = safe_add(j, step, n)
    if DEBUG:
        print('j is int: {}'.format(isinstance(j,int)))
        print('j1 is int: {}'.format(isinstance(j1,int)))
        print('i0 is int: {}'.format(isinstance(i0,int)))
        print('i0=%d, j=%d, j1=%d, i0=%d' % (i0,j,j1,i0))
    while ImageData[j,i0] < ImageData[j1,i0]:
        j = j1
        j1 = safe_add(j, step, n)
        if DEBUG:
            print('j is int: {}'.format(isinstance(j,int)))
            print('j1 is int: {}'.format(isinstance(j1,int)))
            print('i0 is int: {}'.format(isinstance(i0,int)))
            print('i0=%d, j=%d, j1=%d, i0=%d' % (i0,j,j1,i0))
    return j

def search_near_peak(i0, j0, ImageData, step=3):
    j1 = search_peak(i0, j0, ImageData, step)
    j2 = search_peak(i0, j0, ImageData, -step)
    if DEBUG:
        print('j1=%d, j2=%d' % (j1,j2))
    if ImageData[j1, i0] > ImageData[j2, i0]:
        return j1
    else:
        return j2

def track_peak_forward(i0, j0, ImageData, step=3):
    j = search_near_peak(i0, j0, ImageData, step)
    if DEBUG:
        print('Track points: %d %d' % (i0, j))

    m = ImageData.shape[1]
    track_points = np.zeros(m)
    track_points[i0] = j

    for i in range(i0+1,m):
        j1 = search_near_peak(i, j, ImageData, step)
        track_points[i] = j1
        if DEBUG:
            print('Track points: %d %d' % (i, j1))
        j = j1

    return track_points

def track_peak_backward(i0, j0, ImageData, step=3):
    j = search_near_peak(i0, j0, ImageData, step)
    if DEBUG:
        print('Track points: %d %d' % (i0, j))

    m = ImageData.shape[1]
    track_points = np.zeros(m)
    track_points[i0] = j

    for i in range(i0-1,-1,-1):
        j1 = search_near_peak(i, j, ImageData, step)
        track_points[i] = j1
        if DEBUG:
            print('Track points: %d %d' % (i, j1))
        j = j1

    return track_points


def track_peak(i0, j0, ImageData, step=3):
    j = search_near_peak(i0, j0, ImageData, step)
    if DEBUG:
        print('Track points: %d %d' % (i0, j))

    m = ImageData.shape[1]
    track_points = np.zeros(m)
    track_points[i0] = j

    for i in range(i0+1,m):
        j1 = search_near_peak(i, j, ImageData, step)
        track_points[i] = j1
        if DEBUG:
            print('Track points: %d %d' % (i, j1))
        j = j1

    j = int(track_points[i0])
    for i in range(i0-1,-1,-1):
        j1 = search_near_peak(i, j, ImageData, step)
        track_points[i] = j1
        if DEBUG:
            print('Track points: %d %d' % (i, j1))
        j = j1

    return track_points


class testWin:
    def __init__(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.gvi = ftpy.GrpVelImg('../demo_data/CF.BA.040-BA.116.sac',
                      cfile_format='SAC_CF',
                      gvmin   = 0.5,
                      gvmax   = 4,
                      gvdelta = 0.002,
                      Tmin    = 0.5,
                      Tmax    = 5,
                      Tdelta  = 0.1)

        self.disp = None

        mpl.rcParams['xtick.labelsize'] = 8
        mpl.rcParams['ytick.labelsize'] = 8

        if DEBUG:
            print self.gvi.img.shape
        gvimage = self.ax.imshow(self.gvi.img, origin='lower',
                   extent=self.gvi.extent,
                   cmap='jet', aspect=self.gvi.aspect,
                           )

        self.ax.set_xlabel('Period (s)', fontsize=10)
        self.ax.set_ylabel('Velocity (km/s)', fontsize=10)

        #plt.gca().set_aspect(0.01)

        self.fig.canvas.mpl_connect('button_press_event', self.onclick)

        plt.show()

    def onclick(self, event):
        if event.button == 1:
            x, y = event.xdata, event.ydata
            i0 = int((x-0.5)/0.1)
            j0 = int((y-0.5)/0.002)
            if DEBUG:
                print('i0=%d, j0=%d' % (i0, j0))
            track_points = track_peak(i0, j0, self.gvi.img)
            track_points = track_points * 0.002 + 0.5
            T = 0.5 + np.arange(46)* 0.1
            if DEBUG:
                print(T)
                print(track_points)
            if self.disp is None:
                self.disp, = self.ax.plot(T, track_points, '.')
                self.fig.canvas.draw()
            else:
                self.disp.set_data(T, track_points)
                self.fig.canvas.draw()


if __name__ == '__main__':
    testWin()

