# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 22:42:51 2018

Provide both procedure and OOP style.

Keep the function of one 'function' as tight as possible

@author: xuyihe
"""
import os
import numpy as np
from obspy.geodetics import gps2dist_azimuth

class Disp(object):
    def __init__(self, disp_file):
        data, self.lon1, self.lat1, self.lon2, self.lat2 = \
            readDisp(disp_file)
        self.filename = os.path.basename(disp_file)
        self.T = data[:,0]
        self.disp = data[:,1]
        self.err = data[:,2]
        self.valid = data[:,3]
        self.dist = gps2dist_azimuth(self.lat1, self.lon1,
                                     self.lat2, self.lon2)[0]/1000.0
        mask = self.valid == 0
        self.masked_disp = np.ma.masked_array(self.disp, mask=mask)

    def save(self, newdispfile):
        writeDisp(newdispfile, self.lon1, self.lat1, self.lon2, self.lat2,
                  self.T, self.masked_disp.data, self.err, 
                  ~self.masked_disp.mask)

    def plot(self, ax=None, *args, **kwargs):
        ax = plt.gca() if ax is None else ax
        ax.plot(self.T, self.masked_disp, *args, **kwargs)

def readDisp(disp_file):
    with open(disp_file, 'rU') as df:
        lines = df.readlines()
        lon1, lat1 = lines[0].strip().split()[:2]
        lon2, lat2 = lines[1].strip().split()[:2]
        lon1, lat1, lon2, lat2 = float(lon1), float(lat1), \
                float(lon2), float(lat2)
    data = np.loadtxt(disp_file, usecols=range(4), skiprows=2)
    return data, lon1, lat1, lon2, lat2

def writeDisp(dispfile, lon1, lat1, lon2, lat2, T, disp, err, valid):
    with open(dispfile, 'w') as df:
        df.write('%10.4f %10.4f\n%10.4f %10.4f\n' % (lon1, lat1, lon2, lat2))
        for Ti, dispi, erri, validi in zip(T, disp, err, valid):
            df.write('%10.4f %10.4f %10.4f %2d\n' % (Ti, dispi, erri, validi))

def plot_disp_list(disp_list, ax=None):
    import matplotlib.pyplot as plt
    if ax is None:
        ax = plt.gca()

    if len(disp_list) == 0:
        return

    n = len(disp_list)
    m = len(disp_list[0].T)
    disp_array = np.zeros((n,m))
    mask_array = np.zeros((n,m), dtype=bool)

    for i, disp in enumerate(disp_list):
        disp_array[i,:] = disp.masked_disp.data
        mask_array[i,:] = disp.masked_disp.mask
    masked_disp_array = np.ma.masked_array(disp_array, mask=mask_array)

    T = disp_list[0].T

    ax.plot(T, masked_disp_array.T, 'k', lw=0.5, alpha=0.4)

def plot_dispfile_list(dispfile_list, ax=None):
    import matplotlib.pyplot as plt
    disp_list = []
    for dispfile in dispfile_list:
        disp_list.append(Disp(dispfile))

    plot_disp_list(disp_list, ax)

def get_disp_list(dispfile_list):
    disp_list = []
    for dispfile in dispfile_list:
        disp_list.append(Disp(dispfile))
    return disp_list

def disp2nparray(disp_list):
    if len(disp_list) == 0:
        return

    n = len(disp_list)
    m = len(disp_list[0].T)
    disp_array = np.zeros((n,m))
    mask_array = np.zeros((n,m), dtype=bool)

    for i, disp in enumerate(disp_list):
        disp_array[i,:] = disp.masked_disp.data
        mask_array[i,:] = disp.masked_disp.mask
    masked_disp_array = np.ma.masked_array(disp_array, mask=mask_array)

    T = disp_list[0].T

    return masked_disp_array, T

def median(masked_array, axis=-1):
    shape = masked_array.shape
    assert(len(shape) == 2 and shape >= axis)

    if axis == 0:
        median_array = np.zeros(shape[1])
        for i in range(shape[1]):
            col = masked_array[:,i].flatten()
            data = col.data[~col.mask]
            if len(data) == 0:
                median_array[i] = np.nan
            else:
                #print median_array
                median_array[i] = np.median(data)
        median_array.shape = 1, shape[1]
    elif axis == 1:
        median_array = np.zeros(shape[0])
        for i in range(shape[0]):
            col = masked_array[i,:].flatten()
            data = col.data[~col.mask]
            if len(data) == 0:
                median_array[i] = np.nan
            else:
                median_array[i] = np.median(data)
        median_array.shape = shape[0], 1

    return median_array

def get_disp_array(disp_file_list):
    disp_list = get_disp_list(disp_file_list)
    disp_array, T = disp2nparray(disp_list)
    return disp_array, T

def get_median_disp(disp_file_list):
    disp_list = get_disp_list(disp_file_list)
    disp_array, T = disp2nparray(disp_list)
    median_disp = median(disp_array, axis=0)
    return T, median_disp.flatten()


if __name__ == '__main__':
    # test1
    #demo_Disp_file = 'EGFAnalysisTimeFreq_version_2015\Disper\GDisp.GFcn.KMI-MC01_10-50s_10Mon.dat'
    #import matplotlib.pyplot as plt
    #
    #with open(demo_Disp_file, 'rU') as df:
    #    lines = df.readlines()
    #    lon1, lat1 = lines[0].strip().split()[:2]
    #    lon2, lat2 = lines[1].strip().split()[:2]
    #    lon1, lat1, lon2, lat2 = float(lon1), float(lat1), \
    #            float(lon2), float(lat2)
    #data = np.loadtxt(demo_Disp_file, usecols=range(4), skiprows=2)
    #

    #disp = Disp(demo_Disp_file)    
    #
    #plt.plot(data[:,0], data[:,1], 'k+')
    #disp.plot(plt.gca(), 'ko', mfc='none')
    #plt.xlabel('Period (s)')
    #plt.ylabel('Velocity (km/s)')
    #plt.ylim(2,5)
    #plt.show()

    # test2
    import sys
    import matplotlib.pyplot as plt

    disp_array, T = get_disp_array(sys.argv[1:])
    median_disp = median(disp_array, axis=0)

    plt.plot(T, disp_array.T, 'k', lw=0.5, alpha=0.1)
    plt.plot(T, median_disp.flatten(), 'r', lw=1)
    plt.xlim(0.5,5)
    plt.ylim(0.5,4)
    plt.show()

