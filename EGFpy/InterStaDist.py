#!/usr/bin/env python
#-*-coding:utf-8-*-


#DEBUG = True
DEBUG = False
from collections import Iterable

import numpy as np

from obspy.geodetics import gps2dist_azimuth

class InterStaDist(object):
    def __init__(self, lons=None, lats=None, stas=None, nets=None):
        if isinstance(lons, list):
            lons = np.array(lons)
        if isinstance(lats, list):
            lats = np.array(lats)
        self.lons = lons
        self.lats = lats
        self.stas = stas
        self.nets = nets
        self.n = len(self.lons)

        self.init_dist_array()
        self.init_sta_index()

    def init_dist_array(self):
        if self.lons is None or self.lats is None:
            self.inter_sta_dist_array = None
            return

        assert(len(self.lons)==len(self.lats) and \
              len(self.lons) > 0)
        # build the inter-station distance array, distance in km 
        self.inter_sta_dist_array = np.zeros((self.n, self.n))

        for i in range(self.n):
            for j in range(i+1, self.n):
                dist = gps2dist_azimuth(self.lats[i], self.lons[i],
                                        self.lats[j], self.lons[j])[0]
                self.inter_sta_dist_array[i,j] = dist / 1000.0
                self.inter_sta_dist_array[j,i] = dist / 1000.0

    def init_sta_index(self):
        # build staname to index link
        if self.stas is None:
            self.sta_index = None
        else:
            self.sta_index = dict(zip(self.stas, range(self.n)))

    @classmethod
    def from_file(cls, staloc_file, lon_col=0, lat_col=1,
                  sta_col=None, net_col=None,
                  comment='#', sep=None):
        lons = []
        lats = []
        if sta_col is not None:
            stas = []
        else:
            stas = None
        if net_col is not None:
            nets = []
        else:
            stas = None
        with open(staloc_file, 'rU') as slf:
            for line in slf.readlines():
                if line.strip()[0] == '#':
                    continue

                if sep is None:
                    line_split = line.strip().split()
                else:
                    line_split = line.strip().split(sep)

                lon = float(line_split[lon_col])
                lons.append(lon)
                lat = float(line_split[lat_col])
                lats.append(lat)
                if sta_col is not None:
                    stas.append(line_split[sta_col])
                if net_col is not None:
                    nets.append(line_split[net_col])

        lons = np.array(lons)
        lats = np.array(lats)
        #self.stas = stas
        #self.nets = nets
        #self.n = len(self.lons)
        #self.init_dist_array()
        #self.init_sta_index()
        return cls(lons, lats, stas, nets)


    #------------------------------
    #        get method
    #------------------------------
    def get_sta_index(self, sta):
        if isinstance(sta, str):
            return self.sta_index[sta]
        else:
            indexs = [self.sta_index[s] for s in sta]
            return np.array(indexs)

    def get_sta(self, index):
        if self.stas is None:
            return None
        else:
            return self.stas[index]

    def get_lonlat(self, sta=None):
        if sta is None:
            return self.lons, self.lats
        else:
            sta_index = self.get_sta_index(sta)
            return self.lons[sta_index], self.lats[sta_index]

    def get_adjacent_sta(self, center_sta, search_radius, exclude_self=False):
        center_sta_index = self.get_sta_index(center_sta)
        if DEBUG:
            print('sta: {} sta_index: {} search radius: {}'.format(\
                   center_sta, center_sta_index, search_radius))

        if exclude_self:
            exclude_center_sta = np.ones(self.n, dtype=bool)
            exclude_center_sta[center_sta_index] = False
            adjacent_sta_indexs, = np.where(\
                 (self.inter_sta_dist_array[center_sta_index,:] < search_radius) &\
                 exclude_center_sta)
        else:
            adjacent_sta_indexs, = np.where(\
                 self.inter_sta_dist_array[center_sta_index,:] < search_radius)
        if DEBUG:
            print(self.inter_sta_dist_array[center_sta_index,:])
            print(adjacent_sta_indexs)

        adjacent_sta = [self.get_sta(index) \
              for index in adjacent_sta_indexs]

        return adjacent_sta

if __name__ == '__main__':
    DEBUG = True
    import matplotlib.pyplot as plt

    a = InterStaDist.from_file('staloc_BC', 4, 3, 1, 0)

    adj_sta = a.get_adjacent_sta('114', search_radius=8.0)

    print adj_sta


    plt.plot(*a.get_lonlat(), marker='v', mfc='none', linestyle='none')
    plt.plot(*a.get_lonlat(adj_sta), marker='o', mfc='none', linestyle='none')


    plt.gca().set_aspect(1.0)
    plt.show()





