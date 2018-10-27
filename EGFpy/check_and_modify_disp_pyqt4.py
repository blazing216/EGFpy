#!/usr/bin/env python
#-*-coding: utf-8-*-
"""

@author: xuyihe

History:
    Yihe Xu Mar 2 2018: Initial version is copyed from
            example_BinChuanExperiment/check_and_modify_disp.py
"""


from __future__ import print_function
DEBUG = True
import warnings

import os
import sys
from ConfigParser import ConfigParser
from copy import deepcopy

import EGFpy as ftpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from obspy.geodetics import gps2dist_azimuth

from PyQt4 import QtCore, QtGui

from EGFpy.find_associated_EGF_and_disp import (file2list,
    associate_disp_EGF_list, filter_sta, disp_EGF_list_associated,
    find_sta_pair)
from EGFpy.GVImageGroupWidget import GVImageGroupWidget
from EGFpy.InterStaDist import InterStaDist
from EGFpy.InterStaDistWidget import StationMapWidget

if DEBUG:
    import pdb

class MainWindow(QtGui.QMainWindow):
    def __init__(self, configfile):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle('Check and Modify Program (Version 1.0)')

        self.sta_pair_list, self.disp_list, self.EGF_list = \
                self.get_disp_EGF_list_associated(configfile)

        #self.skip = skip
        self.configfile = configfile
        self.config = ConfigParser()
        self.config.read(self.configfile)

        self.checked_list_file = self.config.get('output', 'disp_checked_list')
        self.accepted_disp_folder = self.config.get('output',
                                                    'accepted_disp_folder')
        if not os.path.isdir(self.accepted_disp_folder):
            os.mkdir(self.accepted_disp_folder)

        #self.disp_checked_dict = None
        self.__read_disp_checked_list()

        self.setupUI()

        self.start()

        self.statusBar().showMessage('Ready!', 2000)

    def __read_disp_checked_list(self):
        if not os.path.isfile(self.checked_list_file):
            sta_pair_list = deepcopy(self.sta_pair_list)
            sta_pair_state_list = ['P'] * len(sta_pair_list)
        else:
            with open(self.checked_list_file, 'rU') as clf:
                sta_pair_list = []
                sta_pair_state_list = []
                for line in clf.readlines():
                    line_split = line.strip().split()
                    sta_pair_list.append(line_split[0])
                    sta_pair_state_list.append(line_split[1])
            if len(sta_pair_list) == 0:
                sta_pair_list = deepcopy(self.sta_pair_list)
                sta_pair_state_list = ['P'] * len(sta_pair_list)

        self.disp_checked_dict = dict(zip(sta_pair_list, sta_pair_state_list))
        if DEBUG:
            pass
            #print(self.disp_checked_dict)

    def __write_disp_checked_list(self):
        sorted_checked_list = sorted(self.disp_checked_dict.items(), \
                                     key=lambda d:d[0])
        with open(self.checked_list_file, 'w') as clf:
            for checked_item in sorted_checked_list:
                clf.write('{} {}\n'.format(checked_item[0], checked_item[1]))
    def __get_disp_count(self, state_type=None):
        state_list = self.disp_checked_dict.values()
        if state_type is None:
            return len(state_list)
        else:
            return state_list.count(state_type)

    def setupUI(self):
        self.main_widget = QtGui.QWidget(self)
        hbox = QtGui.QHBoxLayout(self.main_widget)

        self.embed_widget = QtGui.QWidget(self)
        vbox = QtGui.QVBoxLayout(self.embed_widget)
        if DEBUG:
            print('reading stations')
        self.stamapwidget = StationMapWidget(self, width=3, height=3)
        self.stamapwidget.setSizePolicy(QtGui.QSizePolicy.Maximum,
                                        QtGui.QSizePolicy.Maximum)

        self.countlabel = QtGui.QLabel('Dispersions left: {}/{}'.format(\
                    self.__get_disp_count('P'), self.__get_disp_count()))
        self.savebutton = QtGui.QPushButton('Save')
        self.nextbutton = QtGui.QPushButton('Next')
        self.nextskipbutton = QtGui.QPushButton('Next (skip non-pending)')
        self.prevbutton = QtGui.QPushButton('Prev')
        self.radiusedit = QtGui.QLineEdit()
        self.radiusratiolabel = QtGui.QLabel()
        self.countlabel.setSizePolicy(QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)
        self.savebutton.setSizePolicy(QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)
        self.nextbutton.setSizePolicy(QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)
        self.nextskipbutton.setSizePolicy(QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)
        self.prevbutton.setSizePolicy(QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)
        self.radiusedit.setSizePolicy(QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)
        self.radiusratiolabel.setSizePolicy(QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)
        vbox.addWidget(self.stamapwidget)
        vbox.addWidget(self.countlabel)
        vbox.addWidget(self.savebutton)
        vbox.addWidget(self.nextbutton)
        vbox.addWidget(self.nextskipbutton)
        vbox.addWidget(self.prevbutton)
        vbox.addWidget(self.radiusedit)
        vbox.addWidget(self.radiusratiolabel)
        vbox.addStretch(1)

        #if DEBUG:
        #    print('reading disp and EGF')
        #sta_pair_list, disp_list, EGF_list = disp_EGF_list_associated(self.configfile)
        #if DEBUG:
        #    print('generating gvi group widget')
        self.gvigroupwidget = GVImageGroupWidget(self,
              #sta_pair_list[:10], disp_list[:10], EGF_list[:10],
                                                )

        #self.main_widget.setSizePolicy(QtGui.QSizePolicy.Expanding,
        #                                  QtGui.QSizePolicy.Expanding)

        hbox.addWidget(self.embed_widget)
        hbox.addWidget(self.gvigroupwidget)

        #if DEBUG:
        #    print('set size')

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        self.resize(800,600)

        #self.embed_widget = QtGui.QWidget(self)
        #vbox = QtGui.QVBoxLayout(self.embed_widget)
        #self.stamapwidget = StationMapWidget(self)
        #vbox.addWidget(self.stamapwidget)

        #sta_pair_list, disp_list, EGF_list = disp_EGF_list_associated('cnm.cfg')
        #self.gvigroupwidget = GVImageGroupWidget(self,
        #      sta_pair_list[:10], disp_list[:10], EGF_list[:10])
        #self.gvigroupwidget.setSizePolicy(QtGui.QSizePolicy.Expanding,
        #                                  QtGui.QSizePolicy.Expanding)

        #vbox.addWidget(self.gvigroupwidget)

        #self.embed_widget.setFocus()
        #self.setCentralWidget(self.embed_widget)
        #self.resize(800,600)

    def __sort_sta_pairs(self, sta_pair_list):
        sta_pair_list_sort = deepcopy(sta_pair_list)

        # move pending sta_pair front
        spl_pending = []
        spl_accept_reject = []
        for sp in sta_pair_list_sort:
            if self.disp_checked_dict[sp] == 'P':
                spl_pending.append(sp)
            else:
                spl_accept_reject.append(sp)
        sta_pair_list_sort = spl_pending + spl_accept_reject

        # move center_sta_pair to be the first element
        sta_pair_list_sort.pop(sta_pair_list_sort.index(self.center_sta_pair))
        sta_pair_list_sort.insert(0, self.center_sta_pair)


        return sta_pair_list_sort

    def __cal_intersta_dist(self, sta1, sta2):
        lon1, lat1 = self.stamapwidget.inter_sta_dist.get_lonlat(sta1)
        lon2, lat2 = self.stamapwidget.inter_sta_dist.get_lonlat(sta2)
        return gps2dist_azimuth(lat1, lon1, lat2, lon2)[0] / 1000.0

    def start(self):
        print('updating group velocity image...', end='')
        self.isd = InterStaDist.from_file(self.config.get('auxiliary',
                                                          'sta_loc_file'),
                                   4, 3, 1, 0)
        self.stamapwidget.set_InterStaDist(self.isd)
        self.stamapwidget.plot_init()

        self.sta_pair_list, self.disp_list, self.EGF_list = \
                disp_EGF_list_associated(self.configfile)
        self.disp_dict = dict(zip(self.sta_pair_list, self.disp_list))
        self.EGF_dict = dict(zip(self.sta_pair_list, self.EGF_list))
        self.sta_pair_dict = dict(zip(self.sta_pair_list,
                                      range(len(self.sta_pair_list))))

        self.__update_disp_dict(only_selected=False)
        # plot the first radius
        first_sta_pair = self.sta_pair_list[0]
        self.center_sta_pair = first_sta_pair
        #if self.skip is True:
        #    #sta_pair_done = self.center_sta_pair
        #    index = self.sta_pair_dict[self.center_sta_pair]
        #    for sta_pair in self.sta_pair_list[index+1:]:
        #        if self.disp_checked_dict[sta_pair] == 'P':
        #            break
        #    self.center_sta_pair = sta_pair

        sta1, sta2 = self.center_sta_pair.split('-')
        print(sta1, sta2) 

        # calculate distance between sta1 and sta2
        self.dist_sta12 = self.__cal_intersta_dist(sta1, sta2)
        radius1, radius2 = self.dist_sta12 / 6.0, self.dist_sta12 / 6.0 
        #radius1, radius2 = 5, 5
        sta_group1, sta_group2 = \
                self.stamapwidget.plot_sta_radius(sta1, sta2,
                                                  radius1, radius2)
        self.radiusedit.setText('%.2f' % radius1)
        self.radiusratiolabel.setText('%.3f' % (1.0/6.0))

        # plot Group Velocity Image
        self.sta_pair_selected = find_sta_pair(sta_group1, sta_group2,
                                          self.sta_pair_list)
        self.sta_pair_selected = self.__sort_sta_pairs(self.sta_pair_selected)
        self.disp_selected = [self.disp_dict[sta] for sta in self.sta_pair_selected]

        if DEBUG:
            print('disp selected:')
            print(self.disp_selected)
        #accepted_disp_selected = [self.accepted_disp_dict[sta] \
        #                          for sta in sta_pair_selected]
        self.EGF_selected = [self.EGF_dict[sta] for sta in  self.sta_pair_selected]
        self.state_selected = [self.disp_checked_dict[sta] for sta in \
                          self.sta_pair_selected]
        #self.__update_disp_selected()

        self.gvigroupwidget.set_sta_pair_list(self.sta_pair_selected)
        self.gvigroupwidget.set_disp_list(self.disp_selected)
        #self.gvigroupwidget.set_accepted_disp_list(disp_accepted)
        self.gvigroupwidget.set_EGF_list(self.EGF_selected)
        self.gvigroupwidget.set_sta_pair_state_list(self.state_selected)
        self.gvigroupwidget.initUI()

        for gviwidget in self.gvigroupwidget.gvigroup:
            gviwidget.selectChanged.connect(self.onselectchanged)
            gviwidget.stateChanged.connect(self.onstatechanged)

        # plot sta pair path
        self.sta_pair_path_dict = dict()
        for sta_pair in self.sta_pair_selected:
            sta1, sta2 = sta_pair.split('-')
            lon1, lat1 = self.stamapwidget.inter_sta_dist.get_lonlat(sta1)
            lon2, lat2 = self.stamapwidget.inter_sta_dist.get_lonlat(sta2)

            self.sta_pair_path_dict[sta_pair], = \
                    self.stamapwidget.ax.plot([lon1, lon2], [lat1, lat2],
                                              'lightgray', alpha=0)
            if DEBUG:
                print('{}:({},{}) ({},{})'.format(sta_pair, lon1, lat1,
                                              lon2, lat2))
        self.savebutton.clicked.connect(self.onclicksavebutton)
        self.nextbutton.clicked.connect(self.onclicknextbutton)
        self.nextskipbutton.clicked.connect(self.onclicknextskipbutton)
        self.prevbutton.clicked.connect(self.onclickprevbutton)
        self.radiusedit.returnPressed.connect(self.onradiuseditchanged)

        print('Done.')

        #print sta_group1, sta_group2
    def update(self, sta_pair, radius=None):
        print('updating group velocity image...', end='')
        self.statusBar().showMessage('Updating Group Velocity Image...')
        if radius is None:
            # plot the radius
            self.center_sta_pair = sta_pair
            sta1, sta2 = self.center_sta_pair.split('-')
            print(sta1, sta2)

            # calculate distance between sta1 and sta2
            self.dist_sta12 = self.__cal_intersta_dist(sta1, sta2)
            radius1, radius2 = self.dist_sta12 / 6.0, self.dist_sta12 / 6.0 
            #radius1, radius2 = 5, 5
            sta_group1, sta_group2 = \
                    self.stamapwidget.plot_sta_radius(sta1, sta2,
                                                      radius1, radius2)
            self.radiusedit.setText('%.2f' % radius1)
            self.radiusratiolabel.setText('%.3f' % (1.0/6.0))
            #radius1, radius2 = 5, 5
            #sta_group1, sta_group2 = \
            #        self.stamapwidget.plot_sta_radius(sta1, sta2,
            #                                          radius1, radius2)
        else:
            sta1, sta2 = self.center_sta_pair.split('-')
            print(sta1, sta2)
            sta_group1, sta_group2 = \
                    self.stamapwidget.plot_sta_radius(sta1, sta2,
                                                      radius, radius)

        # plot Group Velocity Image
        self.sta_pair_selected = find_sta_pair(sta_group1, sta_group2,
                                          self.sta_pair_list)
        self.sta_pair_selected = self.__sort_sta_pairs(self.sta_pair_selected)
        self.disp_selected = [self.disp_dict[sta] for sta in self.sta_pair_selected]

        self.EGF_selected = [self.EGF_dict[sta] for sta in  self.sta_pair_selected]
        self.state_selected = [self.disp_checked_dict[sta] for sta in \
                          self.sta_pair_selected]

        self.gvigroupwidget.set_sta_pair_list(self.sta_pair_selected)
        self.gvigroupwidget.set_disp_list(self.disp_selected)
        #self.gvigroupwidget.set_accepted_disp_list(disp_accepted)
        self.gvigroupwidget.set_EGF_list(self.EGF_selected)
        self.gvigroupwidget.set_sta_pair_state_list(self.state_selected)
        self.gvigroupwidget.updateUI()

        for gviwidget in self.gvigroupwidget.gvigroup:
            gviwidget.selectChanged.connect(self.onselectchanged)
            gviwidget.stateChanged.connect(self.onstatechanged)

        # plot sta pair path
        self.sta_pair_path_dict = dict()
        for sta_pair in self.sta_pair_selected:
            sta1, sta2 = sta_pair.split('-')
            lon1, lat1 = self.stamapwidget.inter_sta_dist.get_lonlat(sta1)
            lon2, lat2 = self.stamapwidget.inter_sta_dist.get_lonlat(sta2)

            self.sta_pair_path_dict[sta_pair], = \
                    self.stamapwidget.ax.plot([lon1, lon2], [lat1, lat2],
                                              'lightgray', alpha=0)
            if DEBUG:
                print('{}:({},{}) ({},{})'.format(sta_pair, lon1, lat1,
                                              lon2, lat2))
        #self.savebutton.clicked.connect(self.onclicksavebutton)
        #self.nextbutton.clicked.connect(self.onclicknextbutton)
        #self.prevbutton.clicked.connect(self.onclickprevbutton)
        #self.radiusedit.returnPressed.connect(self.onradiuseditchanged)
        print('Done.')
        self.statusBar().showMessage('Done', 2000)

    def __update_disp_dict(self, only_selected=True):
        if only_selected:
            for sta_pair in self.sta_pair_selected:
                if self.disp_checked_dict[sta_pair] == 'A':
                    disp_filename = os.path.basename(\
                        self.disp_dict[sta_pair])
                    modified_disp_file = os.path.join(\
                        self.accepted_disp_folder, disp_filename)
                    if os.path.isfile(modified_disp_file):
                        self.disp_dict[sta_pair] = modified_disp_file
        else:
            for sta_pair in self.sta_pair_list:
                if self.disp_checked_dict[sta_pair] == 'A':
                    disp_filename = os.path.basename(\
                        self.disp_dict[sta_pair])
                    modified_disp_file = os.path.join(\
                        self.accepted_disp_folder, disp_filename)
                    if os.path.isfile(modified_disp_file):
                        self.disp_dict[sta_pair] = modified_disp_file


    def onstatechanged(self, sta_pair, state):
        if DEBUG:
            print('{}:{}'.format(sta_pair, state))
        self.disp_checked_dict[sta_pair] = state

    def onselectchanged(self, sta_pair, selected):
        if DEBUG:
            print('onselectchanged (receiver): {} {}'.format(\
                    sta_pair, selected))
        if selected:
            self.sta_pair_path_dict[sta_pair].set_alpha(1)
        else:
            self.sta_pair_path_dict[sta_pair].set_alpha(0)
        self.stamapwidget.fig.canvas.draw()

    def onclicksavebutton(self):
        self.__write_disp_checked_list()
        self.__save_accepted_disp()
        self.__update_disp_dict()
        self.countlabel.setText('Dispersions left: {}/{}'.format(\
                self.__get_disp_count('P'), self.__get_disp_count()))
        print('Save')

    def __save_accepted_disp(self):
        self.gvigroupwidget.save_accepted_disp(self.accepted_disp_folder)
        #for sta_pair in self.sta_pair_selected:
        #    if self.disp_checked_dict[sta_pair] == 'A':
        #        dispfn = os.path.basename(self.disp_dict[sta_pair])
        #        o
        #        
        #        self.accepted_disp_folder

        #
        #sorted_checked_list = sorted(self.disp_checked_dict.items(), \
        #                             key=lambda d:d[0])
        #with open(self.checked_list_file, 'w') as clf:
        #    for checked_item in sorted_checked_list:
        #        clf.write('{} {}\n'.format(checked_item[0], checked_item[1]))

    def onclicknextbutton(self):
        self.onclicksavebutton()
        if DEBUG:
            print('Next click')
        sta_pair_done = self.center_sta_pair
        index = self.sta_pair_dict[self.center_sta_pair]
        if index == len(self.sta_pair_list) - 1:
            return
        else:
            sta_pair = self.sta_pair_list[index+1]
            self.update(sta_pair)

    def onclicknextskipbutton(self):
        self.onclicksavebutton()
        if DEBUG:
            print('Next skip click')
        sta_pair_done = self.center_sta_pair
        index = self.sta_pair_dict[self.center_sta_pair]
        if index == len(self.sta_pair_list) - 1:
            return
        else:
            for sta_pair in self.sta_pair_list[index+1:]:
                if self.disp_checked_dict[sta_pair] == 'P':
                    break
            if DEBUG:
                print('update on next click: {}'.format(sta_pair))

            self.update(sta_pair)

    def onclickprevbutton(self):
        self.onclicksavebutton()
        print('Prev click')
        sta_pair_done = self.center_sta_pair
        index = self.sta_pair_dict[self.center_sta_pair]
        if index == 0:
            return
        else:
            sta_pair = self.sta_pair_list[index-1]
            self.update(sta_pair)
        #for sta_pair in self.sta_pair_list[index+1:]:
        #    if self.disp_checked_dict[sta_pair] == 'P':
        #        break
        #if DEBUG:
        #    print('update on next click: {}'.format(sta_pair))


    def onradiuseditchanged(self):
        print('Radius changed')
        radius = float(self.radiusedit.text())
        ratio = radius / self.dist_sta12
        self.radiusratiolabel.setText('%.3f' % ratio)
        self.update(self.center_sta_pair, radius)


    def closeEvent(self, event):
        self.onclicksavebutton()
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
        self.gvigroupwidget.onkeypress(event)

        #if event.key() == QtCore.Qt.Key_F2:
        #    self.gvigroupwidget.set_copy_ninja(True)
        #elif event.key() == QtCore.Qt.Key_F3:
        #    self.gvigroupwidget.set_copy_ninja(False)
        #elif event.key() in [QtCore.Qt.Key_Control,
        #                     QtCore.Qt.Key_Shift,
        #                     QtCore.Qt.Key_Alt]:
        #    for gviwidget in self.gvigroupwidget.gvigroup:
        #        gviwidget.onkeypress(event)
        #else:
        #    for gviwidget in self.gvigroupwidget.gvigroup:
        #        if gviwidget.selected:
        #            gviwidget.onkeypress(event)

    def keyReleaseEvent(self, event):
        if DEBUG:
            print('Key pressed: {}'.format(event.key()))

        self.gvigroupwidget.onkeyrelease(event)

        #if event.key() in [QtCore.Qt.Key_Control,
        #                     QtCore.Qt.Key_Shift,
        #                     QtCore.Qt.Key_Alt]:
        #    for gviwidget in self.gvigroupwidget.gvigroup:
        #        gviwidget.onkeyrelease(event)
        #else:
        #    for gviwidget in self.gvigroupwidget.gvigroup:
        #        if gviwidget.selected:
        #            gviwidget.onkeyrelease(event)

    def get_disp_EGF_list_associated(self, configfile):
        return disp_EGF_list_associated(configfile)
        #config = ConfigParser()
        ##config.read('cnm.cfg')

        #sacflist = file2list(config.get('input', 'cf_list_file'))
        #displist = file2list(config.get('input', 'disp_list_file'))

        #asso_dict = associate_disp_EGF_list(displist, sacflist, 
        #            r'GDisp\.BA\.(\d+)-BA\.(\d+)\.dat',
        #            r'BA\.(\d+)\.\d+\.BA\.(\d+)\..+\.s',
        #            )

        #stalist_BA = file2list('stalist.BA')
        #asso_dict_BA = filter_sta(asso_dict, stalist_BA)

        #sta_pair_list = []
        #dispfilelist = []
        #EGFfilelist = []
        #for sta_pair, disp_EGF_pair in asso_dict_BA.iteritems():
        #    if disp_EGF_pair['disp'] is None or  disp_EGF_pair['EGF'] is None:
        #        continue
        #    sta_pair_list.append(sta_pair)
        #    dispfilelist.append(disp_EGF_pair['disp'])
        #    EGFfilelist.append(disp_EGF_pair['EGF'])

        #return dispfilelist, EGFfilelist

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    win = MainWindow(sys.argv[1])
    win.show()
    sys.exit(app.exec_())

#class DispCheckModify(object):
#    def __init__(self, sta_pair_list, dispfilelist, EGFfilelist, 
#                 config, skip_checked=False):
#        assert(len(dispfilelist) == len(EGFfilelist) and \
#            len(dispfilelist) > 0)
#
#        #assert(os.path.isfile('cnm.cfg'))
#
#        #self.config = ConfigParser()
#        #self.config.read('cnm.cfg')
#        self.config = config
#
#        self.newdispfolder = self.config.get('output',
#                                             'checked_disp_folder')
#        self.skip_checked = skip_checked
#        self.index = 0
#        self.max_index = len(dispfilelist) - 1
#
#        self.sta_pair_list = sta_pair_list
#        self.dispfilelist = dispfilelist
#        self.EGFfilelist = EGFfilelist
#
#        self.dispfile_checked_list = ftpy.file2list(\
#            self.config.get('output', 'disp_checked_list'))
#        if self.dispfile_checked_list is None:
#            self.dispfile_checked_list = []
#
#        staloc_file =  self.config.get('auxiliary', 'sta_loc_file')
#        if staloc_file is None:
#            self.staloc = np.array([[np.nan], [np.nan]])
#        else:
#            self.staloc = []
#            with open(staloc_file, 'rU') as slf:
#                for line in slf.readlines():
#                    if line[0] == '#' or line.strip().split()[0]=='BP':
#                        continue
#                    lat, lon = line.strip().split()[3:5]
#                    lat = float(lat)
#                    lon = float(lon)
#                    self.staloc.append([lon, lat])
#            self.staloc = np.array(self.staloc)
#
#        self.fig, self.axs = plt.subplots(1,2,figsize=(8,4))
#        self.ax = self.axs[0]
#        self.mapax = self.axs[1]
#
#    def draw(self):
##        self.dispfile = self.dispfilelist[self.index]
##        self.EGFfile = self.EGFfilelist[self.index]
##        
##        self.disp = ftpy.Disp(self.dispfile)
#        self.dispfile = self.dispfilelist[self.index]
#        self.EGFfile = self.EGFfilelist[self.index]
#        print self.dispfile, self.EGFfile
#        self.disp = ftpy.Disp(self.dispfile)
#        self.newdispfile = os.path.join(self.newdispfolder, 
#                                    self.disp.filename)
#        if self.skip_checked:
#            #while os.path.isfile(self.newdispfile) or \
#            #        self.newdispfile in self.dispfile_checked_list:
#            while self.newdispfile in self.dispfile_checked_list:
#                print self.newdispfile, 'skipped.'
#                self.index += 1
#                
#                self.dispfile = self.dispfilelist[self.index]
#                self.EGFfile = self.EGFfilelist[self.index]
#                self.disp = ftpy.Disp(self.dispfile)
#                self.newdispfile = os.path.join(self.newdispfolder, 
#                                            self.disp.filename)
#        
#        self.gvi = ftpy.GrpVelImg(self.EGFfile, 
#                      cfile_format=self.config.get('input', 'cf_format'), 
#                      gvmin=0.5, gvmax=4,
#                      Tmin=0.5, Tmax=5, Tdelta=0.1)
#        self.bline = None
#        self.b = self.gvi.T[0]
#        self.e = self.gvi.T[-1]
#        self.eline = None
#        self.accept = True
#        
#        self.modify_state = False
#        
#        
#        self.ax.cla()
#        mpl.rcParams['xtick.labelsize'] = 12
#        mpl.rcParams['ytick.labelsize'] = 12
#        
#        self.ax.imshow(self.gvi.img, origin='lower', extent=self.gvi.extent, 
#                   cmap='jet', aspect=self.gvi.aspect)
#                   
#        self.mask_origin = self.disp.masked_disp.mask.copy()
#        self.displine_origin = self.ax.plot(self.disp.T, self.disp.masked_disp,
#                                        '.', color='gray', lw=2)
#        self.displine, = self.ax.plot(self.disp.T, self.disp.masked_disp, '.w', lw=2)
#        
#        self.ax.set_xlabel('Period (s)', fontsize=14)
#        self.ax.set_ylabel('Velocity (km/s)', fontsize=14)
#        
#        self.info_text = self.ax.text(0,1, sta_pair_list[self.index] + '\n' +\
#            'interstation distance: %.1f km \n' % self.disp.dist + \
#            'n: next disp  b: previous \n' + \
#            'left_click: min T accepted  right_click: max T \n' + \
#            'ctrl+left_click: change dispersion point manually\n',
#            ha='left', va='bottom', transform=self.ax.transAxes,
#            fontsize=8)
#            
#        
#        self.accept_text = self.ax.text( 1,1, 'Accept', ha='right', va='bottom',
#                                    transform=self.ax.transAxes)
#    #    fig = plt.figure(1, figsize=(6,6))
#    #    ftpy.overlay_disp_EGF(disp_EGF_pair['disp'], disp_EGF_pair['EGF'], 
#    #                      cfile_format='SAC_CF', gvmin=0.5, gvmax=4,
#    #                      Tmin=0.5, Tmax=5, Tdelta=0.1)
#        #plt.title(self.sta_pair_list[self.index])
#
#        self.mapax.cla()
#        self.mapax.plot(self.staloc[:,0], self.staloc[:,1], 'v', mfc='none',
#                        mec='black')
#        self.mapax.plot([self.disp.lon1, self.disp.lon2], 
#                        [self.disp.lat1, self.disp.lat2], '-v',
#                        mfc='lightgreen', mec='black', color='black')
#        self.fig.canvas.draw()
#    
#    def start(self):
#        
#        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
#        self.fig.canvas.mpl_connect('key_press_event', self.onkeypress)
#        self.fig.canvas.mpl_connect('key_release_event', self.onkeyrelease)
#        plt.show()
#     
#        
#    def onclick(self, event):
#        if event.button == 1:
#            if self.modify_state:
#                x, y = event.xdata, event.ydata
#                xi = np.argmin(np.abs(x - self.disp.T))
#                self.disp.masked_disp.data[xi] = y
#                
#                #self.displine.set_data(self.disp.T, self.disp.masked_disp)
#            else:
#                self.b = event.xdata
#                if self.b is None: # out of axes
#                    self.b = self.gvi.T[0]
#                if self.bline is None:
#                    self.bline, = self.ax.plot([self.b, self.b], [0, 1], 
#                                         'w', lw=2, 
#                                         transform=self.ax.get_xaxis_transform())
#                else:
#                    self.bline.set_data([self.b,self.b], [0,1])
#        elif event.button == 3:
#            self.e = event.xdata
#            if self.e is None: # out of axes
#                self.e = self.gvi.T[-1]
#            if self.eline is None:
#                self.eline, = self.ax.plot([self.e, self.e], [0, 1], 
#                                     'y', lw=2, 
#                                     transform=self.ax.get_xaxis_transform())
#            else:
#                self.eline.set_data([self.e,self.e], [0,1])
#                
#        self.newmask = (self.disp.T < self.b) | (self.disp.T > self.e) \
#                | self.mask_origin
##        print self.disp.T
##        print self.b, self.e
##        print self.mask_origin
##        print self.newmask
#        self.disp.masked_disp.mask = self.newmask
##        print self.disp.masked_disp.mask
#        self.displine.set_data(self.disp.T, self.disp.masked_disp)
#        
#        self.fig.canvas.draw()
#    
#    def onkeypress(self, event):
#        if event.key == 'r':
#            self.accept = False
#            self.accept_text.set_text('Reject')
#            self.fig.canvas.draw()
#        elif event.key == 'a':
#            self.accept = True
#            self.accept_text.set_text('Accept')
#            self.fig.canvas.draw()
#        elif event.key == 'n':
#            self.dispfile_checked_list.append(self.newdispfile)
#            ftpy.list2file(self.config.get('output', 'disp_checked_list'), 
#                           self.dispfile_checked_list)
#            if self.accept:
#                self.disp.save(self.newdispfile)
#                
#                print self.newdispfile, 'created.'
#                
#            self.index += 1
#            if self.index <= self.max_index:
#                print self.index, self.max_index
##                    self.dispfile = self.dispfilelist[self.index]
##                    self.EGFfile = self.EGFfilelist[self.index]
#                self.draw()
#            else:
#                plt.close(self.fig)
#        elif event.key == 'b':
#            self.index -= 1
#            if self.index < 0:
#                self.index += 1
#            else:
#                print self.dispfile_checked_list[-1]
#                self.dispfile_checked_list.pop()
#                self.draw()
#        elif event.key == 'control':
#            self.modify_state = True
#
#    def onkeyrelease(self, event):
#        if event.key == 'control':
#            self.modify_state = False
#
## back up code    
#    #config = ConfigParser()
#    #config.read('cnm.cfg')
#    #
#    #sacflist = ftpy.file2list(config.get('input', 'cf_list_file'))
#    #displist = ftpy.file2list(config.get('input', 'disp_list_file'))
#    #
#    #asso_dict = ftpy.associate_disp_EGF_list(displist, sacflist, 
#    #            r'GDisp\.BA\.(\d+)-BA\.(\d+)\.dat',
#    #            r'BA\.(\d+)\.\d+\.BA\.(\d+)\..+\.s',
#    #            )
#    #
#    #stalist_BA = ftpy.file2list('stalist.BA')
#    #asso_dict_BA = ftpy.filter_sta(asso_dict, stalist_BA)
#    #
#    #sta_pair_list = []
#    #dispfilelist = []
#    #EGFfilelist = []
#    #for sta_pair, disp_EGF_pair in asso_dict_BA.iteritems():
#    #    if disp_EGF_pair['disp'] is None or  disp_EGF_pair['EGF'] is None:
#    #        continue
#    #    sta_pair_list.append(sta_pair)
#    #    dispfilelist.append(disp_EGF_pair['disp'])
#    #    EGFfilelist.append(disp_EGF_pair['EGF'])
#    #prog = DispCheckModify(sta_pair_list, dispfilelist, EGFfilelist,
#    #                       skip_checked=True, config=config)
#    #
#    #
#    #
#    #prog.draw()
#    #prog.start()
#
#
#class MainWindow(QtGui.QMainWindow):
#    def __init__(self):
#        QtGui.QMainWindow.__init__(self)
#        self.setWindowTitle('Check and Modify Program (Version 0.1)')
#
#    def closeEvent(self, ce):
#        self.close()

