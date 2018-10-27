#!/usr/bin/env python
#-*-coding:utf-8-*-

#DEBUG = True
DEBUG = False

import warnings

from copy import deepcopy

import os
import sys
from ConfigParser import ConfigParser
from find_associated_EGF_and_disp import (file2list,
    associate_disp_EGF_list, filter_sta)


from PyQt4 import QtCore, QtGui
from GVImageWidget import DispCheckModifyWidget

class GVImageSubWidget(DispCheckModifyWidget):
    def __init__(self, parent=None, **kwargs):
        DispCheckModifyWidget.__init__(self, parent, **kwargs)
        self.fig.canvas.mpl_connect('button_press_event', self.onpress_emit)

    def onpress_emit(self, event):
        self.emit(QtCore.SIGNAL('mouse_click_sub'), event)


class GVImageGroupWidget(QtGui.QWidget):
    def __init__(self, parent, sta_pair_list=None, disp_list=None,
                 EGF_list=None, sta_pair_state_list=None):
        QtGui.QWidget.__init__(self, parent)

        self.sta_pair_list = sta_pair_list
        self.disp_list = disp_list
        self.EGF_list = EGF_list
        self.sta_pair_state_list = sta_pair_state_list
        self.copy_ninja_state = False

        self.initUI()


    def initUI(self):
        if self.disp_list is None or self.EGF_list is None:
            self.gvigroup = None
            return
        else:
            assert(len(self.disp_list) == len(self.EGF_list))
            self.ndisp = len(self.disp_list)
            self.gvigroup = []

            self.main_widget = QtGui.QWidget(self)
            grid = QtGui.QGridLayout(self.main_widget)

            #updated_disp_list = self.__get_updated_disp_list()
            #updated_disp_list = deepcopy(self.disp_list)
            #for dispfile in update

            #self.sta_pair_state_dict = dict(zip(self.sta
            if DEBUG:
                print('GVImageGroupWidget: in initUI')
                print('self.disp_list = {}'.format(self.disp_list))

            for i, sta_pair, dispfile, EGFfile, state in zip(range(len(self.sta_pair_list)),
                                                      self.sta_pair_list,
                                                      self.disp_list,
                                                      #updated_disp_list,
                                                      self.EGF_list,
                                                      self.sta_pair_state_list):
                #gviwidget = DispCheckModifyWidget(self,
                #    width=5, height=4, dpi=100)
                if DEBUG:
                    print dispfile, EGFfile
                gviwidget = GVImageSubWidget(self,
                    width=3, height=3, dpi=100)
                gviwidget.set_gv_para(0.5, 4, 0.002)
                gviwidget.set_T_para(0.5, 5, 0.1)
                gviwidget.set_EGF_file(EGFfile)
                gviwidget.set_cfile_format('SAC_CF')
                gviwidget.set_Disp_file(dispfile)
                gviwidget.set_sta_pair(sta_pair)
                gviwidget.set_accept_state(state)
                gviwidget.plot_init()
                gviwidget.connect()

                gviwidget.setSizePolicy(QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)

                grid.addWidget(gviwidget, i//3, i%3)
                self.gvigroup.append(gviwidget)

            self.scroll = QtGui.QScrollArea()
            self.scroll.setWidget(self.main_widget)

            vbox = QtGui.QVBoxLayout()
            vbox.addWidget(self.scroll)


            #grid.addWidget(gviwidget)
            #for i in range(1):
                #self.gvigroup[i].setSizePolicy(QtGui.QSizePolicy.Expanding,
                #                               QtGui.QSizePolicy.Expanding)
                #grid.addWidget(self.gvigroup[i])

            #self.grid = QtGui.QGridLayout(self)
            #ncol = 3
            #for i in range(self.ndisp):
            #    self.gvigroup[i].setSizePolicy(QtGui.QSizePolicy.Expanding,
            #                                   QtGui.QSizePolicy.Expanding)
            #    self.grid.addWidget(self.gvigroup[i], i//3, i%3)

            #self.set_copy_ninja(self)
            self.setLayout(vbox)

        self.setSizePolicy(QtGui.QSizePolicy.Expanding,
                           QtGui.QSizePolicy.Expanding)

    #    for gviwidget in self.gvigroup:
            #gviwidget.selectChanged.connect(self.onselectchanged)
    #        gviwidget.stateChanged.connect(self.onstatechanged)

    #def onstatechanged(self, sta_pair, state):
        #self.sta_pair_state_dict[sta_pair] = state

    def updateUI(self):
        if self.disp_list is None or self.EGF_list is None:
            self.gvigroup = None
            return
        else:
            assert(len(self.disp_list) == len(self.EGF_list))
            self.ndisp = len(self.disp_list)
            self.gvigroup = []

            if self.main_widget is not None:
                self.main_widget.deleteLater()

            self.main_widget = QtGui.QWidget(self)
            grid = QtGui.QGridLayout(self.main_widget)

            #updated_disp_list = self.__get_updated_disp_list()
            for i, sta_pair, dispfile, EGFfile, state in zip(range(len(self.sta_pair_list)),
                                                      self.sta_pair_list,
                                                      #updated_disp_list,
                                                      self.disp_list,
                                                      self.EGF_list,
                                                      self.sta_pair_state_list):
                #gviwidget = DispCheckModifyWidget(self,
                #    width=5, height=4, dpi=100)
                if DEBUG:
                    print dispfile, EGFfile
                gviwidget = GVImageSubWidget(self,
                    width=3, height=3, dpi=100)
                gviwidget.set_gv_para(0.5, 4, 0.002)
                gviwidget.set_T_para(0.5, 5, 0.1)
                gviwidget.set_EGF_file(EGFfile)
                gviwidget.set_cfile_format('SAC_CF')
                gviwidget.set_Disp_file(dispfile)
                gviwidget.set_sta_pair(sta_pair)
                gviwidget.set_accept_state(state)
                gviwidget.plot_init()
                gviwidget.connect()

                gviwidget.setSizePolicy(QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)

                grid.addWidget(gviwidget, i//3, i%3)
                self.gvigroup.append(gviwidget)

            #scroll = QtGui.QScrollArea()
            self.scroll.setWidget(self.main_widget)

            #vbox = QtGui.QVBoxLayout()
            #vbox.addWidget(scroll)

            #self.setLayout(vbox)

        #self.setSizePolicy(QtGui.QSizePolicy.Expanding,
        #                   QtGui.QSizePolicy.Expanding)

    #def __get_updated_disp_list(self):
    #    updated_disp_list = []
    #    for dispfile in self.disp_list:
    #        filename = os.path.basename(dispfile)
    #        if filename in self.accepted_disp_dict.keys():
    #            updated_disp_list.append(\
    #                self.accepted_disp_dict[filename])
    #        else:
    #            updated_disp_list.append(\
    #                dispfile)
    #    return updated_disp_list

    def save_accepted_disp(self, accepted_disp_folder):
        for gviwidget in self.gvigroup:
            if gviwidget.get_accept_state() == 'A':
                gviwidget.save_disp(accepted_disp_folder)

    def set_copy_ninja(self, state=None):
        if state is not None:
            self.copy_ninja_state = state

        if self.gvigroup is None:
            return
        if self.copy_ninja_state == True:
            for gviwidget in self.gvigroup:
                self.connect(gviwidget, QtCore.SIGNAL('mouse_click_sub'),
                             self.groupclick)
        else:
            for gviwidget in self.gvigroup:
                self.disconnect(gviwidget, QtCore.SIGNAL('mouse_click_sub'),
                               self.groupclick)

    def set_disp_list(self, disp_list):
        self.disp_list = disp_list

    def set_accepted_disp_list(self, accepted_disp_list):
        self.accepted_disp_list = accepted_disp_list
        self.accepted_disp_dict = dict([\
            (os.path.basename(f), f) \
            for f in self.accepted_disp_list])

    def set_EGF_list(self, EGF_list):
        self.EGF_list = EGF_list

    def set_sta_pair_list(self, sta_pair_list):
        self.sta_pair_list = sta_pair_list
    #def copy_ninja(self, state):
    #    if self.gvigroup is None:
    #        return
    #    if state == True:
    #        for gviwidget in self.gvigroup:
    #            self.connect(gviwidget, QtCore.SIGNAL('mouse_click_sub'),
    #                         self.groupclick)
    #    else:
    #        for gviwidget in self.gvigroup:
    #            self.disconnect(gviwidget, QtCore.SIGNAL('mouse_click_sub'),
    #                           self.groupclick)

    def set_sta_pair_state_list(self, state_list):
        self.sta_pair_state_list = state_list

    def groupclick(self, event):
        if DEBUG:
            print('xdata={}, ydata={}'.format(event.xdata, event.ydata))

        for gviwidget in self.gvigroup:
            gviwidget.onclick(event)

    def onkeypress(self, event):
        if DEBUG:
            print('Key pressed: {}'.format(event.key()))
            print('event.key(): {}'.format(type(event.key())))
            #print('GVImage selected: {}'.format(\
                #self.main_cnm_figure.selected))

        if event.key() == QtCore.Qt.Key_F2:
            self.set_copy_ninja(True)
        elif event.key() == QtCore.Qt.Key_F3:
            self.set_copy_ninja(False)
        elif event.key() in [QtCore.Qt.Key_Control,
                             QtCore.Qt.Key_Shift,
                             QtCore.Qt.Key_Alt]:
            for gviwidget in self.gvigroup:
                gviwidget.onkeypress(event)
        else:
            for gviwidget in self.gvigroup:
                if gviwidget.selected:
                    gviwidget.onkeypress(event)

    def onkeyrelease(self, event):
        if DEBUG:
            print('Key pressed: {}'.format(event.key()))
        if event.key() in [QtCore.Qt.Key_Control,
                             QtCore.Qt.Key_Shift,
                             QtCore.Qt.Key_Alt]:
            for gviwidget in self.gvigroup:
                gviwidget.onkeyrelease(event)
        else:
            for gviwidget in self.gvigroup:
                if gviwidget.selected:
                    gviwidget.onkeyrelease(event)

    def enterEvent(self, event):
        self.setFocus()

def get_disp_EGF_asso_list():
    config = ConfigParser()
    config.read('cnm.cfg')
    
    sacflist = file2list(config.get('input', 'cf_list_file'))
    displist = file2list(config.get('input', 'disp_list_file'))
    
    asso_dict = associate_disp_EGF_list(displist, sacflist, 
                r'GDisp\.BA\.(\d+)-BA\.(\d+)\.dat',
                r'BA\.(\d+)\.\d+\.BA\.(\d+)\..+\.s',
                )
    
    stalist_BA = file2list('stalist.BA')
    asso_dict_BA = filter_sta(asso_dict, stalist_BA)
    
    sta_pair_list = []
    dispfilelist = []
    EGFfilelist = []
    for sta_pair, disp_EGF_pair in asso_dict_BA.iteritems():
        if disp_EGF_pair['disp'] is None or  disp_EGF_pair['EGF'] is None:
            continue
        sta_pair_list.append(sta_pair)
        dispfilelist.append(disp_EGF_pair['disp'])
        EGFfilelist.append(disp_EGF_pair['EGF'])

    return sta_pair_list, dispfilelist, EGFfilelist


class MainWindow(QtGui.QMainWindow):
    def __init__(self):
        try:
            QtGui.QMainWindow.__init__(self)
            self.setWindowTitle('Check and Modify Program (Version 0.5)')

            sta_pair_list, disp_list, EGF_list = get_disp_EGF_asso_list()
            if DEBUG:
                print(sta_pair_list[:10], disp_list[:10], EGF_list[:10])
            #self.main_widget = GVImageGroupWidget(self)
            #print sta_pair_list[:10]
            self.main_widget = GVImageGroupWidget(self,
                  sta_pair_list[:10], disp_list[:10], EGF_list[:10])

            #self.main_widget.set_disp_list(disp_list[:10])
            #self.main_widget.set_EGF_list(EGF_list[:10])
            #self.main_widget.setupUI()

            #l = QtGui.QVBoxLayout(self.main_widget)

            self.main_widget.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                           QtGui.QSizePolicy.Expanding)


            self.main_widget.setFocus()
            self.setCentralWidget(self.main_widget)

            self.statusBar().showMessage('Ready!', 2000)
        except RuntimeWarning:
            print('Found it')

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

        self.main_widget.onkeypress(event)

        #if event.key() == QtCore.Qt.Key_F2:
        #    self.main_widget.set_copy_ninja(True)
        #elif event.key() == QtCore.Qt.Key_F3:
        #    self.main_widget.set_copy_ninja(False)
        #elif event.key() in [QtCore.Qt.Key_Control,
        #                     QtCore.Qt.Key_Shift,
        #                     QtCore.Qt.Key_Alt]:
        #    for gviwidget in self.main_widget.gvigroup:
        #        gviwidget.onkeypress(event)
        #else:
        #    for gviwidget in self.main_widget.gvigroup:
        #        if gviwidget.selected:
        #            gviwidget.onkeypress(event)



    def keyReleaseEvent(self, event):
        if DEBUG:
            print('Key pressed: {}'.format(event.key()))

        self.main_widget.onkeyrelease(event)

        #if event.key() in [QtCore.Qt.Key_Control,
        #                     QtCore.Qt.Key_Shift,
        #                     QtCore.Qt.Key_Alt]:
        #    for gviwidget in self.main_widget.gvigroup:
        #        gviwidget.onkeyrelease(event)
        #else:
        #    for gviwidget in self.main_widget.gvigroup:
        #        if gviwidget.selected:
        #            gviwidget.onkeyrelease(event)

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    win = MainWindow()
    if DEBUG:
        print isinstance(win, QtGui.QWidget)
        print dir(win)
    win.show()
    sys.exit(app.exec_())

