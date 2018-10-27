# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 23:28:58 2018

@author: xuyihe
"""
DEBUG=False
#DEBUG=False
import os
import sys
import re
from collections import OrderedDict
from ConfigParser import ConfigParser

import matplotlib.pyplot as plt
import matplotlib as mpl

from Disp import Disp
from GVImage import GrpVelImg#, readCF, taper_gv
#import numpy as np


def overlay_disp_EGF(Dispfile, EGFfile, **kwargs):
    gvi = GrpVelImg(EGFfile, **kwargs)
    disp = Disp(Dispfile)
    
    #plt.figure(figsize=(6,6))
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    
    plt.imshow(gvi.img, origin='lower', extent=gvi.extent, 
               cmap='jet', aspect=gvi.aspect)
    plt.plot(disp.T, disp.masked_disp, 'w', lw=2)
    
    plt.xlabel('Period (s)', fontsize=14)
    plt.ylabel('Velocity (km/s)', fontsize=14)
    
def file2list(filename):
    if not os.path.isfile(filename):
        return None
    with open(filename, 'rU') as f:
        filelist = [line.strip() for line in f.readlines()]
    return filelist

def list2file(filename, lst):
    with open(filename, 'w') as f:
        f.write('\n'.join(lst) + '\n')
    
def associate_disp_EGF_list(displist, EGFlist, 
                            re_pattern_disp=r'GDisp\.GFcn\.([^-]+)-([^_]+)_.*\.dat',
                            re_pattern_EGF=r'GFcn\.([^-]+)-([^_]+)_.*\.dat'):
    disp_dict = OrderedDict()
    EGF_dict = OrderedDict()
    associate_dict = OrderedDict()
    for dispf in displist:
        m_disp = re.match(re_pattern_disp, os.path.basename(dispf))
        if DEBUG:
            if m_disp:
                print('match:', re_pattern_disp, os.path.basename(dispf),
                      m_disp.groups())
            else:
                print('not match:', re_pattern_disp, os.path.basename(dispf))
                sys.exit()
        if m_disp:
            sta_pair = m_disp.group(1) + '-' + m_disp.group(2)
            disp_dict[sta_pair] = dispf

    for EGFf in EGFlist:
        m_EGF = re.match(re_pattern_EGF, os.path.basename(EGFf))
        if DEBUG:
            if m_EGF:
                print('match:', re_pattern_EGF, os.path.basename(EGFf),
                     m_EGF.groups())
            else:
                print('not match:', re_pattern_EGF, os.path.basename(EGFf))
                sys.exit()
        if m_EGF:
            sta_pair = m_EGF.group(1) + '-' + m_EGF.group(2)
            EGF_dict[sta_pair] = EGFf

    if DEBUG: 
        print disp_dict.items()[:10]
        print EGF_dict.items()[:10]
    #print disp_dict.keys()[:10]
    #sys.exit()    
    sorted_sta_pair = sorted(list(set(disp_dict.keys() + EGF_dict.keys())))
    if DEBUG:
        pass
        #print 'sorted sta pair: {}'.format(sorted_sta_pair)
    associate_dict = OrderedDict.fromkeys(sorted_sta_pair)
    
    for sta_pair in associate_dict:
        associate_dict[sta_pair] = dict()
        if sta_pair in disp_dict:
            associate_dict[sta_pair]['disp'] = disp_dict[sta_pair]
        else:
            associate_dict[sta_pair]['disp'] = None
        if sta_pair in EGF_dict:
            associate_dict[sta_pair]['EGF'] = EGF_dict[sta_pair]
        else:
            associate_dict[sta_pair]['EGF'] = None
        #print associate_dict['KMI-MC15']['disp']

    return associate_dict

def filter_sta(associate_dict, stalist):
    filtered_dict = OrderedDict()
    for sta_pair in associate_dict:
        sta1, sta2 = sta_pair.split('-')
        if sta1 in stalist and sta2 in stalist:
            if DEBUG:
                print sta1, sta2
            filtered_dict[sta_pair] = associate_dict[sta_pair]
    return filtered_dict

def associate_dict_from_config(configfile):
    config = ConfigParser()
    config.read(configfile)

    sacflist = file2list(config.get('input', 'cf_list_file'))
    displist = file2list(config.get('input', 'disp_list_file'))
    re_pattern_cf = config.get('input', 're_pattern_cf')
    re_pattern_disp = config.get('input', 're_pattern_disp')

    if DEBUG:
        print(re_pattern_disp)
        print(re_pattern_cf)

    asso_dict = associate_disp_EGF_list(displist, sacflist,
                re_pattern_disp, re_pattern_cf,
                #r'GDisp\.BA\.(\d+)-BA\.(\d+)\.dat',
                #r'BA\.(\d+)\.\d+\.BA\.(\d+)\..+\.s',
                )
    if DEBUG:
        print asso_dict.items()[:10]
    return asso_dict

def disp_EGF_list_associated(configfile):
    asso_dict = associate_dict_from_config(configfile)

    config = ConfigParser()
    config.read(configfile)
    stalist_BA = file2list(config.get('input', 'stalist'))

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

def find_sta_pair(sta_group1, sta_group2, sta_pair_list):
    sta_pair_selected = []
    for sta1 in sta_group1:
        for sta2 in sta_group2:
            sta12 = '{}-{}'.format(sta1,sta2)
            sta21 = '{}-{}'.format(sta2,sta1)
            if sta12 in sta_pair_list:
                sta_pair_selected.append(sta12)
            if sta21 in sta_pair_list:
                sta_pair_selected.append(sta21)
    return sta_pair_selected


if __name__ == '__main__':
#    import matplotlib.pyplot as plt
#    import matplotlib as mpl
#    EGFfile = r'EGFAnalysisTimeFreq_version_2015\EGFs\GFcn.KMI-MC01_10-50s_10Mon.dat'
#    Dispfile = r'EGFAnalysisTimeFreq_version_2015\Disper\GDisp.GFcn.KMI-MC01_10-50s_10Mon.dat'
#    Dispfile1 = r'EGFAnalysisTimeFreq_version_2015\tempdisp\GDisp.GFcn.KMI-MC01_10-50s_10Mon.dat'

    displist = file2list('Disp.lst')
    EGFlist = file2list('CF.lst')

    associate_dict = associate_disp_EGF_list(displist, EGFlist,
                re_pattern_disp=r'GDisp\.BA\.(\d+)-BA\.(\d+)\.dat',
                re_pattern_EGF=r'BA\.(\d+)\.\d+\.BA\.(\d+)\..+\.s',
                )

    associate_dict_config = associate_dict_from_config('cnm.cfg')
    #sys.exit()

    if DEBUG:
        adkeys = associate_dict.keys()
        adckeys = associate_dict_config.keys()
        print('len(adkeys), {}, len(adckeys), {}'.format(\
                len(adkeys), len(adckeys)))
        for i in range(len(adkeys)):
            if adkeys[i] != adckeys[i]:
                print('Key not equal: {} {} {}'.format(\
                    i, adkeys[i], adckeys[i]))
                #sys.exit()
        #pass


    plt.figure(figsize=(6,6))

    for sta_pair, disp_EGF_file in associate_dict.iteritems():
        if disp_EGF_file['disp'] is None or disp_EGF_file['EGF'] is None:
            continue
        print disp_EGF_file['disp'], disp_EGF_file['EGF']
        overlay_disp_EGF(disp_EGF_file['disp'], disp_EGF_file['EGF'],
                         cfile_format='SAC_CF', 
                        gvmin=0.5, gvmax=4, gvdelta=0.002,
                        Tmin=0.5, Tmax=5, Tdelta=0.1)
        plt.title(sta_pair)
        #figfilename = sta_pair+'.jpg'
        #plt.savefig(figfilename, dpi=100)
        #plt.clf()
        #print figfilename, 'created.'
        plt.show()
        sys.exit()
