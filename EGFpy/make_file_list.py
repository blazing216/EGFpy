#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 2018

@author: xuyihe
"""

from __future__ import division, print_function, unicode_literals
import os
import sys

#path = '/home/xuyh/Research/BinchuanArray/ncf/disper/' + \
#        'cluster_analysis/FTDisp/alldisp'

usage = 'Usage: python %s dir > demo.lst' % sys.argv[0]

if len(sys.argv) != 2:
    print(usage)
    sys.exit()

path = os.path.abspath(sys.argv[1])
path = unicode(path)
assert(isinstance(path, unicode))

for p, ds, fs in os.walk(path):
    if len(fs) == 0: continue

    for f in fs:
        line = os.path.join(p, f)
        print(line.encode('utf-8'))

#with open('CF.lst', 'w') as cf:
#    for p, ds, fs in os.walk(ur'D:\2018 Project I 江西实验\ANT\stack'):
#        if len(fs) > 0:
#            for f in fs:
#                line = os.path.join(p, f) + '\n'
#                cf.write(line.encode('gbk'))
#
##%%
#sacflist = ftpy.file2list('CF.lst')
#displist = ftpy.file2list('Disp.lst')
#
#asso_dict = ftpy.associate_disp_EGF_list(displist, sacflist, 
#            r'GDisp\.JX\.(\d+)-JX\.(\d+)\.dat',
#            r'JX\.(\d+)\.\d+\.JX\.(\d+)\..+\.sac',
#            )
##print asso_dict.items()[0]
##asso_dict = 
#
##demo_sac_file = r'EGF_winlink/JX.10921.01.JX.21305.01.00012/JX.10921.01.JX.21305.01.00012.2017.310.TO.2017.365.N.00035.BHZ.BHZ.Z.Z.stack.sac'
##%%
##print os.path.exists(demo_sac_file)
#import matplotlib
#matplotlib.use('agg')
#for sta_pair, disp_EGF_pair in asso_dict.iteritems():
#    #disp_EGF_pair = asso_dict.items()[0][1]
#    if disp_EGF_pair['disp'] is None or  disp_EGF_pair['EGF'] is None:
#        continue
#    figfile = sta_pair + '.jpg'
#    figpath = os.path.join('cmp_figs', figfile)
#    if os.path.exists(figpath):
#        continue
#    #fig = plt.figure(1, figsize=(6,6))
#    ftpy.overlay_disp_EGF(disp_EGF_pair['disp'], disp_EGF_pair['EGF'], 
#                      cfile_format='SAC_CF', gvmin=0.5, gvmax=4,
#                      Tmin=0.5, Tmax=5, Tdelta=0.1)
#    plt.title(sta_pair)
#    plt.savefig(figpath, dpi=100)
#    print figfile
#    plt.close('all')
#    #gc.collect()
#    #sys.exit()
#                      
#                
