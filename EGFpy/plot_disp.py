#!/usr/bin/env python
#-*-coding:utf-8-*-

from __future__ import print_function, division, unicode_literals
import sys

usage = 'Usage:'
if len(sys.argv) == 1:
    print(usage)
    sys.exit()
else:
    argc = len(sys.argv)
    options = dict(outpufig=None)
    args = []
    i = 1
    while i < argc:
        arg = sys.argv[i]
        if arg[0] == '-':
            if arg[1] == 'o':
                assert(i+1 < argc)
                options['outputfig'] = sys.argv[i+1]
                i += 1
            else:
                print(usage)
                sys.exit()
        else:
            args.append(arg)
        i += 1

from EGFpy.Disp import Disp, plot_dispfile_list, plot_disp_list
import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(6,4), tight_layout=True)
ax = fig.add_subplot(111)
disp_list = []
for dispfile in args:
    disp = Disp(dispfile)
    #mask = (disp.masked_disp.data >=4.0) | (disp.masked_disp.mask)
    #disp.masked_disp = np.ma.masked_array(disp.masked_disp.data,
    #                                      mask=mask)
    disp_list.append(disp)
plot_disp_list(disp_list, ax=ax)
plt.xlabel('Period (s)')
plt.ylabel('Velocity (km/s)')
plt.xlim([0.5, 5])
plt.ylim([0.5, 4])
if options['outputfig'] is not None:
    plt.savefig(options['outputfig'], dpi=300)
plt.show()

