#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Created on Wed Apr 18 22:42:51 2018

Provide both procedure and OOP style.

Keep the function of one 'function' as tight as possible

@author: xuyihe
"""

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

import matplotlib.pyplot as plt
from EGFpy.Disp import get_disp_array, median


fig = plt.figure(figsize=(6,4), tight_layout=True)
ax = fig.add_subplot(111)
disp_array, T = get_disp_array(args)
median_disp = median(disp_array, axis=0)

ax.plot(T, disp_array.T, 'k', lw=0.5, alpha=0.1)
ax.plot(T, median_disp.flatten(), 'r', lw=1, label='median')
plt.xlim(0.5,5)
plt.ylim(0.5,4)
plt.xlabel('Period (s)')
plt.ylabel('Velocity (km)')
plt.legend()
if options['outputfig'] is not None:
    plt.savefig(options['outputfig'], dpi=300)
plt.show()

