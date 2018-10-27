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

import numpy as np
import matplotlib.pyplot as plt
from EGFpy.Disp import get_disp_array, median


disp_array, T = get_disp_array(sys.argv[1:])
median_disp = median(disp_array, axis=0)

for t, d in zip(T, median_disp.flatten()):
    if np.isnan(d):
        print('%.2f %.2f' % (t, 0))
    else:
        print('%.2f %.2f' % (t, d))


