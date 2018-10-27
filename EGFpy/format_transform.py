#!/usr/bin/env python
#-*-coding:utf-8-*-

from __future__ import print_function, division, \
        unicode_literals
import os
import sys
import numpy as np
from obspy import read

def sac2CF(sacfile):
    st = read(sacfile)
    tr = st[0]
    sachd = tr.stats.sac
    data = tr.data
    t = tr.times() + sachd.b

    assert(sachd.npts % 2 == 1)

    data_pos = data[sachd.npts//2:]
    data_neg = data[:sachd.npts//2+1]
    t_pos = np.arange(sachd.npts//2+1) * sachd.delta
    data_part = np.vstack((t_pos, data_neg[::-1], data_pos)).T

    header = np.array([[sachd.evlo, sachd.evla, 0.0],\
                       [sachd.stlo, sachd.stla, 0.0]])

    output = np.vstack((header, data_part))

    CF_filename = change_ext(sacfile)
    np.savetxt(CF_filename, output)

def change_ext(sacfile):
    sacfile_split = os.path.basename(sacfile).split('.')
    CF_filename = 'CF.' + '.'.join(sacfile_split[:-1])\
                   + '.dat'
    return CF_filename


if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.exit()
    for sacfile in sys.argv[1:]:
        sac2CF(sacfile)


