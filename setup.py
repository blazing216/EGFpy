#!/usr/bin/env python
#-*-coding:utf-8-*-

from distutils.core import setup

setup(name='EGFpy',
      version='0.1',
      description='Python packages of validating Emperical Green' + \
        ' Functions in ambient noise tomograph',
      author='Yihe Xu',
      author_email='xuyihe@cea-igp.ac.cn',
      packages=['EGFpy'],
      scripts=['EGFpy/format_transform.py', 'EGFpy/plot_disp.py',
              'EGFpy/check_and_modify_disp_pyqt4.py',
              'EGFpy/make_file_list.py',
              'EGFpy/plot_median_disp.py', 'EGFpy/get_median_disp.py'],
     )
