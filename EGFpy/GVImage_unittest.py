# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 13:25:27 2018

@author: xuyihe
"""

import unittest
import numpy as np

from GVImage import cos_taper

class TestGVImage(unittest.TestCase):
    def test_cos_taper(self):
        N = 10
        t = np.linspace(0, np.pi/2, N)
        f1 = np.ones_like(t)
        f1_cos_taper = cos_taper(f1, 0, N-1)
        f1_true_result = np.cos(t)
        self.assertTrue(np.allclose(f1_cos_taper, f1_true_result))
        
    def test_cos_taper_inv(self):
        N = 10
        t = np.linspace(0, np.pi/2, N)
        f1 = np.ones_like(t)
        f1_cos_taper = cos_taper(f1, N-1, 0)
        f1_true_result = np.cos(t[::-1])
        self.assertTrue(np.allclose(f1_cos_taper, f1_true_result))
        
    def test_cos_taper_assert(self):
        N = 10
        t = np.linspace(0, np.pi/2, N)
        f1 = np.ones_like(t)
        with self.assertRaises(AssertionError):
            cos_taper(f1, 0, 0)
        with self.assertRaises(AssertionError):
            cos_taper(f1, N-1, N-1)
        with self.assertRaises(AssertionError):
            cos_taper(f1, 0, N)
        with self.assertRaises(AssertionError):
            cos_taper(f1, -1, N-1)
            
    def test_cos_taper_part_taper(self):
        N = 100
        n1 = np.random.randint(0, N-1)
        M = np.random.randint(1, N-1)
        n2 = min(n1 + M, N-1)
        
        f1 = np.random.rand(N)
        f1_cos_taper = cos_taper(f1, n1, n2)
        f1_cos_taper_manual = f1.copy()
        f1_cos_taper_manual[n1:n2+1] *= \
            np.cos(np.arange(n2-n1+1, dtype=float)/(n2-n1)*np.pi/2)
        f1_cos_taper_manual[n2+1:] = 0.0
        self.assertTrue(np.allclose(f1_cos_taper, f1_cos_taper_manual))
        
    def test_cos_taper_part_taper_inv(self):
        N = 100
        n1 = np.random.randint(1, N-1)
        M = np.random.randint(1, N-1)
        n2 = max(n1 - M, 0)
        
        f1 = np.random.rand(N)
        f1_cos_taper = cos_taper(f1, n1, n2)
        f1_cos_taper_manual = f1.copy()
        f1_cos_taper_manual[n2:n1+1] *= \
            np.sin(np.arange(n1-n2+1, dtype=float)/(n1-n2)*np.pi/2)
        f1_cos_taper_manual[:n2] = 0.0
        self.assertTrue(np.allclose(f1_cos_taper, f1_cos_taper_manual))
        
        
            

        
        

if __name__ == '__main__':
    unittest.main()