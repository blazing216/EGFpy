# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 09:55:51 2018

@author: xuyihe
"""
import warnings
warnings.filterwarnings('always')

import numpy as np
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth
from scipy.signal import hilbert
from scipy.interpolate import interp1d, pchip
import obspy

# TODO(xuyihe, 20180418): change 'readCF' to return t, data, 
# lon1, lat1, lon2, lat2 and provide a 'loc2dist' function to compute
# distance
def readCF(CFfile, return_staloc=True):
    data = np.loadtxt(CFfile, usecols=range(3),
                      skiprows=2)
    t = data[:,0]
    #WinWave = data[:,1]
    WinWave = 0.5 * (data[:,1] + data[:,2])
    if return_staloc is True:
        with open(CFfile, 'rU') as cf:
            cf.seek(0)
            lines = cf.readlines()[:2]
            lon1, lat1 = lines[0].strip().split()[:2]
            lon2, lat2 = lines[1].strip().split()[:2]
        lon1, lat1, lon2, lat2 = \
            float(lon1), float(lat1), float(lon2), float(lat2)
        dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)[0] / 1000.0
        return t, WinWave, dist
    else:
        return t, WinWave

def readSAC(sacfile, return_staloc=True):
    st = obspy.read(sacfile)
    tr = st[0]
    sachd = tr.stats.sac
    data = tr.data
    t = tr.times() + sachd.b
    assert(sachd.npts % 2 == 1)
    pos_branch = slice((sachd.npts-1)/2, sachd.npts)
    neg_branch = slice(0, (sachd.npts+1)/2)
    t = t[pos_branch]
    WinWave = 0.5 * (data[pos_branch] + data[neg_branch][::-1])
    
    if return_staloc is True:
        return t, WinWave, sachd.dist
    else:
        return t, WinWave
    
    
def get_alpha(dist):
    """get alpha for guassian narrow-band filter 
    from an empirical relation of alpha versus distance.
    Longer the distance is, larger the alpha is.
    Distance is in km
    """
    dist_ref  = [0, 100, 250, 500, 1000, 2000, 4000, 20000]
    alpha_ref = [5, 8,   12,  20,  25,   35,   50,   75]
    return np.interp(dist, dist_ref, alpha_ref)
    
def nextpow2(a):
    """return an integer n so that 2**(n-1) < a <= 2**n
    """
    return int(np.ceil(np.log2(a)))
    
def cos_taper(data, n1, n2):
    """taper the data to reduce it gradually to zero.
    taper from n1 (1.0) to n2 (0.0). 
    n1 could be larger than n2 to make a left taper
    
    data_taper[n1] == data[n1]
    data_taper[n2] == 0
                                                    k     pi
    data_taper[n1+1:n2] = data[n1+1:n2] * np.cos(-------* -- ), k=1,..,n2-n1-1
                                                  n2-n1   2
    or                                      
                                                    k     pi
    data_taper[n1:n2+1] = data[n1:n2+1] * np.cos(-------* -- ), k=0,..,n2-n1
                                                  n2-n1   2
    """
    if isinstance(data, list):
        data = np.array(data)
    data_taper = data.copy()
    N = len(data_taper)
    
    assert(isinstance(n1, int) and isinstance(n2,int) and \
        n1 >=0 and n1 < N \
        and n2 >=0 and n2 < N \
        and n1 != n2)
    
    if n2 > n1:
        k = np.arange(n2-n1+1, dtype=float)
        data_taper[n1:n2+1] *= np.cos(k / (n2-n1) * np.pi/2)
        data_taper[n2+1:] = 0.0
    else: # n2 < n1 because we exclude n1==n2 case in the assertation above
        k = np.arange(n1-n2+1, dtype=float)
        data_taper[n2:n1+1] *= np.cos(k / (n1-n2) * np.pi/2 )[::-1]
        data_taper[:n2] = 0.0
    
    return data_taper

def taper_gv(data, delta, dist, b=0.0,
             gvmin=2, gvmax=5, #gvdelta=0.1,
             taper_length=20.0):
    """taper waveform based on the target group velocity range.
    signals arrived before dist/gvmax and after dist/gvmin are reduced.
         dist/gvmax-taper_length   dist/gvmax   dist/gvmin   dist/gvmin+taper_length
                  |                      |        |                      |
    beyond        v                      v        v                      v       beyond
     0 ...        0    <-- cos_taper     1,  ..., 1  -->  cos_taper  --> 0       0
     """
    assert(gvmax > gvmin)
    fs = 1.0 / delta
    npts = len(data)
    
    n1 = int(np.round(fs * (dist/gvmax - taper_length)))
    n2 = int(np.round(fs * dist/gvmax))
    n3 = int(np.round(fs * dist/gvmin))
    n4 = int(np.round(fs * (dist/gvmin + taper_length)))
    
    if n3 >= npts-1:
        n3 = npts-2
        n4 = npts-1
        e = b + (npts-1) * delta
        gvmin = np.ceil(10.0*dist/e) * 10.0
    assert(n2 < npts-1)
    n1 = max(0,n1)
    #n3 = min(npts-2, n3)
    n4 = min(npts-1, n4)
    
    data = cos_taper(data, n2, n1)
    data = cos_taper(data, n3, n4)
    #print n1, n2, n3, n4
    
    return data, gvmin, gvmax 
    
def gv_image(data, delta, dist, b=0.0, 
             gvmin=2.0, gvmax=5.0, gvdelta=0.1,
             Tmin=5.0, Tmax=50.0, Tdelta=1.0):
    # T for period, gv for group velocity
    # The two array determine the resolution of the final group velocity image
    T = np.arange(Tmin, Tmax+0.5*Tdelta, Tdelta, dtype=float)
    nT = len(T)
    gv = np.arange(gvmin, gvmax+0.5*gvdelta, gvdelta, dtype=float)
    ngv = len(gv)

    # Transform waveform to frequency domain using FFT
    # 
    # ensure the sampling rate in frequency domain >= 1024
    fs = 1.0 / delta
    npts = len(data)    
    nfft = 2**nextpow2(max(npts, 1024*fs))
    D = np.fft.fft(data, nfft) # Spectrum of data
    # sampling points in frequency domain
    f = np.arange(0,nfft/2+1, dtype=float)/nfft * fs  
    
    # parameters of Gaussian narrow band filter
    alpha = get_alpha(dist)
    
    # for storing envelope after each Gaussian filter
    env_img = np.zeros((nT, npts)) 
    
    # group velocity points corresponds to time domain sampling points
    # used for interpolate envelope to group velocity image
    t = b + np.arange(0, npts) * delta
    pgv = dist / t[::-1] 
    
    # for storing group velocity image
    gv_img = np.zeros((ngv, nT)) 
    
    for i, Ti in enumerate(T):
        # construct Gaussian filter
        fc = 1/Ti;                
        Hf = np.exp(-alpha*(f - fc)**2/fc**2);
        
        # compute analytic signal of data
        # data_anal = data + j * hilbert(data)
        Da = np.zeros(nfft, dtype=complex); # for analytic spectrum of data
        Da[0:(nfft/2+1)] = D[0:(nfft/2+1)] * Hf
        assert(nfft%2 == 0)
        Da[1:(nfft/2)] *= 2.0
        data_anal = np.fft.ifft(Da, nfft)
        
        # compute envelope of data
        # data_env = |data_anal|
        data_env = np.abs(data_anal)
        env_img[i,:] = data_env[0:npts]
        
        # construct group velocity image
        # normalized for each frequency band
        maxamp = np.max(env_img[i,:])
        #f_interp = pchip(pgv[:-1], env_img[i,:0:-1] / maxamp)
        f_interp = interp1d(pgv[:-1], env_img[i,:0:-1] / maxamp, 'cubic',
                            bounds_error=False)
        gv_img[:,i] = f_interp(gv)
    return gv_img

# TODO(xuyihe, 20180420): modified code to make return_staloc work properly
def EGF_from_file(cfile, return_staloc=True, cfile_format='EGF'):
    if cfile_format.upper() == 'EGF':
        t, EGF, dist = readCF(cfile, return_staloc)
    elif cfile_format.upper() == 'SAC_CF':
        #print type(cfile)
        t, CF, dist = readSAC(cfile.decode('gbk'), return_staloc)
        EGF = hilbert(CF).imag
    
    if return_staloc:
        return t, EGF, dist
    else:
        return t, EGF

def win_EGF_from_file(cfile, gvmin=2, gvmax=5, taper_length=20.0, cfile_format='EGF'):
    t, EGF, dist = EGF_from_file(cfile, cfile_format=cfile_format)
    delta = t[1] - t[0]
    
    win_EGF, gvmin_m, gvmax_m = \
        taper_gv(EGF, delta, dist, b=0.0,
                 gvmin=gvmin, gvmax=gvmax, #gvdelta=gvdelta,
                 taper_length=taper_length)
    return t, win_EGF, dist#, gvmin_m, gvmax_m

def gv_image_from_EGF(t, EGF, dist, gvmin=2, gvmax=5, gvdelta=0.1, 
                       Tmin=5, Tmax=50, Tdelta=1,
                       taper_length=20.0):
    EGF /= np.abs(EGF).max()
    delta = t[1] - t[0]
    
    win_EGF, gvmin_m, gvmax_m = \
        taper_gv(EGF, delta, dist, b=0.0,
                 gvmin=gvmin, gvmax=gvmax, #gvdelta=gvdelta,
                 taper_length=taper_length)
    #n4 = int(np.round((dist/gvmin + taper_length)/delta))
    #if n4 < len(win_EGF)-1:      
        #print n4
        #win_EGF_trunc = np.zeros(n4+1)
        #win_EGF_trunc = win_EGF[:(n4+1)]
    gv_img = gv_image(win_EGF, delta, dist, b=0.0,
                      gvmin=gvmin_m, gvmax=gvmax_m, gvdelta=gvdelta,
                      Tmin=Tmin, Tmax=Tmax, Tdelta=Tdelta)
    
    T = np.arange(Tmin, Tmax+0.5*Tdelta, Tdelta)
    gv = np.arange(gvmin_m, gvmax_m+0.5*gvdelta, gvdelta)   
    #print gvmin_m, gvmax_m           
    #return GrpVelImg(gv_img, T, gv)
    return gv_img, T, gv                       

def gv_image_from_EGF_file(cfile, cfile_format='EGF', **kwargs):
    t, EGF, dist = EGF_from_file(cfile, cfile_format=cfile_format)    
    gv_img, T, gv = gv_image_from_EGF(t, EGF, dist, **kwargs)
    return gv_img, T, gv

class GrpVelImg(object):
    def __init__(self, img_or_EGFfile, cfile_format, **kwargs):
        if isinstance(img_or_EGFfile, str):
            self.t, self.EGF, self.dist = EGF_from_file(img_or_EGFfile,
                                                        cfile_format=cfile_format)
            img, T, gv = gv_image_from_EGF(self.t, self.EGF, self.dist, **kwargs)
            #img, T, gv = gv_image_from_EGF_file(img_or_EGFfile, **kwargs)
        else:
            img = img_or_EGFfile
            assert(T is not None and gv is not None)
        self.img = img
        self.T = T
        self.gv = gv
        half_dT = 0.5 * (T[1] - T[0])
        half_dgv = 0.5 * (gv[1] - gv[0])
        self.extent = [T[0]-half_dT, T[-1]+half_dT, \
            gv[0]-half_dgv, gv[-1]+half_dgv]
        self.cmap = 'jet'
        # scaling the y units by 0.6 / (gv_len/T_len), 
        # to show the image in a 3:5 scale
        T_len = T[-1] - T[0]
        gv_len = gv[-1] - gv[0]
        self.aspect = 0.6 * T_len / gv_len
    
    def plot(self, ax=None, extent=None, cmap=None, #with_colorbar=False,
             **kwargs):
        ax = plt.gca() if ax is None else ax
        extent = self.extent if extent is None else extent
        cmap = self.cmap if cmap is None  else cmap
        
        img = ax.imshow(self.img, origin='lower', 
                  extent=extent, cmap=cmap, **kwargs)
        for c in [1.0, 2.0, 3.0]:
            x = self.dist / (2.0*c)
            ax.plot([x,x],[0,1],'y--', transform=ax.get_xaxis_transform())
        ax.set_aspect(self.aspect)
        return img
#        if with_colorbar:
#            plt.colorbar()


    
defaultCFfile = 'EGFAnalysisTimeFreq_version_2015\EGFs\GFcn.KMI-MC01_10-50s_10Mon.dat'

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    GroupVImage_class = GrpVelImg(defaultCFfile)
    
    t, EGF, dist = readCF(defaultCFfile)
    SampleT = t[1]-t[0]
    WinMinV = 2
    WinMaxV = 5
    #
    # compute group velocity image using subfunctions
    #
    WinWave_taper_gv, gvmin_m, gvmax_m = \
        taper_gv(EGF, SampleT, dist, b=0.0,
                 gvmin=2, gvmax=5, #gvdelta=(WinMaxV-WinMinV)/29.0,
                 taper_length=20.0)
    GroupVImage_subfun = \
        gv_image(WinWave_taper_gv, SampleT, dist, b=0.0,
                 gvmin=gvmin_m, gvmax=gvmax_m, gvdelta=(gvmax_m-gvmin_m)/29.0,
                 Tmin=5, Tmax=50, Tdelta=1,
                 )

    PtNum = len(EGF)
    fs = 1.0 / (t[1]-t[0]) # sampling rate
    # make window
    g_WinMaxPtNum = int(np.round(fs*dist/WinMinV))
    g_WinMinPtNum = int(np.round(fs*dist/WinMaxV))
    if g_WinMaxPtNum > PtNum-1:
        g_WinMaxPtNum = PtNum-2
        WinMinV = np.ceil(10*dist/t[-1])/10;
        print 'Min velocity reset to ', WinMinV
    
    
    Window = np.zeros(PtNum)
    Window[g_WinMinPtNum:g_WinMaxPtNum+1] = 1.0
    TaperNum = int(np.round(20/SampleT)) # 20 seconds Taper?
    if TaperNum > g_WinMinPtNum:
        TaperNum1 = g_WinMinPtNum
        Window[:g_WinMinPtNum] = np.sin(0.5*np.pi*np.arange(0,TaperNum1)/TaperNum1);
    else:
        Window[(g_WinMinPtNum-TaperNum):g_WinMinPtNum] = np.sin(0.5*np.pi*np.arange(0,TaperNum)/TaperNum)
        #print g_WinMinPtNum-TaperNum, g_WinMinPtNum
    
    if (g_WinMaxPtNum + TaperNum) < PtNum:
        Window[(g_WinMaxPtNum+1):(g_WinMaxPtNum+TaperNum+1)] = np.sin(0.5*np.pi*np.arange(TaperNum-1,-1,-1)/TaperNum);
        #print g_WinMaxPtNum, g_WinMaxPtNum+TaperNum
    else:
        TaperNum2 = PtNum - (g_WinMaxPtNum+1);
        Window[(g_WinMaxPtNum+1):PtNum] = np.sin(0.5*np.pi*np.arange(TaperNum2-1,-1,-1)/TaperNum2);
    
    WinWave = EGF * Window  

    alpha = get_alpha(dist)
    
    TPoint = np.linspace(5, 50, 46)    
    
    NumCtrT = len(TPoint)

    nfft = 2**nextpow2(max(PtNum,1024*fs))# ensure the frequency sampling rate > 1024
    xxfft = np.fft.fft(WinWave, nfft)
    fxx = np.linspace(0,nfft/2, nfft/2+1)/nfft*fs
    
    
    EnvelopeImage = np.zeros((NumCtrT, PtNum))
    # calculate envelope using alternative approach 
    # multiple 2 on positive frequency components
    EnvelopeImage_alt = np.zeros((NumCtrT, PtNum)) 
    for i in range(NumCtrT):
    
        CtrT = TPoint[i]
        fc = 1/CtrT;                
        Hf = np.exp(-alpha*(fxx - fc)**2/fc**2);
        yyfft = np.zeros(nfft, dtype=complex);
        yyfft[0:(nfft/2+1)] = xxfft[0:(nfft/2+1)]*Hf
        yyfft[(nfft/2+1):nfft] = np.conjugate(yyfft[nfft/2-1:0:-1])
        yy = np.fft.ifft(yyfft, nfft).real
        filtwave = np.abs(hilbert(yy))
        EnvelopeImage[i,:] = filtwave[0:PtNum]
        
        even_nfft = nfft%2 == 0
        assert(even_nfft)
        yyfft_alt = np.zeros(nfft, dtype=complex)
        yyfft_alt[0:(nfft/2+1)] = xxfft[0:(nfft/2+1)] * Hf
        yyfft_alt[1:(nfft/2)] *= 2.0
        yy_alt = np.fft.ifft(yyfft_alt, nfft)
        
        filtwave_alt = np.abs(yy_alt)
        EnvelopeImage_alt[i,:] = filtwave_alt[0:PtNum]
        
        assert(np.allclose(filtwave, filtwave_alt))
    
    AmpS_T = np.max(EnvelopeImage, axis=1)
    AmpS_T_alt = np.max(EnvelopeImage_alt, axis=1)
    VImgPt = np.linspace(WinMinV, WinMaxV, 30)
    GroupVImage = np.zeros((len(VImgPt), NumCtrT))   
    GroupVImage_alt = np.zeros((len(VImgPt), NumCtrT)) 
    TravPtV = dist / t[::-1]
    for i in range(NumCtrT):
        f = interp1d(TravPtV[:-1], EnvelopeImage[i,:0:-1] / AmpS_T[i], 'cubic')
        GroupVImage[:,i] = f(VImgPt)
        f_alt = interp1d(TravPtV[:-1], EnvelopeImage_alt[i,:0:-1] / AmpS_T_alt[i], 'cubic')
        GroupVImage_alt[:,i] = f_alt(VImgPt)
    
    plt.figure(figsize=(4,9))
    plt.subplot(5,1,1)
    plt.imshow(GroupVImage, origin='lower', cmap='jet')
    plt.xlabel('Period (s)')
    plt.ylabel('Velocity (km/s)')
    plt.colorbar(shrink=0.5)
    
    plt.subplot(5,1,2)
    plt.imshow(GroupVImage_subfun, origin='lower', cmap='jet')
    plt.xlabel('Period (s)')
    plt.ylabel('Velocity (km/s)')
    plt.colorbar(shrink=0.5)
    
    plt.subplot(5,1,3)
    plt.imshow(GroupVImage_subfun - GroupVImage, origin='lower', cmap='jet')
    plt.xlabel('Period (s)')
    plt.ylabel('Velocity (km/s)')
    plt.colorbar(shrink=0.5)
    
    plt.subplot(5,1,4)
    #plt.imshow(GroupVImage_class, origin='lower', cmap='jet')
    GroupVImage_class.plot()
    plt.xlabel('Period (s)')
    plt.ylabel('Velocity (km/s)')
    plt.colorbar(shrink=0.5)
    
#    plt.subplot(5,1,5)
#    plt.imshow(GroupVImage_helperfun - GroupVImage_subfun, origin='lower', cmap='jet')
#    plt.xlabel('Period (s)')
#    plt.ylabel('Velocity (km/s)')
#    plt.colorbar(shrink=0.5)
    
    plt.show()
    #print np.abs(GroupVImage_alt - GroupVImage).sum()
    #print np.abs(GroupVImage_subfun - GroupVImage).sum()
