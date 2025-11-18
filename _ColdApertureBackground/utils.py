# calc thermal background on HISPEC ATC
# just considering photons from a warm cryostat window 
# without a set f/#
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import pandas as pd
from datetime import date

from astropy.modeling.models import BlackBody
from astropy import units as u

font = {'size'   : 14}
matplotlib.rc('font', **font)

plt.ion()

from matplotlib.ticker import (AutoMinorLocator)

plt.ion()
font = {'size'   : 16}
matplotlib.rc('font', **font)
plt.rcParams['font.size'] = '14'
plt.rcParams['font.family'] = 'sans'
plt.rcParams['axes.linewidth'] = '1.3'
fontname = 'DIN Condensed'


def tophat(x,l0,lf,throughput_in,throughput_out):
    ion = np.where((x > l0) & (x<lf))[0]
    bandpass = np.zeros_like(x) + throughput_out
    bandpass[ion] = throughput_in
    return bandpass

def load_blocking_filter(wave):
    """
    load rebeccas blocking filter profile
    interpolates onto input wave (nm)

    inputs
    ------
    wave - wavelength array

    outputs
    --------
    blocking filter array interpolated onto wave in fractional units
    """
    #fx,fy = np.loadtxt('./BlockingFilter/rebeccas_blocking_filter.TXT',skiprows=20).T
    #fx,fy = np.loadtxt('./BlockingFilter/asahi_filter_7deg.TXT').T
    fx,fy = np.loadtxt('./BlockingFilter/asahi_filter_f6_0AOI.TXT').T
    f = interpolate.interp1d(fx[::-1]*u.nm,fy[::-1],bounds_error=False,fill_value=0)

    return f(wave)/100 # convert from percent



def load_window_emissivity(wave):
    """
    load emissivity of 5mm window which is close to thickness
    we may have
    https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=15605
    I already converted to emissivity in T2_Data.xlsx where I copied the data from
    This does not include reflectance

    old infrasil file may include reflectance so stopped using it

    interpolates onto wave (nm) and returns emissivity of 5mm infrasil window
    
    inputs
    ------
    wave - wavelength array

    outputs
    --------
    emissivity
    """
    fx,fy = np.loadtxt('./emissivity/T2_Data_WGN105T2_infrasil_emissivity.txt').T
    f = interpolate.interp1d(fx*u.nm,fy,bounds_error=False,fill_value=0)

    return f(wave)







