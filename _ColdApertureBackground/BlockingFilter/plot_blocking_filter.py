# calc thermal background on HISPEC ATC
# just considering photons from a warm cryostat window 
# without a set f/#
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt

import pandas as pd
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


def plot_asahi_filters():
    """
    Compare rebecca's blocking filter data vs the asahi data

    rebecca's look like it hit some measurement floor
    """
    fx0,fy0 = pd.read_csv('./asahi_filter_f4_0AOI.TXT',delimiter='\t').values.T
    fx45,fy45 =pd.read_csv('./asahi_filter_f4_45AOI.TXT',delimiter='\t').values.T
    #fx10,fy10 =pd.read_csv('./asahi_filter_f4_10AOI.TXT',delimiter='\t').values.T
    fx25,fy25 =pd.read_csv('./asahi_filter_f4_25AOI.TXT',delimiter='\t').values.T

    plt.figure()
    plt.semilogy(fx0,fy0,label='0 deg AOI')
    plt.semilogy(fx25,fy25,label='25 deg AOI')
    plt.plot(fx45,fy45,label='45 deg AOI')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Throughput [%]')
    plt.title('Blocking Filter Throughput')
    plt.legend(loc=1,fontsize=12)
    plt.grid()
    plt.axvline(1780,c='gray',ls='--',lw=0.8)
    plt.axvline(980,c='gray',ls='--',lw=0.8)
    plt.subplots_adjust(left=0.15,bottom=0.15)
    plt.savefig('Asahi_filter_angles.png')



if __name__=='__main__':
    plot_asahi_filters()

