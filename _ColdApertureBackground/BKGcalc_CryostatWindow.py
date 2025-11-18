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
import glob

from astropy.modeling.models import BlackBody
from astropy import units as u

from utils import tophat,load_blocking_filter

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



def sky_flux(band):
    """
    http://astroweb.case.edu/ssm/ASTR620/mags.html

    outputs:
    --------
    photon flux
    """
    if band=="J":
        dl_l = 0.16
        zp = 1600
        mag = 16.1
    #mag = 16.1 # per asec**2 https://www.gemini.edu/observing/telescopes-and-sites/sites
    if band=='H':
        dl_l = 0.23
        zp = 1080
        mag = 13.8
    if band=='K':
        dl_l = 0.23
        zp = 670
        mag = 14.8

    phot_per_s_m2_per_Jy = 1.51*10**7 # convert to phot/s/m2 from Jansky

    phot_per_s_m2 = dl_l * zp * 10**(-0.4*mag) * phot_per_s_m2_per_Jy
    area      = 76 #m2
    thrput    = 0.25 # to tracking camera
    plt_scale = 0.0106 # asec/pixel

    sky_bkg = phot_per_s_m2 * area * thrput * plt_scale**2

    return sky_bkg



def compare_rebecca_asahi_filters(wave):
    """
    Compare rebecca's blocking filter data vs the asahi data

    rebecca's look like it hit some measurement floor
    """
    fx,fy   = np.loadtxt('./BlockingFilter/rebeccas_blocking_filter.TXT',skiprows=20).T
    fx2,fy2 = np.loadtxt('./BlockingFilter/asahi_filter_7deg.TXT').T

    plt.figure()
    plt.semilogy(fx,fy,label='Rebeccas')
    plt.plot(fx2,fy2,label='Asahi')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Throughput [%]')
    plt.title('Blocking Filter Comp')
    plt.legend()




def load_window_emissivity(wave):
    """
    load emissivity of 5mm window which is close to thickness
    we may have
    https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=15605
    I already converted to emissivity in T2_Data.xlsx where I copied the data from
    This does not include reflectance

    old infrasil file may include reflectance so stopped using it

    interpolates onto wave (nm) and returns emissivity of 5mm infrasil window
    """
    fx,fy = np.loadtxt('./emissivity/T2_Data_WGN105T2_infrasil_emissivity.txt').T
    f = interpolate.interp1d(fx*u.nm,fy,bounds_error=False,fill_value=0)

    return f(wave)

def compute_window_bb_radiation():
    """
    compute number of photons per wavelength being
    emitted from the window itself not considering
    what is seen by a pixel

    this is so Mitsuko can scale the zemax NSQ model
    so the window brightness is scaled correctly 
    to the star brightness

    Notes
    - do not include blocking filter

    """
    D_lens = 25.4*u.mm # mm lens diameter
    f      = 96*u.mm # mm focal length
    QE     = 0.9
    p      = 18*u.um # pixel pitch
    area   = np.pi * (0.5*D_lens)**2 # 
    # pixel_omega = 1.13**2 * np.pi * p**2 / f**2 * u.radian**2

    wave  = u.nm * np.arange(700,2500,0.01) # constrain to relevant lambda range
    window_emissivity = load_window_emissivity(wave)
    blocking_filter   = load_blocking_filter(wave)
    #track_filter     = tophat(wave.value,1170, 1327,1,0) # J 1170, 1327, K 1990, 2460
    #track_filter      = tophat(wave.value,1990, 2460,1,0) # J 1170, 1327, K 1990, 2460 H 1490, 1780]
    track_filter      = tophat(wave.value,1490, 1780,1,0)

    bbtemp_fxn  = BlackBody(277*u.K, scale=1.0 * u.erg / (u.micron * u.s * u.cm**2 * u.arcsec**2)) 
    solid_angle = 4 * np.pi * u.radian**2
    bb   = area * solid_angle * bbtemp_fxn(wave)

    bb_spec_dens = track_filter * QE * window_emissivity * bb.to(u.photon/u.s/u.nm, equivalencies=u.spectral_density(wave))
    bb_spec_dens_block =  blocking_filter * bb_spec_dens

    # integrate over bandpass
    thermal = np.trapz(bb_spec_dens,wave)
    thermal_block = np.trapz(bb_spec_dens_block,wave)

    # plot
    plt.figure()
    plt.semilogy(wave, bb_spec_dens,label='%s'%np.round(thermal))
    #plt.plot(wave, bb_spec_dens_block,label='%s'%np.round(thermal_block))
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Thermal Flux [ph/nm/s]')
    plt.subplots_adjust(left=0.15,bottom=0.15)
    plt.grid()
    plt.title('Window Radial Emission')
    plt.legend()


def cold_aperture_calc(f_num=8,ploton=False):
    """
    """
    # load and interp blocking filter
    wave  = u.nm * np.arange(400,2600,0.01)
    
    blocking_filter   = load_blocking_filter(wave)
    window_emissivity = load_window_emissivity(wave)
    QE = tophat(wave.value,600,2600,0.9,0) # p of h2rg
    
    temp = 277 * u.K
    pixel_size = 18 * u.micron
    #f_num = 8
    #texp = 0.1  * u.second

    area_times_omega = u.radian**2 * 1.13**2 * np.pi**2 * pixel_size**2 / 4 /f_num**2
    bbtemp_fxn  = BlackBody(temp, scale=1.0 * u.erg / (u.micron * u.s * u.cm**2 * u.arcsec**2)) 
    bb   = area_times_omega.to(u.cm**2 * u.arcsec**2) * bbtemp_fxn(wave)

    bb_spec_dens = bb.to(u.photon/u.s/u.nm, equivalencies=u.spectral_density(wave))
    thermal_spectrum = QE * window_emissivity * blocking_filter * bb_spec_dens
    thermal = np.trapz(thermal_spectrum,wave)

    if ploton:
        fig, ax = plt.subplots(3,1,figsize=(4,8),sharex=True,sharey=False)
        plt.subplots_adjust(left=0.25,hspace=0)
        ax[0].plot(wave, bb_spec_dens,label='F# = %s'%f_num)
        ax[1].plot(wave, QE,label='QE')
        ax[1].plot(wave,blocking_filter,label='Blocking\nFilter')
        ax[2].plot(wave,window_emissivity)

        ax[0].set_ylabel('Emission [ph/nm/s]')
        ax[1].set_ylabel('Transmission')
        ax[2].set_ylabel('Window Emissivity')

        ax[2].set_xlabel('Wavelength [nm]')

        ax[0].legend(fontsize=10,loc='best')
        ax[1].legend(fontsize=8,loc='best')

        ax[0].grid(True)
        ax[1].grid(True)
        ax[2].grid(True)

        ax[0].set_title('Background: %s '%thermal.round(2))
        plt.savefig('background_components.png')

    return bb_spec_dens, thermal_spectrum, thermal



if __name__=='__main__':
    f_nums = np.arange(4,35)
    bkgs   = []
    for f_num in f_nums:
        _,_,thermal = cold_aperture_calc(f_num=f_num,ploton=False)
        bkgs.append(thermal.value)

    plt.figure()
    plt.plot(f_nums,bkgs)
    plt.grid(True)
    plt.title('Cold Aperture F/# Trade Off')
    plt.xlabel('F Number')
    plt.ylabel('Thermal Background [ph/s/pix]')
    plt.savefig('bkgs_vs_Fnum.png')



