# calc background for HISPEC on tracking camera
# with mounts and optics in line before ATC window
# assume set F/#s based on model
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import pandas as pd
from datetime import date
import astropy.units as u
from astropy.modeling.models import BlackBody

from utils import load_window_emissivity, load_blocking_filter, tophat

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




def lens_thermal_calc(ploton=False):
    """
    """
    # load and interp blocking filter
    wave  = u.nm * np.arange(400,2600,0.01)
    
    blocking_filter   = load_blocking_filter(wave)
    QE = tophat(wave.value,600,2600,0.9,0) # sensitivity of h2rg
    #mitsuko_blocking_filter = tophat(wave.value,600,1780,1,0.05**2) # matches what mitsuko assumed in models
    
    temp = 293 * u.K # temperature assumed for all optics
    pixel_size = 18 * u.micron
    names  = ['window', 'lens1', 'lens2'] # took out filter and fold mirror - wait for mitsuko zemax analysis
    f_nums = [6,4,9] # off axis: [] on axis: [6,4,9]
    l0,lf  = 1400,1650 # JH band?
    thicknesses = [6,10,10]

    #blackbody spectrum
    bbtemp_fxn  = BlackBody(temp, scale=1.0 * u.erg / (u.micron * u.s * u.cm**2 * u.arcsec**2)) 
    silica_emissivity = 0.05#load_window_emissivity(wave)

    # total thermal per optic
    spectra = {}
    all_thermal = []
    for i,name in enumerate(names):        
        area_times_omega = u.radian**2 * 1.13**2 * np.pi**2 * pixel_size**2 / 4 /f_nums[i]**2
        bb   = area_times_omega.to(u.cm**2 * u.arcsec**2) * bbtemp_fxn(wave)
        bb_spec_dens = bb.to(u.photon/u.s/u.nm, equivalencies=u.spectral_density(wave))

        thermal_spectrum = QE * silica_emissivity * blocking_filter * bb_spec_dens 
        thermal = np.trapz(thermal_spectrum,wave)
        spectra[name] = thermal_spectrum
        all_thermal.append(thermal)

    if ploton:
        fig, ax = plt.subplots(3,1,figsize=(4,8),sharex=True,sharey=False)
        plt.subplots_adjust(left=0.25,hspace=0)
        ax[0].plot(wave, bb_spec_dens,label='F# = %s'%6)
        ax[1].plot(wave, QE,label='QE')
        ax[1].plot(wave,blocking_filter,label='Blocking\nFilter')
        ax[2].plot(wave,silica_emissivity)

        ax[0].set_ylabel('Emission [ph/nm/s]')
        ax[1].set_ylabel('Transmission')
        ax[2].set_ylabel('Window Emissivity')

        ax[2].set_xlabel('Wavelength [nm]')

        ax[0].legend(fontsize=10,loc='best')
        ax[1].legend(fontsize=8,loc='best')

        ax[0].grid(True)
        ax[1].grid(True)
        ax[2].grid(True)

        ax[0].set_title('Background: %s '%thermal.round())
        plt.savefig('background_components.png')

    return bb_spec_dens, spectra, all_thermal


def mount_thermal_calc(ploton=False):
    """
    computes thermal background for assumptions 
    relevant to the optics mounts
    """
    # load and interp blocking filter
    wave  = u.nm * np.arange(400,2700,0.01)
    
    blocking_filter   = load_blocking_filter(wave)
    QE = tophat(wave.value,600,2600,0.9,0) # sensitivity of h2rg
    #mitsuko_blocking_filter = tophat(wave.value,600,1780,1,0.05**2) # matches what mitsuko assumed in models
    
    temp = 293 * u.K # temperature assumed for all optics
    pixel_size = 18 * u.micron
    names  = ['lens1', 'lens2'] # lens 2 will pass through lens 1, not accounted for here
    f_nums_lens = [4,9]#[8.1, 9.9] # lens 1 made big enough to not matter
    f_nums_mount  = [3.9,7.1]#[8,9.4]
    metal_emissivity = 0.1 # gold polished mirror 1-reflectivity, best case

    l0,lf  = 1400,1650 # JH band?

    #blackbody spectrum
    bbtemp_fxn  = BlackBody(temp, scale=1.0 * u.erg / (u.micron * u.s * u.cm**2 * u.arcsec**2)) 
    
    # total thermal per optic
    spectra = {}
    all_thermal = []
    for i,name in enumerate(names):
        area_times_omega_mount = u.radian**2 * 1.13**2 * np.pi**2 * pixel_size**2 / 4 /(f_nums_mount[i]**2)
        area_times_omega_lens = u.radian**2 * 1.13**2 * np.pi**2 * pixel_size**2 / 4 /(f_nums_lens[i]**2)
        area_times_omega = area_times_omega_mount - area_times_omega_lens
        bb   = area_times_omega.to(u.cm**2 * u.arcsec**2) * bbtemp_fxn(wave)
        bb_spec_dens = bb.to(u.photon/u.s/u.nm, equivalencies=u.spectral_density(wave))

        thermal_spectrum = QE * metal_emissivity * blocking_filter * bb_spec_dens
        thermal = np.trapz(thermal_spectrum,wave)
        spectra[name] = thermal_spectrum
        all_thermal.append(thermal)

    return bb_spec_dens, spectra, all_thermal


if __name__=='__main__':
    bb_spec_dens, spectra, all_thermal = lens_thermal_calc(ploton=False)
    bb_spec_dens2, spectra2, all_thermal2 = mount_thermal_calc(ploton=False)

    print(all_thermal)
    print(all_thermal2)


