# 
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import pandas as pd
from datetime import date

font = {'size'   : 14}
matplotlib.rc('font', **font)

plt.ion()

from matplotlib.ticker import (AutoMinorLocator)

plt.ion()
font = {'size'   : 14}
matplotlib.rc('font', **font)
plt.rcParams['font.size'] = '14'
plt.rcParams['font.family'] = 'sans'
plt.rcParams['axes.linewidth'] = '1.3'
fontname = 'DIN Condensed'


def load_brown_dwarfs_AB():
    """
    load adam burgasser's file
    """
    xl = pd.ExcelFile('/Users/ashbake/Documents/Research/Projects/HISPEC/Tests/FEI_TrackingCamera/data/ucd_sheet_teff.xlsx')
    #xl.sheet_names
    first_sheet = xl.sheet_names[0]
    df = xl.parse(first_sheet)
    #df.head()
    teff = df['teff']
    hmag = df['H_2MASS']

    return  hmag.values,teff.values

def load_planets():
    planets_filename = '/Users/ashbake/Documents/Research/_DATA/populations/exoplanets/uncontroversial_planets_PS_2023.06.29_11.49.57.csv'
    planet_data =  pd.read_csv(planets_filename,delimiter=',',comment='#')
    # add brown dwarfs!
    hmags = planet_data['sy_hmag']
    teffs = planet_data['st_teff']
    rvamps = planet_data['pl_rvamp']
    method = planet_data['discoverymethod']
    return hmags,teffs,method,rvamps


if __name__=='__main__':
	#########
    """
    plot planet populations with regions of AO mode and coupling performance

    show with and without a pyramid mode - have to edit which plotting

    To Do: 
        - rerun coupling curves once gary extends tip/tilt grid
        - rerun with denser mag grid
        - label ao regions and double check correct orientation is shown
        - add brown dwarfs to plot
        - do again with contours as peak SNR in a 15 min exposure
    must be careful to pick right ao_modes array since not recorded in header
    """
    band='H'
    ao_modes  = ['LGS_100J_130','SH']
    best_ao_mode_arr   = np.load('/Users/ashbake/Documents/Research/Projects/HISPEC/Tests/FEI_TrackingCamera/data/ao_modes/best_ao_mode_Hmag_temps_1000,1500,2300,3000,3600,4200,5800_modes_%s.npy'%ao_modes)
    magarr             = np.load('/Users/ashbake/Documents/Research/Projects/HISPEC/Tests/FEI_TrackingCamera/data/ao_modes/best_ao_mode_Hmag_mags.npy')
    #ao_modes  = ['LGS_100H_130','SH']
    #best_ao_mode_arr= np.load('./output/ao_modes/best_ao_mode_Hmag_temps_1000,1500,2300,3000,3600,4200,5800_modes_[LGS_100J_130_SH.npy')
    #best_ao_coupling_arr  =  np.load('./output/ao_modes/best_ao_coupling_%smag_[1,15,2]_temps_1000,1500,2300,3000,3600,4200,5800.npy'%band)

    #best2_ao_mode_arr      =  np.load('./output/ao_modes/best2_ao_mode_Hmag_[1,15,2]_temps_1000,1500,2300,3000,3600,4200,5800.npy')
    #best2_ao_coupling      =  np.load('./output/ao_modes/best2_ao_coupling_%smag_[1,15,2]_temps_1000,1500,2300,3000,3600,4200,5800.npy'%band)
    #ao_modes2   = ['SH','LGS_100H_130']

    temp_arr   = [1000,1500,2300,2700,3000,3600,4200,5800]
    extent     = (np.min(temp_arr),np.max(temp_arr),np.min(magarr),np.max(magarr))

    # plot contours of it
    hmags, teffs,method,rvamps = load_planets()
    i_image = np.where(method=='Imaging')[0]
    i_ms= np.where(rvamps<3)[0]
    #hmags_bd, teffs_bd = load_brown_dwarfs()
    hmags_bd, teffs_bd = load_brown_dwarfs_AB()

    def fmt(x):
        "format labels for contour plot"
        return str(np.round(100*x)) + '%'

    ##########
    fig, ax = plt.subplots(1,1,figsize=(6,5))
    ax.imshow(best_ao_mode_arr.T,aspect='auto',origin='lower',\
            interpolation='None',cmap='gray',\
            extent=extent,vmin=0,vmax=2)

    ax.scatter(teffs,hmags,marker='.',s=65,c='m',ec='purple',alpha=0.8,label='Confirmed Planet Hosts')
    #ax.scatter(teffs[i_image],hmags[i_image],marker='s',s=40,c='brown',ec='k',alpha=1,label='Imaged Planets')
    ax.scatter(teffs_bd,hmags_bd,marker='d',s=32,c='darkcyan',ec='b',alpha=.8,label='Brown Dwarfs')
    #ax.scatter(hmags_bd,teffs_bd,marker='.',c='c',alpha=0.5,label='Brown Dwarfs')
    ax.set_ylabel('H Mag')
    ax.set_xlabel('T$_{eff}$ (K)')
    ax.set_title('HISPEC AO Mode Landscape')
    ax.set_ylim(0,16)
    ax.set_xlim(1000,5800)
    #ax.legend(fontsize=10,loc=1)
    plt.subplots_adjust(bottom=0.15,hspace=0.1,left=0.15,right=0.85)

    # add region for cred2 interrupted science region
    teff_limit_jhgap = [1000,1500,2300,3000,4200,6000]
    mode = 'h2rg'
    if mode=='cred2':
        mag_limit_h_cred2  =     [14.5,13.9,13.9,13.9,15,15]
        mag_limit_jhgap_cred2  = [12.3,12.6,13,13.6,14.2,15]
        mag_limit_jhplus_cred2 = [15.2,14.7,14.7,14.8,14.9,15]
        mag_limit_j_cred2      = [14.6,13.7,13.6,13.6,14.4,15]
        mag_limit_jhgap_TT_cred2   = [12.7,12.9,13.8,14.1,14.9,15]
        ax.plot(teff_limit_jhgap,mag_limit_h_cred2,'-',c='orange',lw=2,label='H')
        ax.plot(teff_limit_jhgap,mag_limit_j_cred2,'-',c='pink',lw=2,label='J')
        ax.plot(teff_limit_jhgap,mag_limit_jhplus_cred2,'-',c='limegreen',lw=2,label='J+H')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_cred2,'-',c='yellow',lw=2,label='JH Gap')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_TT_cred2,'--',c='yellow',lw=2,label='JH Gap + TT Star')
        #ax.set_title('Uninterrupted Science Limits')
        ax.set_title('CRED2 Tracking')
    if mode=='cred2_2': # with different read noise values
        mag_limit_yjh_rn30  =  [15.3,14.7, 14.7,14.9, 15.1, 16] # cred2
        mag_limit_yjh_rn25  =  [15.5,14.8,14.8,14.95,15.2,15.8]
        mag_limit_yjh_rn20  =  [15.6,14.9,14.9,15.1,15.3,15.9]
        mag_limit_jhgap_rn30  = [12.3, 12.55, 13.05, 13.6, 14.2,14.9]
        mag_limit_jhgap_rn25  = [12.6,12.7,13.3,13.6,14.2,14.9]
        mag_limit_jhgap_rn20  = [12.7,12.8,13.55,13.8,14.4,15]                
        ax.plot(teff_limit_jhgap,mag_limit_yjh_rn30,'-',c='orange',lw=2,label='yJH RN30')
        ax.plot(teff_limit_jhgap,mag_limit_yjh_rn25,'--',c='orange',lw=2,label='yjH RN25')
        ax.plot(teff_limit_jhgap,mag_limit_yjh_rn20,'-.',c='orange',lw=2,label='yJH RN20')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_rn30,'-',c='red',lw=2,label='JHgap RN30')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_rn25,'--',c='red',lw=2,label='JHgap RN25')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_rn20,'-.',c='red',lw=2,label='JHgap RN20')        
        #ax.set_title('Uninterrupted Science Limits')
        ax.set_title('CRED2 - Various RNs')
    if mode=='h2rg':
        mag_limit_jhgap_h2rg   = [12.5,12.7,13.2,13.7,14.2,15]
        mag_limit_j_h2rg       = [15.5,14.7,14.7,14.6,16,16]
        #mag_limit_jhgap_TTstar_h2rg   = [13.5,13.8,13.6,14.9,14.9,16]
        mag_limit_jplus_h2rg   = [15.5,14.8,14.9,14.9,16,16]
        mag_limit_jhplus_h2rg  = [15.8,15.6, 15.6, 15.6,16,16]
        mag_limit_h_h2rg       = [15.5,15.6,15.6,15.6,16,16]
        mag_limit_k_h2rg       = [15.5,15.5,15.6,15.6,16,16]
        mag_limit_jhgap_TT_h2rg    = [13.6,13.8,14.6,14.8,15,16]
        ax.plot(teff_limit_jhgap,mag_limit_h_h2rg,'-',c='limegreen',lw=2,label='H')
        #ax.plot(teff_limit_jhgap,mag_limit_jplus_h2rg,'--',c='pink',lw=2,label='Jplus')
        ax.plot(teff_limit_jhgap,mag_limit_j_h2rg,'-',c='pink',lw=2,label='J')
        #ax.plot(teff_limit_jhgap,mag_limit_jhplus_h2rg,'-',c='limegreen',lw=2,label='J+H')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_h2rg,'-',c='yellow',lw=2,label='JHgap')
        #ax.plot(teff_limit_jhgap,mag_limit_jhgap_TT_h2rg,'--',c='yellow',lw=2,label='JHgap + Bright TT Star')
        #ax.plot(teff_limit_jhgap,mag_limit_k_h2rg,'--',c='purple',lw=2,label='K')
        ax.set_title('H2RG Tracking')
    if mode=='TTstar':
        mag_limit_jhgap_perfect    = [13.7,13.8,14.6,14.8,15,16]
        mag_limit_jhgap_TT_h2rg    = [13.6,13.8,14.6,14.8,15,16]
        mag_limit_jhgap_h2rg       = [13.1,13.3,13.8,14.3,15,16]
        mag_limit_jhgap_TT_cred2   = [12.7,12.9,13.8,14.1,14.9,15]
        mag_limit_jhgap_cred2      = [12.3,12.6,13,13.6,14.2,15]
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_TT_h2rg,'--',c='c',lw=2,label='H2RG with TT star')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_TT_cred2,'--',c='red',lw=2,label='C-RED2 with TT star')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_h2rg,'-',c='c',lw=2,label='H2RG')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_cred2,'-',c='red',lw=2,label='C-RED2')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_perfect,'-.',c='yellow',lw=2,label='Perfect Detector')
        ax.set_title('Interrupted Science')
    if mode=='Hband': # updated
    	#mag_limit_jplus_cred2     = [14.6,13.9,13.6,13.6,14.4,15]
		#mag_limit_jplus_h2rg      = [15.5,14.8,14.9,14.9,16,16]
        mag_limit_h_cred2_TT      = [15.1,14.8,14.8,14.8,15,16]
        mag_limit_yjh_cred2_TT    = [16,15.6,15.6,15.8,16,16]
        mag_limit_h_h2rg       = [15.5,15.6,15.6,15.6,16,16]
        mag_limit_h_cred2      = [14.5,13.9,13.9,13.9,15,16]
        mag_limit_jhplus_cred2 = [15.2,14.7,14.7,14.8,14.9,16]
        ax.plot(teff_limit_jhgap,mag_limit_jhplus_cred2,'-',c='pink',lw=2,label='J+H CRED2')
        ax.plot(teff_limit_jhgap,mag_limit_yjh_cred2_TT,'--',c='pink',lw=2,label='y+J+H CRED2 + TT')
        ax.plot(teff_limit_jhgap,mag_limit_h_cred2_TT,'--',c='red',lw=2,label='H CRED2 + TT')
        ax.plot(teff_limit_jhgap,mag_limit_h_h2rg,'-',c='c',lw=2,label='H H2RG')
        ax.plot(teff_limit_jhgap,mag_limit_h_cred2,'-',c='red',lw=2,label='H CRED2')
        ax.set_title('H Band Tracking')
    if mode=='uninterrupted': # updated
        mag_limit_jhgap_perfect   = [13.7,13.8,14.6,14.8,15,16]        
        mag_limit_jplus_perfect   = [15.8,15.5,15.4,15.55,16,16]
        mag_limit_h_perfect       = [15.8,15.3,15.6,15.6,16,16]
        ax.plot(teff_limit_jhgap,mag_limit_h_perfect,'-',c='orange',lw=2,label='H')
        ax.plot(teff_limit_jhgap,mag_limit_jplus_perfect,'-',c='limegreen',lw=2,label='Jplus')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_perfect,'--',c='yellow',lw=1.2,label='JHgap')
        ax.set_title('H2RG Tracking - TT Star')
    if mode=='KAO': # old keck! 1000, 1500 2300 3000, for h2rg, didnt check out cred2
        mag_limit_J_KAO        = [14.05,13.40,13.30,13.4,15,15]        
        mag_limit_H_KAO        = [14.65,14.60,14.70,14.6,15,15]
        mag_limit_JHgap_KAO    = [12.20,12.30,12.80,13.1,15,15]
        mag_limit_jhgap_h2rg   = [12.5,12.7,13.2,13.7,16,16]
        mag_limit_j_h2rg       = [15.5,14.7,14.7,14.6,16,16]
        mag_limit_h_h2rg       = [15.5,15.6,15.6,15.6,16,16]
        ax.plot(teff_limit_jhgap,mag_limit_J_KAO,'-',c='orange',lw=2,label='J')
        ax.plot(teff_limit_jhgap,mag_limit_H_KAO,'-',c='r',lw=2,label='H')
        ax.plot(teff_limit_jhgap,mag_limit_JHgap_KAO,'-',c='c',lw=2,label='JHgap')
        ax.plot(teff_limit_jhgap,mag_limit_j_h2rg,'--',c='orange',lw=1.5,label='HAKA')
        ax.plot(teff_limit_jhgap,mag_limit_h_h2rg,'--',c='r',lw=1.5)
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_h2rg,'--',c='c',lw=1.5)
        ax.set_title('H2RG Tracking - No HAKA')
    ax.plot([1000,6000],[15,15],'-',c='gray',lw=1)
    #ax.plot(teff_limit_jhgap,mag_limit_j_h2rg,'.-',c='cyan',lw=1)
    #ax.plot(teff_limit_jhgap,mag_limit_j_cred2,'.-',c='red',lw=1)

    ax.set_ylim(10,16)
    ax.legend(fontsize=10,loc=4)

    ntotal_bd = len(np.where((hmags_bd < 15) & (teffs_bd >1000) & (teffs_bd < 5800))[0])
    ntotal_pl = len(np.where((hmags < 15) & (teffs >1000) & (teffs < 5800))[0])
    nsub_bd_h2rg = len(np.where((hmags_bd < 13.5) & (teffs_bd >1000) & (teffs_bd < 5800))[0])
    nsub_bd_cred = len(np.where((hmags_bd < 12.5) & (teffs_bd >1000) & (teffs_bd < 5800))[0])
    nsub_bd_nocando = len(np.where((hmags_bd < 15) &(hmags_bd > 14.75) & (teffs_bd >1000) & (teffs_bd < 5800))[0])

    plt.savefig('./output/ao_mode_landscape_brown_dwarfs_trackingcutoff_%s.png'%mode,dpi=300)

