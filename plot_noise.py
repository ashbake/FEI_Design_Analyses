# calc signal to noise
# max for the calcium H&K
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


def load_data(camera,ao_mode,track_band,teff,datapath):
	"""
	"""
	# pick right ao mode
	if '%s' in ao_mode:
		if track_band=='H': ao_mode_filled = ao_mode%track_band
		else: ao_mode_filled = ao_mode%'J'
	else: ao_mode_filled=ao_mode

	savename = 'cendata_%s_ao_%s_band_%s_teff_%s_%smag_fieldr_%sarcsec.npy' %(camera,ao_mode_filled,track_band,teff,'H',0.0)
	try: 
		exptimes = np.load(datapath + 'exptimes.npy')
		magarr   = np.load(datapath + 'magarr.npy')
		centroid = np.load(datapath + 'centroid_%s'%(savename))
		fwhm = np.load(datapath + 'fwhm_%s'%savename)
		signal  = np.load(datapath + 'signal_%s'%savename)
		noise  = np.load(datapath + 'noise_%s'%savename)
		strehl = np.load(datapath + 'strehl_%s'%savename)
		snr = np.load(datapath + 'snr_%s'%savename)
		skyinst_bg = np.load(datapath + 'skyinstbg_%s'%savename)
	except: return

	return exptimes,magarr,centroid,fwhm, signal, noise, strehl, centroid, snr, skyinst_bg


def plot_noise(exptimes,magarr,fwhm, signal, noise, skyinst_bg,RN=12, DC=315):
	"""
	"""
	#RN = 12
	npix = np.pi * (fwhm/2)**2

	fig, ax = plt.subplots(1,1,figsize=(6,5),sharex=True,sharey=False)

	#for	itexp in np.arange(5):
	itexp = 3
	shotnoise = np.sqrt(signal[:,itexp])
	darknoise = np.sqrt(DC * exptimes[itexp] * npix)
	readnoise = np.sqrt(npix) * RN
	skyinstnoise = np.sqrt(skyinst_bg[:,itexp])

	ax.plot(magarr,noise[:,itexp],'k-',lw=3)
	ax.semilogy(magarr,shotnoise,'-')
	ax.semilogy(magarr,np.sqrt(shotnoise**2 + readnoise**2),'-.',label='+ Read Noise')
	ax.semilogy(magarr,np.sqrt(shotnoise**2 + readnoise**2 + darknoise**2),'--',label='+ Dark')
	ax.semilogy(magarr,np.sqrt(shotnoise**2 + readnoise**2 + darknoise**2 + skyinstnoise**2),'-.',label='+ Bkg')

	ax.legend()

def plot_noise_ratio(track_band,camera,exptimes,magarr,fwhm, signal, noise, skyinst_bg,RN=12, DC=315):
	"""
	"""
	#RN = 12
	npix = np.pi * (fwhm/2)**2

	fig, ax = plt.subplots(1,1,figsize=(6,5),sharex=True,sharey=False)

	#for	itexp in np.arange(5):
	itexp = 3
	shotnoise = np.sqrt(signal[:,itexp])
	darknoise = np.sqrt(DC * exptimes[itexp] * npix)
	readnoise = np.sqrt(npix) * RN
	skyinstnoise = np.sqrt(skyinst_bg[:,itexp])

	ax.semilogy(magarr,readnoise/shotnoise,'-.',label='Read Noise')
	ax.semilogy(magarr,darknoise/shotnoise,'--',label='Dark Noise')
	ax.semilogy(magarr,skyinstnoise/shotnoise,'-.',label='Inst. & Sky Bkg')

	ax.legend()
	ax.set_ylabel('Fraction to Photon Noise')
	ax.set_xlabel('H Mag')
	plt.subplots_adjust(left=0.15)
	plt.grid()
	plt.ylim(1e-5,3)
	plt.title('Tracking Band: %s \n Camera: %s'%(track_band,camera))

def plot_noise_two(track_band,magarr,noise1,noise2,label1,label2,itexp = 3):
	"""
	"""
	#RN = 12
	npix = np.pi * (fwhm/2)**2

	fig, ax = plt.subplots(1,1,figsize=(6,5),sharex=True,sharey=False)

	ax.semilogy(magarr,noise1,'-.',label=label1)
	ax.semilogy(magarr,noise2,'--',label=label2)

	ax.legend()
	ax.set_ylabel('Noise [e-]')
	ax.set_xlabel('H Mag')
	plt.subplots_adjust(left=0.15)
	plt.grid()
	#plt.ylim(1e-5,3)
	plt.title('Tracking Band: %s'%(track_band))

def plot_snr(exptimes,magarr,fwhm, signal, noise, RN=12, DC=315):
	"""
	"""
	tracking_requirement_pixel = 0.47
	npix = (np.pi * (fwhm)**2)
	readnoise =  np.sqrt(npix * RN**2)

	fig, ax = plt.subplots(1,1,figsize=(6,5),sharex=True,sharey=False)

	for	itexp in np.arange(5):
		p = ax.semilogy(magarr,signal[:,itexp]/noise[:,itexp],'-',label=str(exptimes[itexp]) + 's')
		ax.semilogy(magarr,np.sqrt(signal[:,itexp]),'--',c=p[0].get_color())		

	ax.set_xlim(2,15)
	ax.set_ylim(1,1e8)
	plt.grid()
	snr_needed = fwhm / tracking_requirement_pixel /np.pi
	ax.semilogy(magarr,snr_needed,'k--',lw=2)		
	ax.legend()
	return ax

def plot_data(exptimes,magarr,fwhm, signal, noise, centroid, RN=12, DC=315):
	"""
	"""
	#RN = 12
	npix = (np.pi * (fwhm)**2)
	readnoise =  np.sqrt(npix * RN**2)

	fig, ax = plt.subplots(3,1,figsize=(6,5),sharex=True,sharey=False)
	ax[0].plot(magarr,fwhm)

	for	itexp in np.arange(5):
		darknoise = np.sqrt(DC * exptimes[itexp] * npix)
		newnoise  = np.sqrt(strehl[:,itexp] * signal[:,itexp] + darknoise**2 + readnoise**2)

		ax[1].plot(magarr,noise[:,itexp],'-.')
		ax[1].semilogy(magarr,newnoise)
		ax[1].semilogy(magarr,readnoise,'k--')
		ax[1].semilogy(magarr,darknoise,'k--')
		ax[1].semilogy(magarr,np.sqrt((np.pi * (fwhm*0+np.min(fwhm)/2)**2) * RN**2),'g--')
		ax[1].semilogy(magarr,np.sqrt(((fwhm*0 +1)**2) * RN**2),'r--')
		
		#ax[2].plot(magarr,signal[:,itexp]/noise[:,itexp])
		SNR = strehl[:,itexp] * signal[:,itexp]/newnoise
		centroiderr = (1/np.pi) * fwhm/SNR
		ax[2].semilogy(magarr,centroiderr,label=exptimes[itexp])
		ax[2].semilogy(magarr,centroid[:,itexp],'--',label=exptimes[itexp])
	
	ax[2].legend()
	ax[2].semilogy(magarr,centroiderr*0 + 0.47)

def load_brown_dwarfs_AB():
    """
    load adam burgasser's file
    """
    xl = pd.ExcelFile('/Users/ashbake/Documents/Research/Projects/HISPEC/Tests/TrackingCamera/data/ucd_sheet_teff.xlsx')
    #xl.sheet_names
    first_sheet = xl.sheet_names[0]
    df = xl.parse(first_sheet)
    #df.head()
    teff = df['teff']
    hmag = df['H_2MASS']

    return  hmag.values,teff.values

def load_planets():
    planets_filename = '/Users/ashbake/Documents/Research/Projects/HISPEC/SNR_calcs/data/populations/confirmed_uncontroversial_planets_2023.03.08_14.19.56.csv'
    planet_data =  pd.read_csv(planets_filename,delimiter=',',comment='#')
    # add brown dwarfs!
    hmags = planet_data['sy_hmag']
    teffs = planet_data['st_teff']
    rvamps = planet_data['pl_rvamp']
    method = planet_data['discoverymethod']
    return hmags,teffs,method,rvamps

def plot_planets_ao_modes_strehl():
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
    best_ao_mode_arr   =  np.load('/Users/ashbake/Documents/Research/Projects/HISPEC/SNR_calcs/output/ao_modes/best_ao_mode_Hmag_temps_1000,1500,2300,3000,3600,4200,5800_modes_%s.npy'%ao_modes)
    magarr             = np.load('/Users/ashbake/Documents/Research/Projects/HISPEC/SNR_calcs//output/ao_modes/best_ao_mode_Hmag_mags.npy')
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
    camera = 'Hband'
    if camera=='cred2':
        mag_limit_h_cred2  =     [14.5,13.9,13.9,13.9,15,15]
        mag_limit_jhgap_cred2  = [12.3,12.6,13,13.6,14.2,15]
        mag_limit_jhplus_cred2 = [15.2,14.7,14.7,14.8,14.9,15]
        mag_limit_j_cred2      = [14.6,13.7,13.6,13.6,14.4,15]
        mag_limit_jhgap_TT_cred2   = [12.7,12.9,13.8,14.1,14.9,15]
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_cred2,'-',c='yellow',lw=2,label='JH Gap')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_TT_cred2,'--',c='yellow',lw=2,label='JH Gap + TT Star')
        ax.plot(teff_limit_jhgap,mag_limit_h_cred2,'-',c='orange',lw=2,label='H')
        ax.plot(teff_limit_jhgap,mag_limit_j_cred2,'-',c='limegreen',lw=2,label='J')
        ax.plot(teff_limit_jhgap,mag_limit_jhplus_cred2,'-',c='red',lw=2,label='J+H')
        #ax.set_title('Uninterrupted Science Limits')
        ax.set_title('CRED2 Tracking')
    if camera=='cred2_2':
        mag_limit_yjh_cred2  =     [14.5,13.9,13.9,13.9,15,15]
        mag_limit_yjh_cred2_rn  =     [14.5,13.9,13.9,13.9,15,15]
        mag_limit_yjh_cred2  =     [14.5,13.9,13.9,13.9,15,15]
        mag_limit_jhgap_cred2  = [12.3,12.6,13,13.6,14.2,15]
        mag_limit_jhplus_cred2 = [15.2,14.7,14.7,14.8,14.9,15]
        mag_limit_j_cred2      = [14.6,13.7,13.6,13.6,14.4,15]
        mag_limit_jhgap_TT_cred2   = [12.7,12.9,13.8,14.1,14.9,15]
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_cred2,'-',c='yellow',lw=2,label='JH Gap')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_TT_cred2,'--',c='yellow',lw=2,label='JH Gap + TT Star')
        ax.plot(teff_limit_jhgap,mag_limit_h_cred2,'-',c='orange',lw=2,label='H')
        ax.plot(teff_limit_jhgap,mag_limit_j_cred2,'-',c='limegreen',lw=2,label='J')
        ax.plot(teff_limit_jhgap,mag_limit_jhplus_cred2,'-',c='red',lw=2,label='J+H')
        #ax.set_title('Uninterrupted Science Limits')
        ax.set_title('CRED2 Tracking')
    if camera=='h2rg':
        mag_limit_jhgap_h2rg   = [12.5,12.7,13.2,13.7,14.2,15]
        mag_limit_j_h2rg       = [15.5,14.7,14.7,14.6,16,16]
        #mag_limit_jhgap_TTstar_h2rg   = [13.5,13.8,13.6,14.9,14.9,16]
        mag_limit_jplus_h2rg   = [15.5,14.8,14.9,14.9,16,16]
        mag_limit_jhplus_h2rg  = [15.8,15.6, 15.6, 15.6,16,16]
        mag_limit_h_h2rg       = [15.5,15.6,15.6,15.6,16,16]
        mag_limit_k_h2rg       = [15.5,15.5,15.6,15.6,16,16]
        mag_limit_jhgap_TT_h2rg    = [13.6,13.8,14.6,14.8,15,16]
        ax.plot(teff_limit_jhgap,mag_limit_h_h2rg,'-',c='orange',lw=2,label='H')
        ax.plot(teff_limit_jhgap,mag_limit_jplus_h2rg,'--',c='pink',lw=2,label='Jplus')
        ax.plot(teff_limit_jhgap,mag_limit_j_h2rg,'-',c='pink',lw=2,label='J')
        ax.plot(teff_limit_jhgap,mag_limit_jhplus_h2rg,'-',c='limegreen',lw=2,label='J+H')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_h2rg,'-',c='yellow',lw=2,label='JHgap')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_TT_h2rg,'--',c='yellow',lw=2,label='JHgap + Bright TT Star')
        ax.plot(teff_limit_jhgap,mag_limit_k_h2rg,'--',c='purple',lw=2,label='K')
        ax.set_title('H2RG Tracking')
    if camera=='TTstar':
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
    if camera=='Hband': # updated
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
    if camera=='uninterrupted': # updated
        mag_limit_jhgap_perfect   = [13.7,13.8,14.6,14.8,15,16]        
        mag_limit_jplus_perfect   = [15.8,15.5,15.4,15.55,16,16]
        mag_limit_h_perfect       = [15.8,15.3,15.6,15.6,16,16]
        ax.plot(teff_limit_jhgap,mag_limit_h_perfect,'-',c='orange',lw=2,label='H')
        ax.plot(teff_limit_jhgap,mag_limit_jplus_perfect,'-',c='limegreen',lw=2,label='Jplus')
        ax.plot(teff_limit_jhgap,mag_limit_jhgap_perfect,'--',c='yellow',lw=1.2,label='JHgap')
        ax.set_title('H2RG Tracking - TT Star')
    ax.plot([1000,6000],[15,15],'-',c='gray',lw=1)
    #ax.plot(teff_limit_jhgap,mag_limit_j_h2rg,'.-',c='cyan',lw=1)
    #ax.plot(teff_limit_jhgap,mag_limit_j_cred2,'.-',c='red',lw=1)

    ax.set_ylim(5,16)
    ax.legend(fontsize=10,loc=4)

    ntotal_bd = len(np.where((hmags_bd < 15) & (teffs_bd >1000) & (teffs_bd < 5800))[0])
    ntotal_pl = len(np.where((hmags < 15) & (teffs >1000) & (teffs < 5800))[0])
    nsub_bd_h2rg = len(np.where((hmags_bd < 13.5) & (teffs_bd >1000) & (teffs_bd < 5800))[0])
    nsub_bd_cred = len(np.where((hmags_bd < 12.5) & (teffs_bd >1000) & (teffs_bd < 5800))[0])
    nsub_bd_nocando = len(np.where((hmags_bd < 15) &(hmags_bd > 14.75) & (teffs_bd >1000) & (teffs_bd < 5800))[0])

        plt.savefig('/Users/ashbake/Documents/Research/Projects/HISPEC/SNR_calcs/output/trackingcamera/ao_mode_landscape_brown_dwarfs_trackingcutoff.png',dpi=300)
    else:
        plt.savefig('./output/ao_modes/ao_mode_landscape_brown_dwarfs.png')

if __name__=='__main__':
	#########
	# plot
	tracking_requirement_pixel = 0.47
	teff       = 2300
	track_band = 'H'#['Jplus','JHgap','H']
	camera     = 'h2rg'
	datapath   = './output/_run_20230427_%s/'%camera
	ao_mode    = 'LGS_STRAP_130' if teff>=3000 else 'LGS_100%s_130'#'LGS_STRAP_130'#
	RN = 12#30
	DC = 0.8#315

	out = load_data(camera,ao_mode,track_band,teff,datapath)
	exptimes,magarr,centroid,fwhm,signal,noise,strehl, centroid,snr,skyinst_bg = out
	centroid2 = np.pi *fwhm[0]/signal

	plot_noise(exptimes,magarr,fwhm, signal, noise, skyinst_bg, RN=RN,DC=DC)
	#plot_noise_ratio(track_band,camera, exptimes,magarr,fwhm, signal, noise, skyinst_bg, RN=RN,DC=DC)

	ax  = plot_snr(exptimes,magarr,fwhm, signal, noise, RN=RN,DC=DC)
	ax.set_title('%s %s band' % (camera, track_band))
	
	plot_data(exptimes,magarr,fwhm, signal, noise,centroid, RN=RN,DC=DC)


	##########
	itexp =2

	out1 = load_data('cred2',ao_mode,track_band,teff,'./output/_run_20230427_cred2/')
	out2 = load_data('h2rg',ao_mode,track_band,teff,'./output/_run_20230427_h2rg/')
	noise1 = np.sqrt(np.pi * (out1[3]/2)**2) * 25
	noise2 = np.sqrt(np.pi * (out2[3]/2)**2) * 12
	plot_noise_two(track_band,out1[1],noise1,noise2,'RN cred2','RN h2rg',itexp = 3)
	
	noise1 = np.sqrt(np.pi * (out1[3]/2)**2 * 230 * exptimes[itexp])
	noise2 = np.sqrt(np.pi * (out2[3]/2)**2 * 0.8 * exptimes[itexp])
	plot_noise_two(track_band,out1[1],noise1,noise2,'DN cred2','DN h2rg',itexp = 3)


	plot_noise_two(track_band,out1[1],out1[3],out2[3],'FWHM cred2','FWHM h2rg',itexp = 3)
	plot_noise_two(track_band,out1[1],out1[6],out2[6],'strehl cred2','strehl h2rg',itexp = 3)
	plot_noise_two(track_band,out1[1],3.14 * (out1[3]/2)**2,3.14 * (out2[3]/2)**2,'npix cred2','npix h2rg',itexp = 3)


