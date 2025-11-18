# plot planet and brown dwarf population and
# show limits for each tracking band
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
import pandas as pd

font = {'size'   : 14}
matplotlib.rc('font', **font)

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
	xl = pd.ExcelFile('/Users/ashbake/Documents/Research/_DATA/populations/browndwarfs/ucd_sheet_teff.xlsx')
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
	# plot magnitude limit results for each teff and field radius
	band = 'H'
	best_ao_mode_arr   = np.load("/Users/ashbake/Documents/Research/Projects/HISPEC/Tests/FEI_TrackingCamera/_tracking_performance/_postPDR_run/data/ao_modes/best_ao_mode_Hmag_temps_1000,1500,2300,3000,3600,4200,5800_modes_['LGS_100J_130', 'SH'].npy")
	magarr             = np.load('/Users/ashbake/Documents/Research/Projects/HISPEC/Tests/FEI_TrackingCamera/_tracking_performance/_postPDR_run/data/ao_modes/best_ao_mode_Hmag_mags.npy')

	temp_arr   = [1000,1500,2300,3000,3600,4200,5800]
	extent     = (np.min(temp_arr),np.max(temp_arr),np.min(magarr),np.max(magarr))

	# plot contours of it
	hmags, teffs,method,rvamps = load_planets()
	i_image = np.where(method=='Imaging')[0]
	i_ms= np.where(rvamps<3)[0]
	hmags_bd, teffs_bd = load_brown_dwarfs_AB()

	def fmt(x):
		"format labels for contour plot"
		return str(np.round(100*x)) + '%'

	##########
	fig, ax = plt.subplots(1,1,figsize=(6,5))
	#ax.imshow(best_ao_mode_arr.T,aspect='auto',origin='lower',\
	#		interpolation='None',cmap='gray',\
	#		extent=extent,vmin=0,vmax=2) # THIS ASSUMES TEMP ARRAY IS REGULAR WHEN ITS NOT
	X,Y=np.meshgrid(temp_arr, magarr)
	ax.pcolor(X,Y,best_ao_mode_arr.T,cmap='nipy_spectral_r',shading='auto')

	ax.scatter(teffs,hmags,marker='.',s=65,c='m',ec='purple',alpha=0.8,label='Confirmed Planet Hosts')
	#ax.scatter(teffs[i_image],hmags[i_image],marker='s',s=40,c='brown',ec='k',alpha=1,label='Imaged Planets')
	ax.scatter(teffs_bd,hmags_bd,marker='d',s=32,c='darkcyan',ec='b',alpha=.8,label='Brown Dwarfs')
	#ax.scatter(hmags_bd,teffs_bd,marker='.',c='c',alpha=0.5,label='Brown Dwarfs')
	ax.set_ylabel('H Mag')
	ax.set_xlabel('T$_{eff}$ (K)')
	ax.set_title('HISPEC AO Mode Landscape')
	ax.set_ylim(0,15.9)
	ax.set_xlim(1000,5800)
	plt.subplots_adjust(bottom=0.15,hspace=0.1,left=0.15,right=0.85)

	# load mag limits of tracking bands
	teff                  = [1000,1500,2300,3000,3600, 4200,5800] 
	maglimits_2ishmm      = np.load('./output/centroid_arrs_PWV2.28mm/_plots/maglimits.npy',allow_pickle=True).item()
	maglimits_1mm         = np.load('./output/centroid_arrs_PWV1mm/_plots/maglimits.npy',allow_pickle=True).item()
	maglimits_offaxis     = np.load('./output/centroid_arrs_rad3/_plots/maglimits.npy',allow_pickle=True).item()
	track_bands       = list(maglimits.keys())
	
	# plot JHgap
	colors = ['cyan','red','limegreen']
	track_band = 'JHgap'
	#ax.plot(teff,maglimits_1mm[track_band],'-',c=colors[0],lw=2,label=track_band + '(PWV=1mm)')
	#ax.plot(teff,maglimits_2ishmm[track_band],'--',c=colors[1],lw=2,label=track_band + '(PWV=2.3mm)')
	#ax.plot(teff,maglimits_offaxis[track_band],'-',c=colors[2],lw=2,label=track_band + '(PWV=1mm, 3" offaxis)')

	# plot rest of bands
	colors2 = ['cyan','orange','magenta', 'yellow']
	for i,track_band in enumerate(track_bands):
		ax.plot(teff,maglimits_offaxis[track_band],'--',c=colors2[i],lw=2,label=track_band + '(PWV=2.3mm)')

	ax.set_title('HISPEC APIC JH Gap Tracking')

	ax.plot([1000,6000],[15,15],'-',c='gray',lw=1)

	ax.legend(fontsize=10,loc=4)
	plt.savefig('./output/ao_mode_landscape_trackingcutoff.png',dpi=300)

	ntotal_bd = len(np.where((hmags_bd < 15) & (teffs_bd >1000) & (teffs_bd < 5800))[0])
	ntotal_pl = len(np.where((hmags < 15) & (teffs >1000) & (teffs < 5800))[0])
	nsub_bd_h2rg = len(np.where((hmags_bd < 13.5) & (teffs_bd >1000) & (teffs_bd < 5800))[0])
	nsub_bd_cred = len(np.where((hmags_bd < 12.5) & (teffs_bd >1000) & (teffs_bd < 5800))[0])
	nsub_bd_nocando = len(np.where((hmags_bd < 15) &(hmags_bd > 14.75) & (teffs_bd >1000) & (teffs_bd < 5800))[0])


