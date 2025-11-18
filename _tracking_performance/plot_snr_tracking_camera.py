# calc signal to noise
# max for the calcium H&K
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from datetime import date

font = {'size'   : 14}
matplotlib.rc('font', **font)

sys.path.append('/Users/ashbake/Documents/Research/ToolBox/specsim/')

from specsim.objects import load_object
from specsim.load_inputs import fill_data, load_filter,get_band_mag
from specsim.functions import *

plt.ion()

from matplotlib.ticker import (AutoMinorLocator)

plt.ion()
font = {'size'   : 14}
matplotlib.rc('font', **font)
plt.rcParams['font.size'] = '14'
plt.rcParams['font.family'] = 'sans'
plt.rcParams['axes.linewidth'] = '1.3'
fontname = 'DIN Condensed'

def get_mag_limit_per_exptime(exptimes,magarr,centroid):
	"""
	Computes magnitude limit at which tracking requirement is met
	by interpolating centroid error vs magnitude curve

	inputs
	------
	exptimes: array
		must be in seconds
	magarr: array
		magnitude array
	centroid: array
		centroid error in pixels in len(magarr) x len(exptimes) array
	
	outputs
	-------
	mag_requirement: array
		magnitude limit for each exposure time
	"""
	tracking_requirement_pixel = 0.4 # pixels

	mag_requirement = np.zeros((len(exptimes)))
	for j, exptime in enumerate(exptimes):
		# interpolate curve on higher density grid
		ifit = np.where((~np.isnan(centroid[:,j])) & (centroid[:,j] !=0)& (centroid[:,j] <2))[0]
		if len(ifit) < 3: mag_requirement[j]=np.nan; continue
		magdense = np.arange(np.min(magarr[ifit]),np.max(magarr[ifit]),0.01)
		interp   = interpolate.interp1d(magarr[ifit],centroid[:,j][ifit],kind='quadratic',bounds_error=False,fill_value='extrapolate')
		cendense = interp(magdense)
		try:
			ireq1 = np.where(cendense > tracking_requirement_pixel)[0][0]
			ireq2 = np.where(cendense < tracking_requirement_pixel)[0][-1]
			ireq=ireq1 if np.abs(ireq1-ireq2) <2 else ireq2
		except:
			ireq = 0
		mag_requirement[j] = magdense[ireq]
		if ireq==0 and np.max(cendense) > tracking_requirement_pixel : mag_requirement[j]=np.nan
		if ireq==0 and np.max(cendense) < tracking_requirement_pixel : mag_requirement[j]=np.max(magdense)

		#plt.figure()
		#plt.plot(magarr,centroid[:,j])
		#plt.plot(magarr[ifit],centroid[:,j][ifit])
		#plt.plot(magdense,cendense)
		#plt.title(str(exptime) + '  ' + str(mag_requirement[j]))

	return mag_requirement


def plot_results_mag_req(teff, datapath, track_bands=['JHgap','Jplus','Hplus','JHplus'],field_r=0):
	"""
	"""
	exptimes = np.load(datapath + 'exptimes.npy')
	magarr   = np.load(datapath + 'magarr.npy')

	fig, axs = plt.subplots(1,1,figsize=(7,7),sharex=True,sharey=True)
	for ii,track_band in enumerate(track_bands):
		savename = 'cendata_%s_ao_auto_band_%s_teff_%s_%smag_fieldr_%sarcsec.npy' %('h2rg',track_band,teff,'H',field_r)
		try: centroid = np.load(datapath + 'centroid_%s'%(savename))
		except: continue
		try:	mag_req = get_mag_limit_per_exptime(exptimes,magarr,centroid)
		except: mag_req = np.zeros_like(exptimes)
		axs.semilogx(exptimes,mag_req,label=track_band)

	axs.set_xlabel('Exposure Time [s]')
	plt.subplots_adjust(bottom=0.15,top=0.85,right=0.9)
	axs.axhline(15,c='k',linestyle='--')
	axs.set_ylabel('H Limit')
	axs.grid(color='gray',linestyle='-',alpha=0.2)
	axs.set_ylim(11,17)
	axs.fill_between(exptimes,15,y2=22,facecolor='green',alpha=0.2)
	axs.set_xticks(exptimes)
	axs.set_xticklabels(exptimes,fontsize=12)
	yticks = [10,11,12,13,14,15,16,17,18]
	axs.set_yticks(yticks)
	axs.set_yticklabels(yticks,fontsize=12)
	axs.set_title('Tracking Performance\nT$_{eff}$=' + str(teff) + ' K')
	plt.legend()

	savepath = os.path.join(*datapath.split('/')[0:-1]) + '/_plots/'
	if not os.path.isdir(savepath): os.makedirs(savepath)

	plt.savefig(savepath + 'tracking_H_magnitude_limits_%sK.png'%(teff))



def plot_frame_rates_magnitude(field_r,track_bands=['J','Jplus','H','JHgap'],datapath = './output/centroid_arrs/'):
	"""
	must edit savename and link to folder in centroid/data/ 

	because i move files around to organized
	data in _run_20230221_PDR_offaxis_comparison
	"""
	exptimes = np.load(datapath + 'exptimes.npy')
	magarr   = np.load(datapath + 'magarr.npy')
	maglimits=[11.9,12.9,13.9,14.9,15.9]
	mag_crossover = np.zeros(len(track_bands))

	#colors = ['steelblue','orange','green','purple']
	syms = ['o','s','d','.']
	fig, ax = plt.subplots(1,1,figsize=(6,5),sharex=True,sharey=True)

	for jj,track_band in enumerate(track_bands):
		framerates =np.zeros(len(maglimits))
		for ii,maglimit in enumerate(maglimits):
			savename = 'cendata_h2rg_ao_auto_band_%s_teff_%s_%smag_fieldr_%sarcsec.npy' %(track_band,teff,'H',field_r)
			centroid = np.load(datapath + 'centroid_%s'%(savename))
			try:	mag_req = get_mag_limit_per_exptime(exptimes,magarr,centroid)
			except: mag_req = np.zeros_like(exptimes)
			fwhm = np.load(datapath + 'fwhm_%s'%savename)
			snr  = np.load(datapath + 'snr_%s'%savename)
			exptimes_dense = np.arange(np.min(exptimes),np.max(exptimes),0.00001)
			ftemp = interpolate.interp1d(exptimes,mag_req)
			mag_req_dense = ftemp(exptimes_dense)
			try:
				print(track_band, 1/exptimes_dense[np.where(mag_req_dense>maglimit)][0])
				framerates[ii] =  1/exptimes_dense[np.where(mag_req_dense>maglimit)][0]
			except:
				print(track_band)
				#plt.plot(exptimes,mag_req,syms[ii],c=p[0].get_color())
		p = ax.scatter(maglimits,framerates,marker=syms[0],label=track_band)
		p = ax.plot(maglimits,framerates,'k--',lw=0.8)
		# save frame rate where magnitude crosses 27 Hz
		icrossover = np.where(1/exptimes_dense<27)[0]
		if len(icrossover)>0:
			mag_crossover[jj] = mag_req_dense[icrossover[0]]
		else:
			mag_crossover[jj] = 16

	ax.set_title('T=%sK'%(teff))
	ax.legend(fontsize=10)

	ax.set_xlabel('H Mag')
	plt.subplots_adjust(bottom=0.15,left=0.15,top=0.85,right=0.9)
	ax.set_ylabel('Max Frame Rate (Hz)')
	ax.grid(color='gray',linestyle='--',alpha=0.5,dash_joinstyle='miter')
	ax.set_yscale('log')
	yticks = [10,27,100,500,1000]
	ax.set_yticks(yticks)
	ax.set_yticklabels(yticks,fontsize=14)
	ax.set_ylim(1,1200)
	ax.fill_between([11.7,16.3],0,27,alpha=0.3,color='r',zorder=-100)
	ax.plot([15,15],[1,1200],'k',lw=2)

	savename = 'frame_rate_Teff_%sK.png'%(teff)
	savepath = os.path.join(*datapath.split('/')[0:-1]) + '/_plots/'
	if not os.path.isdir(savepath): os.makedirs(savepath)

	plt.savefig(savepath + savename)

	return mag_crossover




if __name__=='__main__':
	# plot magnitude limit results for each teff and field radius
	
	#### UPDATE ME
	track_bands = ['JHgap','Jplus','Hplus','JHplus']
	datapath = './output/centroid_arrs_rad3/'
	field_r=3
	teffs = np.array([1000,1500,2300,3000,3600,4200,5800])
	###########

	# find mag crossovers and make all plots
	mag_crossover = np.zeros((len(track_bands),len(teffs)))
	for i, teff in enumerate(teffs):
		plot_results_mag_req(teff,datapath,track_bands=track_bands,field_r=field_r)
		mag_crossover[:,i] = plot_frame_rates_magnitude(field_r,track_bands=['JHgap','Jplus','Hplus','JHplus'],datapath = datapath)
		# save magnitude where frame rate dips below 27 Hz
	
	# reformat maglimits for easy use in plot_planets_tracking_performance.py
	maglimits = {}
	for ii,track_band in enumerate(track_bands):
		maglimits[track_band] = []
		for jj,teff in enumerate(teffs):
			maglimits[track_band].append(mag_crossover[ii,jj])
	np.save(datapath + '_plots/maglimits',maglimits)


