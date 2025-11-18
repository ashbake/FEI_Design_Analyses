# calc the positional error for a
# gaussian PSF with noise or no noise
#
# see reference https://www.astro.ucla.edu/~ghezgroup/meeting_papers/fritz_limiting.pdf
# for guidance on theory
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import scipy
import pandas as pd
from datetime import date

from photutils.aperture import ApertureStats, CircularAperture
from photutils.centroids import (centroid_1dg, centroid_2dg,
                                 centroid_com, centroid_quadratic)
from photutils.aperture import aperture_photometry

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

NIM = 20

def twoD_Gaussian(x,y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """
    """
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    return offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))

def add_noise(psf_im,seed):
	"""
	given 2D image, produce a photon + read noise + dark frame for it
	"""
	#seed       = 42
	# shot noise
	rs         = np.random.RandomState(seed)
	shot_noise_frame = rs.poisson(psf_im) - psf_im

	# read noise
	read_noise = 12 # electrons - cds noise
	read_noise_frame = rs.normal(scale=read_noise, size=shot_noise_frame.shape)

	# total noise frame
	total_noise_frame = read_noise_frame + shot_noise_frame

	# total scienc eimage with noise
	electrons_out = total_noise_frame + psf_im

	return electrons_out


def get_peak_cen(data,ploton=False):
    """
    convolve image with  gaussian of expected fwhm to smooth out for peak finding

    inputs
    ------
    data - 2d array
    """
    xtarr,ytarr = np.meshgrid(np.arange(NIM), np.arange(NIM)) # x y arrays corresponding to pixel position of 10x10 sub box
    gaus = twoD_Gaussian(xtarr,ytarr, 1,NIM//2, NIM//2, 4/2.355, 4/2.355, 0, 0)
    data_conv  = scipy.signal.convolve2d(data, gaus, mode='same', boundary='fill', fillvalue=0)

    trim       = 2
    sub_data   = data_conv[trim:-1*trim, trim:-1*trim]

    # fit for centroid
    #x1,y1 = centroid_2dg(data)
    #x2,y2 = centroid_1dg(data)
    x3,y3 = centroid_quadratic(data)
    #x4,y4 = centroid_com(data)
    
    #xcen,ycen = x3+trim, y3+trim
    if ploton:
        plt.figure()
        plt.imshow(data)
        plt.plot(x3,y3,'kx',label='QDR')
        #plt.plot(x1,y1,'kx')
        
    return x3,y3#xcen,ycen


def create_sci_psf(amplitude):
	"""
	make guassian 2d psf and convolve with optical radius for science beam
	and for ghost beam

	Add the two images, add noise if turned on, then compute peak center error
	"""
	# make 2D gaussian and fit center
	pixel_pitch = 18#um/pix

	x,y       = np.meshgrid(np.arange(NIM), np.arange(NIM)) # x y arrays corresponding to pixel position of 10x10 sub box
	#amplitude = peaksnr**2  # matters if add noise
	xo,yo     = NIM//2, NIM//2
	sigma_x, sigma_y = 4.1/2.355,4.1/2.355 # fwhm is 4 pixels in K band for hispec
	theta, offset    = 0, 0

	# perfect PSF no convolution
	psf   = twoD_Gaussian(x,y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)

	return psf



def create_sci_noise_manifestations(sci_frame,ploton=False):
	"""
	create many sci frames with different noise
	"""
	niter = 200
	xcens,ycens= np.zeros(niter),np.zeros(niter)
	psfs = np.zeros((niter,NIM,NIM))
	for seed in np.arange(niter):
		noisy_psf = add_noise(sci_frame,seed=seed)
		xcens[seed],ycens[seed]  = get_peak_cen(noisy_psf,ploton=False)
		psfs[seed,:,:] = noisy_psf

	# compute distance then take std to get error on centroid in pixels
	r = np.sqrt((xcens-NIM//2)**2 + (ycens-NIM//2)**2)
	centroid_error =  np.nanstd(r)

	# compute SNR 
	# 1 - compute snr based on area of flux in radius fwhm/2 compared to std in same area in the noise frame
	error_image = np.std(psfs,axis=0)
	avg_img     = np.mean(psfs,axis=0)

	aper      = CircularAperture([len(avg_img)//2, len(avg_img)//2], 4.1/2)
	ap_stats1 = ApertureStats(avg_img, aper)
	flux = ap_stats1.sum
	ap_stats2 = ApertureStats(error_image**2, aper)
	std       = np.sqrt(ap_stats2.sum) # add in quadrature pixels
	snr       = flux/std

	# 2 - peak SNR
	#snr_peak  = np.max(psf_im) /np.sqrt(np.max(psf_im) + read_noise**2)
	
	# 3 - total but 50% of flux bc of area under fwhm related to npix
	#npix = aper.area # should be 13.2 for 4.1 radius
	#snr_total = 0.5 * np.sum(psf_im) /np.sqrt(0.5 * np.sum(psf_im) + npix * read_noise**2)

	if ploton:
		plt.imshow(psfs[0])
		ap_patches = aper.plot(color='white', lw=2,
                           label='Photometry aperture')
		plt.colorbar()
		plt.title('Example PSF image \n (SNR=%s)'%round(snr))

		plt.figure()
		plt.imshow(error_image)
		ap_patches = aper.plot(color='white', lw=2,
                           label='Photometry aperture')
		plt.colorbar()
		plt.title('Error image')

	return centroid_error, snr 


def get_track_req():
	# compute requirement. requirement is 0.2lambda/D in y band
	fratio = 35
	D      = 10
	pixel_pitch = 18
	platescale_arcsec_um = 206265 / fratio / (D * 10**6) #arc/um
	platescale_arcsec_pix = platescale_arcsec_um * pixel_pitch

	yband_wavelength       = 1020 # nm, center of y band
	tracking_requirement_arcsec = 206265 * 0.2 * yband_wavelength / (D*10**9) 
	tracking_requirement_pixel  = tracking_requirement_arcsec/platescale_arcsec_pix

	return tracking_requirement_pixel

def get_track_exp(snr,fwhm):
	"""
	get tracking expectation
	"""
	return 1/np.pi * fwhm/snr


def calc_peak_to_total_flux_ratio(amplitude):
	"""
	If peak flux is X what is total flux for a 4.1 px FWHM 2D gaussian?

	answer is 0.05 (amp is 1/20th of total flux)

	If flux is X in 4pix diameter aperture  for 4.1 pix FWHM 2D guassian,
	what is the peak flux?
	50% of total flux is in the 4pix diameter aperture - thus answer
	is 1/10th
	"""
	# make 2D gaussian and fit center
	pixel_pitch = 18 #um/pix
	nim         = 30   # size of image

	x,y       = np.meshgrid(np.arange(nim), np.arange(nim)) # x y arrays corresponding to pixel position of 10x10 sub box
	#amplitude = peaksnr**2  # matters if add noise
	xo,yo     = nim//2, nim//2
	fwhm = 4.1 # pixels  - size of PSF
	sigma_x, sigma_y = fwhm/2.355,fwhm/2.355 # fwhm is 4 pixels in K band for hispec
	theta, offset    = 0, 0

	# perfect PSF no convolution
	perfect_psf   = twoD_Gaussian(x,y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)

	# total flux
	total_flux = np.sum(perfect_psf)

	print(amplitude/total_flux)


if __name__=='__main__':
	#load inputs
	amplitudes = np.logspace(2,5,15)
	snrs       = np.zeros(len(amplitudes))
	rs       = np.zeros(len(amplitudes))

	for i,amplitude in enumerate(amplitudes):
		sci_frame     = create_sci_psf(amplitude)
		rs[i],snrs[i] = create_sci_noise_manifestations(sci_frame)
		print(i)

	
	#########
	# PLOTT
	plt.figure()
	plt.loglog(snrs, rs,'-o',label='SNR in 4.1 pix Diameter')

	plt.xlabel('SNR')
	plt.ylabel('Error [pix]')
	
	track_expect = get_track_exp(snrs,fwhm=4.1)
	plt.plot(snrs, track_expect,c='k',ls='-',label='Expectation')
	
	track_req = get_track_req()
	plt.axhline(track_req,c='k',ls='--',label='requirement')

	plt.legend(fontsize=10)
	plt.subplots_adjust(bottom=0.15,left=0.2)
	plt.grid()
	plt.ylim(0.3e-3,1)
	plt.title('Quadratic Centroiding Method')
	plt.savefig('PSFcentroid_error_vs_SNR_nim_%s.png'%NIM)

	#




