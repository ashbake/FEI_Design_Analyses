########################################
# Script to display the PSF position of the CRED2 data
# as analyzed using Reducer_CRED2.py script.
# 
# Script can generate two plots, a simpler 2-axis one
# with just the PSF scatter and a sample frame or a 
# larger one that adds 2 more axes for a histogram, 
# and PSD analysis.
#
# This script should be run on the fits files output by
# the Reducer_CRED2.py file which does the actual PSF finding.
# The finding and plotting were split up to improve usability.
#
#
# NOTE: Updated in January 2023. This is the definitive version as of 01/19/2023
# - Now attempts to recognize bad PSF fits and remove them from the analysis

########################################


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy import fft, signal

#-- Filenames and stuff
# Cred2 Stuff
fl_creddir = f'data_20240417/'
fl_bkgddir = fl_creddir + 'biases/'

specify_bkgd = False    # Flag to denote if specific background should be used. If so, provide bkgd file below
#(OPTIONAL) Pull a more recent BG file
#fl_bkgd = fl_bkgddir+'bias_154047_0400.00_00.00249_01_-40.00000_000_000_000_000.fits'

#--File indices (fits file numbers to plot)
fl2an = np.array([3,5,6,7])      # Should be a list (or numpy array) For single frame, just give a single element: [1]np.arange(168, 168+1)     # Should be a list (or numpy array) For single frame, just have a single element: [1]
#-- Saving option
isSaveFig = True

#-- Misc. analysis parameters
fl_credpre = 'cred2_0'
fl_credpost = '_proc'
PLATESCALE = 7.240      # mas/pix
lamOD = 2.2/(10e6) * 206265*1000            # Lambda over D of Keck in [mas]
plt.ion()


#######################################
#-- Helper function definitions
#######################################
#-- Define a function to load data
def loader(fl_num, pth=fl_creddir, pre=fl_credpre, post='', lev2pull=0):
    flnm = pth+pre+f'{fl_num:03d}'+post
    with fits.open(flnm+'.fits') as f:
        data = f[lev2pull].data
        head = f[lev2pull].header
        
    return (data, head, flnm)
    
#-- Define function to plot the results
def plotter(flnum, ismnzero, ispltmean, ax2tit, isBigPlot=True, isPrintParams=True):
   # Create/format figure and axes
    if isBigPlot:
        fig, axs = plt.subplots(2,2, figsize=(14.3,10))
        ((ax1, ax2), (ax3, ax4)) = axs
    else:
        fig, axs = plt.subplots(1,2, figsize=(9.35, 5.56))
        ax1, ax2 = axs
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.948, top = 0.84)
    # Load the data
    crd_dat, hdr, _ = loader(flnum,post=fl_credpost)
    im_dat, _, flnm = loader(flnum)
    if not specify_bkgd:
        # Find the corresponding background
        fl_bkgd = hdr['BIAS'].split('/')[-1]
        fl_bkgd = fl_bkgddir + fl_bkgd
    # Pull the background
    with fits.open(fl_bkgd) as f:
        bkgd = f[0].data

    #---- Make the basic plots
    # Plot the sample image
    toplot = im_dat[0,1:,:]-bkgd[1:,:]
    toplot[np.isnan(toplot)] = 0
    im = ax1.imshow(toplot, origin='lower')
    # Crop the image to a window around the first PSF coordinate
    coord1 = crd_dat[0,:2]
    pltwin = 5
    ax1.set(xlim=(coord1[1]-pltwin, coord1[1]+pltwin), ylim=(coord1[0]-1-pltwin, coord1[0]-1+pltwin))
    plt.colorbar(im, ax=ax1)
    ax1.scatter(crd_dat[:,1], crd_dat[:,0]-1, c=np.arange(0, crd_dat.shape[0]), marker='+') # -1 to account for tag pixel row which is removed above
    ax1.set_xlabel('[pix]')
    ax1.set_ylabel('[pix]')
    ax1.set_title('Sample Frame')
    # Set rows with bad/failed fits to nan so that means and stds below aren't biased by the bad entries
        # Check for cases where row-value is negative since our INVALID_PSF values are -9999 by default
    inv_psf_rows = crd_dat[:,0] < 0
        # Check for cases where the FWHM is unreasonable. KPIC design is ~4.2pix so check much smaller or larger
    inv_psf_fwhm = np.logical_or(crd_dat[:,2] < 2, crd_dat[:,2] > 20)
        # Now filter those rows out by setting them to nan
    inv_psf = np.logical_or(inv_psf_rows, inv_psf_fwhm)
    crd_dat[inv_psf] = np.nan
    if ismnzero:
        # Shift coordinates so that they are mean zero 
        goal_row, goal_col = crd_dat[:,0].nanmean(), crd_dat[:,1].nanmean()
    else:
        # Shift coordinates so that they are centered w.r.t commanded goal
        goal_col, goal_row = hdr['DISGOALX'], hdr['DISGOALY']
    crd_dat[:,0] -= goal_row 
    crd_dat[:,1] -= goal_col
    # Convert to mas
    crd_dat[:,:2] *= PLATESCALE 
    # Compute Metrics
    rstd, cstd = np.nanstd(crd_dat[:,0]), np.nanstd(crd_dat[:,1])
    rmen, cmen = np.nanmean(crd_dat[:,0]), np.nanmean(crd_dat[:,1])
    # Plot the PSF positions over time
    markerSize = 2
    ax2.scatter(crd_dat[:,1], crd_dat[:,0], markerSize, c=np.arange(0,crd_dat.shape[0]))
    ax2.grid()
    ax2.set_xlabel('X-displacement [mas]')
    ax2.set_ylabel('Y-displacement [mas]')
    ax2.set_title('PSF Position'+ax2tit)
    # Plot some spatial markers for reference
    ax2.axhline(0, c='black')      # Crosshair
    ax2.axvline(0, c='black')      # Crosshair
    if ispltmean:
        # Include a crosshair on the average position as requested
        ax2.scatter(cmen, rmen, 40, c='red', marker='+')
    ellipse = Ellipse(xy=(cmen, rmen), width = 2*cstd, height=2*rstd, edgecolor='r', fc = 'None', lw = 2)
    ax2.add_patch(ellipse)
    ellipse = Ellipse(xy=(0,0), width = 1*lamOD, height=1*lamOD, edgecolor='black', fc = 'None', lw = 1, linestyle='--')
    ax2.add_patch(ellipse)
    ellipse = Ellipse(xy=(0,0), width = 2*lamOD, height=2*lamOD, edgecolor='black', fc = 'None', lw = 1, linestyle='--')
    ax2.add_patch(ellipse)
    axlim = lamOD #mas
    ax2.set_xlim([-axlim, axlim]); ax2.set_ylim([-axlim, axlim]) 
    ax2.set_aspect('equal')

    if isBigPlot:
        #---- Make the remaining plots as requested
        # Plot histogram of T/T Residuals
        _,bins,_ = ax3.hist(x=crd_dat[:,0], bins='auto', color='b', alpha=0.5, rwidth=0.85, label='Y-axis (STD = %0.2f mas)'%rstd)
        _,bins,_ = ax3.hist(x=crd_dat[:,1], bins=bins, color='r', alpha=0.5, rwidth=0.85, label='X-axis (STD = %0.2f mas)'%cstd)
        ax3.set_xlabel('Tip/Tilt Jitter [mas]')
        ax3.set_ylabel('Samples')
        ax3.set_title('KPIC Tip-Tilt Jitter at CRED2')
        ax3.legend()

        # Get the timestamp data needed for frequency analysis
        tstamps, _, _ = loader(flnum, lev2pull=1)
        tstamps -= tstamps[0]
        
        # Since the welch function doesn't like NaN's, set them as equal to the previous value
        inv_psf_inds = np.where(inv_psf)[0]
        # Do replacement iteratively (one-by-one) to avoid issues from consequitive bad frames
        for inv_ind in inv_psf_inds:
            if inv_ind == 0:
                # If the first frame is bad, then we can't replace w/ previous so let's try with the next frame
                crd_dat[inv_ind] = crd_dat[inv_ind+1]
                continue
            crd_dat[inv_ind] = crd_dat[inv_ind-1]

        # Get PSD using Welch's method
        sample_rate = 1/np.mean(np.diff(tstamps))   # [Hz] sampling from average timedelta
            # Set welch window as ~1/4 of total datapoints and round to closest lower power of 2
        welch_NperSeg = 1 << int(np.log2(crd_dat.shape[0]/4))
        rFreq, rPSD = signal.welch(crd_dat[:,0], sample_rate, nperseg=welch_NperSeg, scaling='density')
        cFreq, cPSD = signal.welch(crd_dat[:,1], sample_rate, nperseg=welch_NperSeg, scaling='density')
        
        # Plot the PSD
        ax4.loglog(rFreq, rPSD, '-o', label='Y Jitter')
        ax4.loglog(cFreq, cPSD, '-o', label='X Jitter')
        ax4.set_xlabel('Frequency [Hz]')
        ax4.set_ylabel('PSD [$mas^2$/Hz]')
        ax4.set_title('PSD - Win Size = %d - Samp Rate = %d [Hz]'%(welch_NperSeg, sample_rate))
        ax4.legend()

    # Compute number of invalid frames removed
    ninv = np.sum(inv_psf)
    print('File #{:d}: Frames = {:d}  |  FPS = {:0.1f}  |  Tint = {:0.6f}  |  NDR = {:d} | N_inv = {:d}'.format(flnum, crd_dat.shape[0], hdr['FPS'], hdr['TINT'], hdr['NDR'], ninv))

    return fig, axs, {'rstd': rstd, 'cstd': cstd, 'rmen':rmen, 'cmen':cmen, 'goal_row':goal_row, 'goal_col':goal_col, 'flnm':flnm, 'ninv':ninv}, hdr


for ind, fl in enumerate(fl2an):
    fig, axs, vals, hdr = plotter(flnum=fl, ismnzero=False, ispltmean=True, ax2tit='(Goal-centered)')
        
    # Format suptitle title
    tit= 'File #%d - Goal=(%0.2f, %0.2f) [pix]\n'%(fl, vals['goal_col'], vals['goal_row'])
    if hdr['ISTRK'] or hdr['TRKGAIN'] == 0:
        tit+= 'Closed-Loop\n'
    else:
        tit+= 'Open-Loop\n'
    tit+= 'Algo #%d - Gain=%5.2f\n'%(hdr['TRKALGO'], hdr['TRKGAIN'])
    tit+= 'NAvg %d - %d Frames - %d Bad Fit\n'%(hdr['TRKAVGCT'], hdr['NAXIS2'], vals['ninv'])
    tit+= '%0.3f [ms] Tint - %0.1d FPS\n'%(hdr['TINT']*1000, hdr['FPS'])
    tit+= 'RSTD = %0.2f, CSTD = %0.2f [mas]\n'%(vals['rstd'], vals['cstd'])
    tit+= 'RMean = %0.2f, CMean = %0.2f [mas]'%(vals['rmen'], vals['cmen'])
    plt.suptitle(tit)

    if isSaveFig:
        fig.savefig(vals['flnm']+'.png')
