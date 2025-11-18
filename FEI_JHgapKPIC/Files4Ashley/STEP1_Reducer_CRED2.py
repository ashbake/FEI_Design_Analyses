##########################################
# Script to process CRED2 image cubes and determine the PSF position
# in each each frame of the cube. This is to reduce the data once so that
# we can then plot at will without having to re-run reduction every time.
#
# The PSF finding is done with the JW1DG algo which has proven most accurate
# in multiple simulations
#
# The output of this script can be plotted with the CRED2_Jitter_Plotter.py
# script.
#
#
# NOTE: Updated in January 2023. This is the definitive version as of 01/19/2023
# - Now has option to save figures instead of show during runtime
# - Flux calculation is done using encircled flux and being careful about 
#   consistency in encircled area so that flux values are comparable between frames
# - Minor bug fix in indexing of INVALID_PSF for exception handling
# - Summary values reported at end now exclude failed fit frames
##########################################
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import time
import os,sys

#-- KPIC Imports
sys.path.append('./Files4Ashley/')
import PSF_Algos

startTime = time.time()

#-- Filenames and stuff
# Cred2 Stuff (modify paths as needed)
fl_creddir = f'data_20240417/'
fl_bkgddir = fl_creddir + 'biases/'
fl_credpre = 'cred2_0'

#--File indices (fits file numbers to analyze)
#fl_cred2an = np.array([8,10,14])      # Should be a list (or numpy array) For single frame, just give a single element: [1]
fl_cred2an = np.array([3,5,6,7])      # Should be a list (or numpy array) For single frame, just give a single element: [1]


specify_bkgd = False   # Flag to denote if specific background should be used. If so, provide bkgd file below
#(OPTIONAL) Pull a more recent BG file
#fl_bkgd = fl_bkgddir+'bias_154047_0400.00_00.00249_01_-40.00000_000_000_000_000.fits'

#-- Option to save results or not
isSave = True
isSavePlot = True   # Option to save plot rather than show it during reduction
savePath = fl_creddir+'figures/'
if isSavePlot and not os.path.exists(savePath):
    # Create save directory if it doesn't exist
    os.makedirs(savePath)


#-- Analysis constants
fluxwin = 5.5            # window size for photometric calc
plate_scale = 7.240      # mas/pix
lamOD = 2.2/(10e6) * 206265*1000            # Lambda over D of Keck [mas]


####### Main code ######
#-- Define a function to load data
def loader(fl_num, pth=fl_creddir, pre=fl_credpre):
    flnm = pth+pre+f'{fl_num:03d}'
    with fits.open(flnm+'.fits') as f:
        data = f[0].data
        head = f[0].header
    return (data, head, flnm)

if not isSavePlot:
    plt.ion()

#-- Iterate through data doing analysis
for ind, fl_cred in enumerate(fl_cred2an):
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(14, 4.4))
    print('Frame %d of %d'%(ind+1, len(fl_cred2an)))
    # Load the data frame
    data, hdr, flnm = loader(fl_cred)
    if not specify_bkgd:
        # Find the corresponding background
        fl_bkgd = os.path.basename(hdr['BIAS'])
        fl_bkgd = fl_bkgddir + fl_bkgd
    # Pull the background
    with fits.open(fl_bkgd) as f:
        bkgd = f[0].data
    
    # Subtract background
    data = data - bkgd
    # Remove row with tag pixels
    data = data[:,1:,:]

    # Perform fit on the first frame to find the rough PSF position
    im = data[0]
    _, _, _, cent_x, cent_y, fwhm_x, fwhm_y = PSF_Algos.get_psf_params(im.copy(), 5)
    # Rename coordinate system since x,y always messes me up... who works in XY when you can work in Row,Col?!?!
    colcen, rowcen = cent_x, cent_y
    fwhm_c, fwhm_r = fwhm_x, fwhm_y
    # Display results
    # Full image for reference
    ax1.imshow(im)
    ax1.set_title("Original")
    # Image cropped around found PSF
    ax2.imshow(im)
    # Add crosshair on found-PSF
    ax2.plot(colcen, rowcen, 'r+')  
    ax2.set_title("Fit - Coords (%0.2f, %0.2f) - Size (%0.2f, %0.2f)"%(colcen,rowcen,fwhm_c,fwhm_r))
    # Actually crop the image display now
    cropwin = 10
    ax2.set_xlim([colcen-cropwin, colcen+cropwin+1]); ax2.set_ylim([rowcen-cropwin,rowcen+cropwin+1])
   
    # Preallocate matrix for final results
        # Rows = frames in data, 
        # Col1 = rowCenter, Col2 = colCenter, Col3 = rowSigma, Col4 = colSigma, Col5 = encircled Flux
    Coords = np.full((data.shape[0], 5), np.nan) 
    # Get mask size fora full, well-centered mask
    masksize = PSF_Algos.gen_circ_mask(data[0].copy(), center=np.array(data[0].shape[::-1])//2, radius=fluxwin).sum()
    # Iterate through frames finding PSF and save coords
    for find in range(data.shape[0]):
        if (find%(data.shape[0]//10)) == 0:
            print('    Analyzing %d of %d'%(find+1, data.shape[0]))
        # Perform fit
        try:
            _, _, _, cent_x, cent_y, fwhm_x, fwhm_y = PSF_Algos.get_psf_params(data[find,:,:].copy(), 5)
            # swap coords again cuz no like x,y 
            colcen, rowcen = cent_x, cent_y     
            fwhm_c, fwhm_r = fwhm_x, fwhm_y
            # get encircled flux 
                # Use integer center coords so that mask size is constant between frames
            mask = PSF_Algos.gen_circ_mask(data[find].copy(), center=(int(colcen), int(rowcen)), radius=fluxwin)
            flx = (mask*data[find]).sum() 
            # if PSF at edge, don't compute flux since sum skewed by missing pixels w.r.t frames where PSF fully in frame
                # (make sure mask is full size for all frames)
            if mask.sum() != masksize:
                flx = -9999.
            # (old square calculation method below for reference)
            # flx = data[find,int(rowcen)-fluxwin:int(rowcen)+fluxwin, int(colcen)-fluxwin:int(colcen)+fluxwin].sum()
        except KeyboardInterrupt as e:
            print('keyboard interrupt detected')
            raise e 
        except Exception as e:
            print('--- Fit failed on %d:%d'%(ind,find))
            # Set the values to be obviosly wrong (using default "wrong" value from PSF fitting lib)
            res = PSF_Algos.INVALID_PSF
            # swap coords again cuz no like x,y 
            colcen, rowcen = res[3], res[4]-1   # -1 to account for tag pixel cropping...     
            fwhm_c, fwhm_r = res[5], res[6]  # (old versions had 6,7 rather than 5,6. That's a bug)
            flx = -9999.
        # Save coords (rowcen + 1 to account for cropped tag row)
        Coords[find, :] = np.array([rowcen+1, colcen, fwhm_r, fwhm_c, flx]).T

    # Make a copy of the coords result to be able to save "raw" coords later
    Coords_orig = Coords.copy()
    # Mask out bad coords entries (from where fit failed)
    badentries = np.where(Coords[:,1] == PSF_Algos.INVALID_PSF[3])
    Coords[badentries,:2] = np.nan

    # Display results summary
    #print(Coords[:,0], Coords[:,1])
    std_vals = (np.nanstd(Coords[:,0]*plate_scale), np.nanstd(Coords[:,1]*plate_scale))
    print('-- STDs: row=%f [mas]    col=%f [mas]'%std_vals)

    #-- Plot 2D displacement
    markerSize = 8
    ax3.scatter(Coords[:,1]-np.nanmean(Coords[:,1]), Coords[:,0]-np.nanmean(Coords[:,0]),markerSize,c=np.linspace(0,1,len(Coords[:,1])))
    ax3.set_xlabel('Displacement [pix]')
    ax3.set_ylabel('Displacement [pix]')
    ax3.set_title('PSF Positions (relative to mean)')
    #- Axis limits (fixed axis size)
    axlim = 5
    ax3.set_xlim([-axlim, axlim]); ax3.set_ylim([-axlim,axlim])
    ax3.set_aspect('equal')
    
    #-- Add title to main figure panel marking which frame this is and some analysis properties
    tit = 'Frame #%d'%fl_cred
    tit += '\nSTDC = {:0.2f}, STDR = {:0.2f} [pix]'.format(std_vals[0],std_vals[1])
    tit += '\nPSF Finding failed on %d frames'%badentries[0].shape
    plt.suptitle(tit)   

    plt.tight_layout()
    
    executionTime = (time.time() - startTime)
    print('Total Runtime so far: %f'%executionTime)
    
    if isSave: 
        #-- Save results into new fits file
        fitsfl = fits.PrimaryHDU(data=Coords_orig, header=hdr)
        fitsfl.writeto(flnm+'_proc.fits')

    if isSavePlot:
        #-- Save plot 
        plt.savefig(savePath+os.path.basename(flnm).split('.')[0]+'_step1plot.png')
        plt.close()
        
