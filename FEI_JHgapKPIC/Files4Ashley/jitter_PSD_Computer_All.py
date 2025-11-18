'''

NOTE:::: NEED TO UPDATE THE DESCRIPTION BELOW, WHICH IS OUT OF DATE

This script loads all the files in the CombineDirs subdir of the 2_VFN_TrackLoopCharac_NoInject dir
It computes the PSD for all of these samples and plots them

It ALSO takes ratios of CL vs. OL PSDs

__NOTES:__
- All saving elements are shielded behind isSave** cases
    -->> therefore it is okay to run this script with isSave** = False and that will not overwrite any saved results

- 'rowcen' and 'colcen' (ie. PSF coordinates) are relative to average position, not the goal
    -->> This means that this analysis is only truly valid for jitter, NOT for drift error
'''
from astropy.io import fits
import numpy as np
from scipy import fft, signal
import matplotlib.pyplot as plt
from glob import glob
import pandas as pd
from copy import deepcopy

plt.ion()
main_fontsize = 18      # Default fontsize for figures
stack_fontsize = 14     # fontsize for figures w/ 2 vertically-stacked subplots
plt.rcParams.update({'font.size': main_fontsize})

#-- Inputs
datapath = 'rawdata/CRED2/'
OLdatapath = datapath     # path to open-loop data
platescale = 7.24   # [mas/pix]

isPlotHist = True
isSaveCSV = False
isPlotOverlayPSD = False
isSaveOverlayPSD = False
isPlotSamplePSD = False  
isSavePSD = False
isDoGain0p3 = False
isSaveGain0p3 = False
isDoGain0p1 = False
isSaveGain0p1 = False
isDoFPS400 = False
isSaveFPS400 = False

flnm_csv = 'summary.csv'
savepath = datapath+'reduced/results/'
psd_figsize = (11,7)
stack_figsize = (11,10)
figdpi = 200

#-- Glob to find all relevant files
flnms = glob(datapath+'*_proc.fits')
flnms = sorted(flnms)

#-- Define function to load data
def dataloader(flnm):
    with fits.open(flnm) as f:
        pztDict = {
            'coords': f[0].data,
            'times': f[1].data,
            'hdr': f[0].header}
        # NOTE: for simplicity "centering" is defined relative to average position NOT goal (see note at top)
        pztDict['rowcen'] = pztDict['coords'][:,0] - pztDict['coords'][:,0].mean()
        pztDict['colcen'] = pztDict['coords'][:,1] - pztDict['coords'][:,1].mean()

    return pztDict

#-- Load data
data = []
for flnm in flnms:
    ddict = dataloader(flnm)
    data.append(ddict)
    ddict['flnm'] = flnm.split('/')[-1].strip('_proc.fits').replace('_',' ')
    ddict['rowstd'] = ddict['rowcen'].std()
    ddict['colstd'] = ddict['colcen'].std()
    ddict['rowmean'] = ddict['rowcen'].mean()
    ddict['colmean'] = ddict['colcen'].mean()
    ddict['nsamp'] = ddict['coords'].shape[0]
    # Create a "simple" title that removes unecessary elements reformats gain element
    ddict['title'] = ddict['flnm'].replace('KPICSrc ', '').replace('MaxTint ', '').replace('0p', '0.') 

#-- Plot histogram
if isPlotHist:
    # Note: briefly change the fontsize so that all text fits in figure
    plt.rcParams.update({'font.size': 12})
    for ddict in data:
        plt.figure()
        _,bins,_ = plt.hist(x=ddict['rowcen'], bins='auto', color='blue', alpha=0.5, rwidth=0.85, label='Row Jitter')
        _,bins,_ = plt.hist(x=ddict['colcen'], bins=bins, color='red', alpha=0.5, rwidth=0.85, label='Col Jitter')
        plt.xlabel('Jitter [pix]')
        plt.ylabel('Occurrences')
        title = ddict['title'] + '\n'
        title += 'RSTD = %0.2f [mas] | CSTD = %0.2f [mas]'%(ddict['rowstd']*platescale, ddict['colstd']*platescale)
        plt.title(title)
        plt.legend()
    # Reset fontsize
    plt.rcParams.update({'font.size': main_fontsize})

#-- Print in csv format
if isSaveCSV:
    with open(savepath+flnm_csv, 'x') as f:
        f.write('filename,rowstd[mas],colstd[mas],rowmean[pix],colmean[pix],nsamp,rowstd[pix],colstd[pix]\n')
        for ddict in data:
            row = [ddict['flnm'].replace(' ','_'),ddict['rowstd']*platescale,ddict['colstd']*platescale,ddict['rowmean'],ddict['colmean'],ddict['nsamp'],ddict['rowstd'],ddict['colstd']]
            row = [str(el) for el in row]
            row = ','.join(row)
            f.write(row+'\n')

#-- Define a function to compute PSD
def PSD_computer(pztDict):
    # Deal with timing elements
    tstamps = pztDict['times'].copy()
    tstamps -= tstamps[0]
    sample_rate = 1/np.mean(np.diff(tstamps))   # [Hz] sampling from average timedelta

    # Get PSD using Welch's method
        # Set welch window as ~1/4 of total datapoints and round to closest lower power of 2
    welch_NperSeg = 1 << int(np.log2(len(tstamps)/4))
    rFreq, rPSD = signal.welch(pztDict['rowcen'], sample_rate, nperseg=welch_NperSeg, scaling='density')
    cFreq, cPSD = signal.welch(pztDict['colcen'], sample_rate, nperseg=welch_NperSeg, scaling='density')

    # Store resuts in dict
    pztDict['srate'] = sample_rate
    pztDict['freq'] = rFreq     # rFreq = cFreq so only use one
    pztDict['rowpsd'] = rPSD
    pztDict['colpsd'] = cPSD

    return pztDict

#-- Compute PSD
for ddict in data:
    ddict = PSD_computer(ddict)

#-- Display PSDs (overlaying all samples on the same plot)
if isPlotOverlayPSD or isSaveOverlayPSD:
    plt.figure(figsize=psd_figsize)
    for ddict in data:
        plt.loglog(ddict['freq'], ddict['rowpsd'], ':', linewidth=1, markersize=2, color='blue')
        plt.loglog(ddict['freq'], ddict['colpsd'], ':', linewidth=1, markersize=2, color='orange')
    plt.legend(['Row', 'Col'])
    plt.title('All 2_ Data')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power [pix$^2$/Hz]')
    plt.ylim(1e-8 ,3e-1)
    plt.xlim(1e-2, 5e3)
    plt.tight_layout()

    if isSaveOverlayPSD:
        plt.savefig(savepath+'AllPSDSamples_Overlaid.png', dpi=figdpi)

    # Display zoom of previous plot
    plt.figure(figsize=psd_figsize)
    for ddict in data:
        plt.loglog(ddict['freq'], ddict['rowpsd'], ':', linewidth=1, markersize=2, color='blue')
        plt.loglog(ddict['freq'], ddict['colpsd'], ':', linewidth=1, markersize=2, color='orange')
    plt.legend(['Row', 'Col'])
    plt.title('All 2_ Data')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power [pix$^2$/Hz]')
    plt.ylim(2e-7, 1e-2)
    plt.xlim(4e0, 2.1e2)
    plt.tight_layout()

    if isSaveOverlayPSD:
        plt.savefig(savepath+'AllPSDSamples_Overlaid_zoom.png', dpi=figdpi)

#-- Average similar samples and display PSDs 
avg_data = []
# NOTE: The logic in this for-loop assumes the samples are in order! (ie. data has been sorted so similar samples are next to each other) 
# NOTE: It also assumes at most 2 samples per dataset! (successfully deals with 400FPS datasets through dedicated filter)
prev_flnm = ''
for sind, ddict in enumerate(data):
    # get filename ommitting sample number (so that it just includes the settings)
    cur_flnm = ' '.join(ddict['flnm'].split(' ')[:-1])
    samp_num = str(ddict['flnm'][-1])

    # check if repeat settings
    same_settings = (prev_flnm == cur_flnm)
    # Prep the prev_flnm for the next itr
    prev_flnm = cur_flnm

    if not same_settings:

        if isPlotSamplePSD or isSavePSD:
            # Unique file settings found so create a new figure
            plt.figure(figsize=psd_figsize)

        # Since unique settings, this must be the start of a new dataset
        set_freqs = []
        set_rowPSD = []
        set_colPSD = []

        # Add list for files that were combined
        set_flnms = []
        # Add inds of files that were combined
        set_inds = []

    if isPlotSamplePSD or isSavePSD:
        # Plot the sample
        plt.loglog(ddict['freq'], ddict['rowpsd'], ':', linewidth=1, markersize=2, color='blue', label='Row - Sample '+samp_num)
        plt.loglog(ddict['freq'], ddict['colpsd'], ':', linewidth=1, markersize=2, color='orange', label='Col - Sample '+samp_num)

    # Add the sample to "set"
    set_freqs.append(ddict['freq'])
    set_rowPSD.append(ddict['rowpsd'])
    set_colPSD.append(ddict['colpsd'])
    set_flnms.append(ddict['flnm'])
    set_inds.append(sind)

    if same_settings or ('400FPS' in ddict['flnm']):
        # This is the second file in the dataset (or is a 400FPS file which have only 1 sample) 
        
        if same_settings:
            # Since second sample in dataset, compute average values (resampling to get clean behavior)
            avg_rowPSD = np.array( [ np.interp(set_freqs[0], set_freqs[1], set_rowPSD[1]) , set_rowPSD[0] ] ).mean(axis=0)
            avg_colPSD = np.array( [ np.interp(set_freqs[0], set_freqs[1], set_colPSD[1]) , set_colPSD[0] ] ).mean(axis=0)
            avg_freq = set_freqs[0]
            if isPlotSamplePSD or isSavePSD:
                # Add them to the plot
                plt.loglog(avg_freq, avg_rowPSD, '-', linewidth=2, color='blue', label='Row - Avg')
                plt.loglog(avg_freq, avg_colPSD, '-', linewidth=2, color='orange', label='Col - Avg')
        else:
            # Set avg_ to the single-sample setting
            avg_rowPSD = set_rowPSD[0]
            avg_colPSD = set_colPSD[0]
            avg_freq = set_freqs[0]

        # Create dictionary of averaged entry
        adict = {
            'cur_flnm'  : cur_flnm,
            'rowpsd'    : avg_rowPSD,
            'colpsd'    : avg_colPSD,
            'freq'      : avg_freq,
            'flnms'     : set_flnms,
            'inds'      : set_inds,
            'title'     : ddict['title'][:-5]
                }
        #- Add key settings to dictionary
        title = adict['title']
        # Gain (Note: assumes all names start with "Gain#.#")
        adict['gain'] = float(adict['title'][4:7])  
        # Noisy Col status
        if 'OnNoisyCol' in title:
            adict['noisy'] = True
        else:
            adict['noisy'] = False
        # Ctrl code type
        if ('CTRLOn' in title) and ('PZTOn' in title):
            adict['mode'] = 'kpicnorm'
        elif ('CTRLOff' in title) and ('PZTOn' in title):
            adict['mode'] = 'MOV'
        elif ('CTRLOff' in title) and ('PZTOff' in title):
            adict['mode'] = 'SVA'
        elif 'OL' in title:
            adict['mode'] = 'OL'
        else:
            raise ValueError('Unrecognized Mode')
        # FPS (Note: assumes maximum 4- and minimum 3-character FPS with space before FPS number
        adict['fps'] = int(title.split('FPS')[0][-4:])

        # Store averaged dict into list of averages
        avg_data.append(adict)
         
        if isPlotSamplePSD or isSavePSD:
            # Now format the figure
            plt.legend()
            title = ddict['title']
            plt.title(title)
            plt.xlabel('Frequency [Hz]')
            plt.ylabel('Power [pix$^2$/Hz]')
            #lolim, _ = plt.xlim()
            #plt.xlim(lolim, uplim)
            plt.ylim(1e-8 ,3e-1)
            plt.xlim(1e-2, 5e3)
            plt.tight_layout()

        print('%d) COMBINED FILES (%d):\n\t'%(len(avg_data),len(set_rowPSD))+set_flnms[0])
        try:
            print('\t'+set_flnms[1])
        except:
            pass

        if isSavePSD:
            plt.savefig(savepath+prev_flnm.replace(' ','_')+'_PSD.png', dpi=figdpi)

#-- Load and average Open-Loop data
#- First do 6000FPS data w/ settings similar to majority of 2_ data
OL_flnms = sorted(glob(OLdatapath+'OffNoisyCol*PZTOn*HEPAOff*_proc.fits')) 
OL_data = []
set_flnms = []
for flnm in OL_flnms:
    # Skip the 0002 sample which had weird time sampling
    if '0002' in flnm:
        continue
    ddict = dataloader(flnm)
    ddict['flnm'] = flnm.split('/')[-1].strip('_proc.fits').replace('_',' ')
    ddict['rowstd'] = ddict['rowcen'].std()
    ddict['colstd'] = ddict['colcen'].std()
    ddict['rowmean'] = ddict['rowcen'].mean()
    ddict['colmean'] = ddict['colcen'].mean()
    ddict['nsamp'] = ddict['coords'].shape[0]
    # Create a "simple" title that removes unecessary elements reformats gain element
    ddict['title'] = ddict['flnm'].replace('KPICSrc ', '').replace('MaxTint ', '').replace('0p', '0.') 

    # Compute PSD
    ddict = PSD_computer(ddict)

    # List filenames that went into this set
    set_flnms.append(ddict['flnm'])

    OL_data.append(ddict)

print('OL 1) COMBINED FILES (%d):\n\t'%(len(set_flnms)) + set_flnms[0])
print('\t' + set_flnms[1])

# Combine the two samples - including resampling (using the first sample as the reference frequency)
avg_rowPSD = np.array([np.interp(OL_data[0]['freq'], OL_data[1]['freq'], OL_data[1]['rowpsd']), OL_data[0]['rowpsd']]).mean(axis=0)
avg_colPSD = np.array([np.interp(OL_data[0]['freq'], OL_data[1]['freq'], OL_data[1]['colpsd']), OL_data[0]['colpsd']]).mean(axis=0)

# Create dictionary of averaged OL entry    (populated with same stuff as other avg_data elements)
adict = {
    'cur_flnm'      : ' '.join(ddict['flnm'].split(' ')[:-1]),
    'rowpsd'        : avg_rowPSD,
    'colpsd'        : avg_colPSD,
    'freq'          : OL_data[0]['freq'],
    'flnms'         : set_flnms,
    'title'         : ddict['title'][:-5],
    'fps'           : int(ddict['title'].split('FPS')[0][-4:])
        }
if 'OnNoisyCol' in adict['title']:
    adict['noisy'] = True
else:
    adict['noisy'] = False
# Add to list of OL samples
OL_avg_data = [adict]

#- Now add the 400 FPS, OffNoisyCol data
# find the entry with the relevant OL data
for ddict in avg_data:
    # loop until data matching the requested settings is found
    if (ddict['mode'] == 'OL') and (not ddict['noisy']) and (ddict['fps'] == 400):
        break
# Clean up dict entries to match other OL_avg_data entry
adict = deepcopy(ddict)
adict.pop('inds')
adict.pop('gain')
adict.pop('mode')
# Add to the list
OL_avg_data.append(adict)


#-- Define function to compare CL to OL operation (and plot/save results)
def OLvCLAnalyzer(fps, gain, isSave):
    #-- Decrease fontsize for this big figure
    plt.rcParams.update({'font.size': stack_fontsize})
    
    #-- Plot samples with 0.3 gain and 6000FPS
    CL_list = []    # list containing the CL samples relevant to this subset
    fig, [axR, axC] = plt.subplots(2,1, sharex = True, sharey = True, figsize=stack_figsize)
    for adict in avg_data:
        if (adict['gain'] != gain) or (adict['fps'] != fps) or (adict['mode'] == 'OL'):
            # This is not a sample matching the requested values, so skip it
            continue
        if adict['noisy']:
            # add suffix to label marking that this is on a noisy col
            labsuf = ' (Noisy Col)'
        else:
            labsuf = ''
        axR.loglog(adict['freq'], adict['rowpsd'], '-', linewidth=2, label=adict['mode'].upper()+labsuf)
        axC.loglog(adict['freq'], adict['colpsd'], '-', linewidth=2, label=adict['mode'].upper()+labsuf)

        # Add sample to the list
        CL_list.append(adict)
    # Overlay open-loop data
      # First find the right OL sample
    aind = [ind for ind, adict in enumerate(OL_avg_data) if adict['fps'] == fps][0]
    OL_dict = OL_avg_data[aind]
    axR.loglog(OL_dict['freq'], OL_dict['rowpsd'], ':', linewidth=2, label='OL Sample')
    axC.loglog(OL_dict['freq'], OL_dict['colpsd'], ':', linewidth=2, label='OL Sample')

    # Format figure
    axR.legend(loc = 'upper right')
    axC.legend(loc = 'upper right')
    axR.set_title('Row Axis')
    axC.set_title('Col Axis')
    fig.suptitle('Tracking Script Assessment\nGain %0.1f FPS %d'%(gain, fps))
    fig.supxlabel('Frequency [Hz]')
    fig.supylabel('Power [pix$^2$/Hz]')
    plt.ylim(1e-8, 3e-1)
    plt.xlim(1e-2, 5e3)
    plt.tight_layout()
    if isSave:
        plt.savefig(savepath+'Gain%0.1f_FPS%d_samples_combined_PSD.png'%(gain,fps), dpi=figdpi)

    #-- Compute and display ratio of CL/OL PSD
    fig, [axR, axC] = plt.subplots(2,1, sharex = True, sharey = True, figsize=stack_figsize)
    for adict in CL_list:
        # resample CL data to OL sampling
        CL_rowPSD = np.interp(OL_dict['freq'], adict['freq'], adict['rowpsd'])
        CL_colPSD = np.interp(OL_dict['freq'], adict['freq'], adict['colpsd'])

        # take ratio
        row_ratio = CL_rowPSD / OL_dict['rowpsd']
        col_ratio = CL_colPSD / OL_dict['colpsd']

        # plot
        if adict['noisy']:
            # add suffix to label marking that this is on a noisy col
            labsuf = ' (Noisy Col)'
        else:
            labsuf = ''
        axR.loglog(OL_dict['freq'], row_ratio, '-', linewidth=2, label=adict['mode'].upper()+labsuf)
        axC.loglog(OL_dict['freq'], col_ratio, '-', linewidth=2, label=adict['mode'].upper()+labsuf)

    # Add horizontal line at y=1 (ie. marking where control loop isn't changing anything)
    axR.axhline(y=1, linestyle='--', color='gray')
    axC.axhline(y=1, linestyle='--', color='gray')
    axR.legend(loc = 'upper left')
    axC.legend(loc = 'upper left')
    axR.set_title('Row Axis')
    axC.set_title('Col Axis')
    fig.suptitle('Tracking Script Assessment\nGain %0.1f FPS %d'%(gain, fps))
    fig.supxlabel('Frequency [Hz]')
    fig.supylabel('PSD Closed-Loop / Open-Loop')
    plt.ylim(1e-4, 7e2)
    plt.xlim(1e-2, 5e3)
    plt.tight_layout()
    if isSave:
        plt.savefig(savepath+'Gain%0.1f_FPS%d_Ratio_PSD.png'%(gain,fps), dpi=figdpi)

    #-- Reset fontsize
    plt.rcParams.update({'font.size': main_fontsize})

#-- Analyze Gain=0.3, FPS=6000 data
if isDoGain0p3 or isSaveGain0p3:
    OLvCLAnalyzer(6000, 0.3, isSaveGain0p3)

#-- Analyze Gain=0.1, FPS=6000 data
if isDoGain0p1 or isSaveGain0p1:
    OLvCLAnalyzer(6000, 0.1, isSaveGain0p3)

#-- Analyze Gain=0.3, FPS=400 data
if isDoFPS400 or isSaveFPS400:
    OLvCLAnalyzer(400, 0.3, isSaveGain0p3)
