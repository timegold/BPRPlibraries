#this script
#is adapted from
#https://docs.datacentral.org.au/help-center/virtual-observatory-examples/ssa-galah-dr3-interactive-spectra-explorer-enhanced-data-central-api/
#adapted by Maria Messineo, 18 July 2023

#Given an ascii file with the GALAH  sobject_id located in column XX (14 in my example)
#it connects to the server load the normalized spectra and saves it in the local disk 
# e.g. figure 141202003001147.eps

#if you uncomment the line
#     #write_spectrum(namedr3,[axB,axV,axR,axI],namedr3)
#then a copy of the 4 fits is also saved in the local disk

from specutils import Spectrum1D
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import matplotlib.ticker as plticker
from astropy.stats import sigma_clip
from astropy.time import Time
import requests
import re
from io import BytesIO
from astropy import units as u
from pyvo.dal.ssa  import search, SSAService

import pyvo
import requests
import wget
import sys
#called by onpick to retrieve and show the spectra of the target of interest


def show_spectrum(name,axes,title): 
    #query the SSA service
    #no position required as we already know the target name from the DataCentral API query
    url = "https://datacentral.org.au/vo/ssa/query"
    service = SSAService(url)
    custom = {}
    custom['TARGETNAME'] = name
    custom['COLLECTION'] = 'galah_dr3'
    results = service.search(**custom)
    df = results.votable.get_first_table().to_table(use_names_over_ids=True).to_pandas()
    filters = ['B','V','R','I']
    colours = ["#085dea","#1e8c09","#cf0000","#640418"]

    #go through each filter and plot its spectrum (if available)
    for idx in range(0,4):
        filt = filters[idx]
        ax = axes[idx]
        #remove any previous spectrum/labels/titles that may have been plotted in a previous call of 
        #the show_spectrum function
        ax.clear()
        #show the title in the first position only
        if(idx == 0):
            ax.set_title(title)
        #select only the spectrum of the filter of interest
        subset = df[(df['band_name'] == filt)].reset_index()
        #give preference to using the continuum normalised spectra 
        print(subset['dataproduct_subtype'])

#0    normalised
#1      combined
#Name: dataproduct_subtype, dtype: object
        
        if(subset[subset['dataproduct_subtype'].isin(['normalised'])].shape[0] > 0):
            subset = subset[(subset['dataproduct_subtype'] == "normalised")].reset_index()
        #only proceed if we have the filter of interest available in the results 
        if(subset.shape[0] > 0):
            #add RESPONSEFORMAT=fits here to ensure we get fits format back
            url= subset.loc[0,'access_url'] + "&RESPONSEFORMAT=fits"
            #load the spectrum
#SPECUTILS  https://specutils.readthedocs.io/en/stable/api/specutils.Spectrum1D.html#
            spec = Spectrum1D.read(BytesIO(requests.get(url).content),format='wcs1d-fits')
            spec.write(name+filt+".fits", format="wcs1d-fits",overwrite=True)
            exptime = subset.loc[0,'t_exptime']
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            loc = plticker.MultipleLocator(base=20.0)
            ax.xaxis.set_major_locator(loc)
            #plot label at last position (caveat: no check if spectra are missing)
            if(idx == 3): 
                ax.set_xlabel("Wavelength ($\mathrm{\AA}$)",labelpad=10)
            #plot the spectrum 
            ax.plot(spec.wavelength,spec.flux,linewidth=LWIDTH,color=colours[idx])
            #adjust the y-scale to best fit the spectrum with some sigma clipping
            nspec = len(spec.spectral_axis)
            clipped = sigma_clip(spec[int(0+0.01*nspec):int(nspec-0.01*nspec)].flux,masked=False,sigma=15)
            ymin = min(clipped).value
            ymax = max(clipped).value 
            xmin = spec.wavelength[0].value
            xmax = spec.wavelength[nspec-1].value
            #add a 1% buffer either side of the x-range
            dx=0.01*(xmax-xmin)
            ax.set_xlim(xmin-dx,xmax+dx)
            #add a 1% buffer either side of the y-range
            dy=0.03*(ymax-ymin)
            ax.set_ylim(ymin-dy,ymax+dy)
        #else:
            #print missing data...
            #print spectrum on local disk
            
            
def write_spectrum(name,axes,title): 
    #query the SSA service
    #no position required as we already know the target name from the DataCentral API query
    url = "https://datacentral.org.au/vo/ssa/query"
    service = SSAService(url)
    custom = {}
    custom['TARGETNAME'] = name
    custom['COLLECTION'] = 'galah_dr3'
    results = service.search(**custom)
    df = results.votable.get_first_table().to_table(use_names_over_ids=True).to_pandas()
    filters = ['B','V','R','I']
    colours = ["#085dea","#1e8c09","#cf0000","#640418"]

    #go through each filter and plot its spectrum (if available)
    for idx in range(0,4):
        filt = filters[idx]
        ax = axes[idx]
        #remove any previous spectrum/labels/titles that may have been plotted in a previous call of 
        #the show_spectrum function
        ax.clear()
        #show the title in the first position only
        if(idx == 0):
            ax.set_title(title)
        #select only the spectrum of the filter of interest
        subset = df[(df['band_name'] == filt)].reset_index()
        #give preference to using the continuum normalised spectra 
        if(subset[subset['dataproduct_subtype'].isin(['normalised'])].shape[0] > 0):
            subset = subset[(subset['dataproduct_subtype'] == "normalised")].reset_index()
        #only proceed if we have the filter of interest available in the results 
        if(subset.shape[0] > 0):
            #add RESPONSEFORMAT=fits here to ensure we get fits format back
            url= subset.loc[0,'access_url'] + "&RESPONSEFORMAT=fits"
            #load the spectrum
#SPECUTILS  https://specutils.readthedocs.io/en/stable/api/specutils.Spectrum1D.html#
            spec = Spectrum1D.read(BytesIO(requests.get(url).content),format='wcs1d-fits')
            spec.write(name+filt+".p.fits", format="wcs1d-fits",overwrite=True)


#*********************************************************************************************
###MARIA"S MODIFICATION
filemy='/Users/mmessine/DUTY/Spectral_Libraries/MASTERLIST/master_galah.tab'
#https://stackoverflow.com/questions/26000336/execute-curl-command-within-a-python-script


#set the figure size to be wide in the horizontal direction
#to accommodate the image cutouts alongside the spectra 
fsize=[8,8]
mpl.rcParams['axes.linewidth'] = 0.7
FSIZE=18
LWIDTH=0.5
LABELSIZE=10

gs = gridspec.GridSpec
fig = plt.figure(figsize=fsize)
gs = gridspec.GridSpec(4,1)
#B
axB = fig.add_subplot(gs[0,0])
#V
axV = fig.add_subplot(gs[1,0])
#R
axR = fig.add_subplot(gs[2,0])
#I
axI = fig.add_subplot(gs[3,0])


f = open(filemy, 'r')
data = []
frames = ()
for line in f:
    line = line.strip()
    columns = line.split()
    namedr3 = columns[14] #;;;column nummer 5-1, zero is counted
    print(namedr3)
# Submit queries;;;this retrieves all neigbours too...clean after
#curl -k -b cookies.txt "https://datacentral.org.au/vo/slink/links?ID=141202003001147&DR=galah_dr3&IDX=0&FILT=V&RESPONSEFORMAT=fits" > test.fits

    #write_spectrum(namedr3,[axB,axV,axR,axI],namedr3)
    show_spectrum(namedr3,[axB,axV,axR,axI],namedr3)
    plt.savefig(namedr3+'.eps',dpi=300,bbox_inches="tight")

    #plt.show()
#called when the user clicks a target

    
def onpick(event):
    ind = event.ind
    pt = event.artist.get_offsets()
    #mouse positions
    mx = event.mouseevent.xdata
    my = event.mouseevent.ydata
    nb = input('Choose a number')
    if(nb == 0): chosen = None
    if(nb == 1): destroy()

 