# Hide warnings
import warnings
warnings.catch_warnings()
warnings.simplefilter("ignore")
#%matplotlib inline
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator
from pylab import *
import numpy as np
import datetime

import pandas as pd
import pathlib
from matplotlib import rcParams


from astropy.table import Table
data = Table.read('synth_SPSS.tab', format='ascii')
data2=(np.array(data))
df = pd.DataFrame(data2)
df.columns = ['Gsyn','Gphotdr3','BPsyn','BPphotdr3','RPsyn','RPphotdr3','COMMENT','flag','COMMENT2']
dfsel=df[ (df['COMMENT'] == 'OK') & (df['flag']==1) ]
dfsel['deltaG']=dfsel['Gsyn']-dfsel['Gphotdr3']
dfsel['deltaBP']=dfsel['BPsyn']-dfsel['BPphotdr3']
dfsel['deltaRP']=dfsel['RPsyn']-dfsel['RPphotdr3']
dfsel['colBPRP']=dfsel['BPphotdr3']-dfsel['RPphotdr3']
print(dfsel)

xx=dfsel['colBPRP'].values
x=list(flatten(xx))
dG=dfsel['deltaG'].values
dBP=dfsel['deltaBP'].values
dRP=dfsel['deltaRP'].values
Gmag=dfsel['Gphotdr3']
BPmag=dfsel['BPphotdr3']
RPmag=dfsel['RPphotdr3']

xywidth=1.1
linesize=1.1
xytitlesize=8
xylenght1=3
xylenght2=2

rc('axes', linewidth=xywidth)
mpl.rcParams['axes.titlesize'] = xytitlesize #?? no effect
mpl.rcParams['axes.labelsize'] = xytitlesize # this are x- and ytitle
mpl.rcParams['axes.labelpad'] = 1 # distance from the axis
mpl.rcParams['figure.labelweight'] = 9 # ????

mpl.rcParams['xtick.major.size'] = xylenght1
mpl.rcParams['xtick.major.width'] =xywidth
mpl.rcParams['xtick.minor.size'] = xylenght2
mpl.rcParams['xtick.minor.width'] = xywidth
mpl.rcParams['ytick.major.size'] = xylenght1
mpl.rcParams['ytick.major.width'] =xywidth
mpl.rcParams['ytick.minor.size'] = xylenght2
mpl.rcParams['ytick.minor.width'] = xywidth
mpl.rcParams['xtick.labelsize'] = xytitlesize  #this are the number labels
mpl.rcParams['ytick.labelsize'] = xytitlesize

fontsize = 16
mpl.rcParams['lines.linewidth'] = linesize
mpl.rcParams['lines.markersize'] = xywidth
mpl.rcParams['grid.color'] = '#E6E6E6'   #light gray
markersize=7

rect_3 = [ 0.15, .10, 0.72, 0.30]
rect_2 = [ 0.15, .40, 0.72, 0.30]
rect_1 = [ 0.15, .70, 0.72, 0.30]
fig = plt.figure(figsize=(4, 4))

box3 = fig.add_subplot(rect_3)
plt.ylabel('\u0394 $G_{RP}$ [mag]')
plt.xlabel('$G_{ BP}-G_{ RP}$ [mag]')
box3.set_xlim([-0.7, 2.5])
box3.set_ylim([-0.05, +0.05])
box3.plot([-0.7, 3.0],[0.,0.], color ="#C0C0C0")
box3.plot([-0.7, 3.0],[0.01,0.01],linestyle='dashed', color ="#C0C0C0")
box3.plot([-0.7, 3.0],[0.03,0.03],linestyle='dashed', color ="#C0C0C0")
box3.plot([-0.7, 3.0],[-0.01,-0.01],linestyle='dashed', color ="#C0C0C0")
box3.plot([-0.7, 3.0],[-0.03,-0.03],linestyle='dashed', color ="#C0C0C0")
box3.tick_params(axis='both', which='major',direction='in', left=True, top=True, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=True)
box3.tick_params(axis='both', which='minor', direction='in',left=True, top=False, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=True)
box3.minorticks_on
box3.xaxis.set_minor_locator(MultipleLocator(1))
box3.yaxis.set_minor_locator(AutoMinorLocator())
box3.xaxis.set_major_locator(MultipleLocator(0.5))
#box3.xaxis.set_major_formatter('{x:.0f}')
box3.scatter(x=x, y=dRP, c=RPmag, cmap="inferno",s=markersize) 
im3=box3.scatter(x=x, y=dRP, c=RPmag, cmap="inferno",s=markersize) 
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(box3)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im3, cax=cax, orientation='vertical',label='$G_{RP}$ [mag]')
N=dRP.size
mu = dRP.mean()
median = np.median(dRP)
sigma = dRP.std()
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr = '\n'.join((
    r'$\mu=%.3f$' % (mu, ),
    r'$\mathrm{median}=%.3f$' % (median, ),
    r'$\sigma=%.3f$' % (sigma, ),
    r'$N=%.i$' % (N, )))
box3.text(0.05, 0.95, textstr, transform=box3.transAxes, fontsize=7,
        verticalalignment='top')




box2 = fig.add_subplot(rect_2)
plt.ylabel('\u0394 $G_{BP}$ [mag]')
box2.set_xlim([-0.7, 2.5])
box2.set_ylim([-0.05, +0.05])
box2.plot([-0.7, 3.0],[0.,0.], color ="#C0C0C0")
box2.plot([-0.7, 3.0],[0.01,0.01],linestyle='dashed', color ="#C0C0C0")
box2.plot([-0.7, 3.0],[0.03,0.03],linestyle='dashed', color ="#C0C0C0")
box2.plot([-0.7, 3.0],[-0.01,-0.01],linestyle='dashed', color ="#C0C0C0")
box2.plot([-0.7, 3.0],[-0.03,-0.03],linestyle='dashed', color ="#C0C0C0")
box2.tick_params(axis='both', which='major',direction='in', left=True, top=True, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=False)
box2.tick_params(axis='both', which='minor', direction='in',left=True, top=False, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=False)
box2.minorticks_on
box2.xaxis.set_minor_locator(MultipleLocator(1))
box2.yaxis.set_minor_locator(AutoMinorLocator())
box2.xaxis.set_major_locator(MultipleLocator(0.5))
box2.scatter(x=x, y=dBP, c=BPmag, cmap="Blues",s=markersize) 
im2=box2.scatter(x=x, y=dBP, c=BPmag, cmap="Blues",s=markersize) 
divider = make_axes_locatable(box2)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im2, cax=cax, orientation='vertical',label='$G_{BP}$ [mag]')
N=dBP.size
mu = dBP.mean()
median = np.median(dBP)
sigma = dBP.std()
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr = '\n'.join((
    r'$\mu=%.3f$' % (mu, ),
    r'$\mathrm{median}=%.3f$' % (median, ),
    r'$\sigma=%.3f$' % (sigma, ),
    r'$N=%.i$' % (N, )))
box2.text(0.05, 0.95, textstr, transform=box2.transAxes, fontsize=7,
        verticalalignment='top')


box1 = fig.add_subplot(rect_1)
plt.ylabel('\u0394 $G$ [mag]')
box1.set_xlim([-0.7, 2.5])
box1.set_ylim([-0.05, +0.05])
box1.plot([-0.7, 3.0],[0.,0.], color ="#C0C0C0")
box1.plot([-0.7, 3.0],[0.01,0.01],linestyle='dashed', color ="#C0C0C0")
box1.plot([-0.7, 3.0],[0.03,0.03],linestyle='dashed', color ="#C0C0C0")
box1.plot([-0.7, 3.0],[-0.01,-0.01],linestyle='dashed', color ="#C0C0C0")
box1.plot([-0.7, 3.0],[-0.03,-0.03],linestyle='dashed', color ="#C0C0C0")
box1.tick_params(axis='both', which='major',direction='in', left=True, top=True, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=False)
box1.tick_params(axis='both', which='minor', direction='in',left=True, top=False, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=False)
box1.minorticks_on
box1.xaxis.set_major_locator(MultipleLocator(0.5))
#box1.xaxis.set_major_formatter('{x:.0f}')
box1.xaxis.set_minor_locator(MultipleLocator(1))
box1.yaxis.set_minor_locator(AutoMinorLocator())
im1=box1.scatter(x=x, y=dG, c=Gmag, cmap="summer",s=markersize) 
divider = make_axes_locatable(box1)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im1, cax=cax, orientation='vertical',label='$G$ [mag]')
N=dG.size
mu = dG.mean()
median = np.median(dG)
sigma = dG.std()
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr = '\n'.join((
    r'$\mu=%.3f$' % (mu, ),
    r'$\mathrm{median}=%.3f$' % (median, ),
    r'$\sigma=%.3f$' % (sigma, ),
    r'$N=%.i$' % (N, )))
box1.text(0.05, 0.95, textstr, transform=box1.transAxes, fontsize=7,
        verticalalignment='top')
#box1.text(0.05, 0.95, textstr, transform=box1.transAxes, fontsize=14,
#        verticalalignment='top', bbox=props)
plt.savefig('plotsynth_SPSS.eps',dpi=1000)

   

    
