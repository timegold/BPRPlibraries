# THIS MACRO READS the file 'BP_SPSS-v046a-invWC-019_TEST.csv'
# and, for every source, it creates two pdf files
# BP_5006921282807193856_Apy.pdf
# BP_5006921282807193856_Bpy.pdf

# 1) change the file name (same format)
# 2) run: python  calibrator_BP.py
# 3) to run it, you must have installed  matplotlib, pylab, numpy, pandas, pathlib 
# 4) if you miss any package, install it: e.g., "pip install pylab"
# WRITTEN by Maria Messineo on July,18th,2023

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

############################
def show_spectrum(name,xx,fObs,eObs,fExp,ebar,func1,file): 

    xywidth=1.1
    linesize=1.1

    xytitlesize=8
    xylenght1=3
    xylenght2=2

    rc('axes', linewidth=xywidth)
    mpl.rcParams['axes.titlesize'] = xytitlesize #?? no effect
    mpl.rcParams['axes.labelsize'] = xytitlesize # this are x- and ytitle
    mpl.rcParams['axes.labelpad'] = 2 # distance from the axis
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

    myymax=max(fExp)
    if (myymax > 1)  : 
       lala=  max((fExp+0.1)*1.05)*100.
    if (myymax > .1 and myymax <= 1 )  : 
       lala=  max((fExp+0.01)*1.05)*100.
    if (myymax > .01 and myymax <= .1)  : 
       lala=  max((fExp+0.001)*1.05)*1000.
    if (myymax > .001 and myymax <= .01)  : 
       lala=  max((fExp+0.001)*1.05)*1000.
    lala=round(lala)
    if (myymax > 1)  : 
       lala=lala/100.
    if (myymax > .1 and myymax <= 1 )  : 
       lala=lala/100.
    if (myymax > .01 and myymax <= .1)  : 
       lala=lala/1000.
    if (myymax > .001 and myymax <= .01)  : 
       lala=lala/1000.
    myymax=lala
    

    rect_2 = [ 0.15, .08, 0.82, 0.21]
    rect_1 = [ 0.15, .29, 0.82, .63]

    fig = plt.figure(figsize=(4, 4))
    box2 = fig.add_subplot(rect_2) #plt.axes(rect_2)
    #plt.ylabel('Residual Flux')
    plt.xlabel('pseudo-wavelength')
    box2.set_xlim([0, 59])
    box2.set_ylim([-0.1, +0.1])
    box2.set_ylim(bottom=-0.09900,top=+0.09900)
    #box2.set_yticks(np.arange(-0.1, 0.1, 0.025))
    box2.xaxis.set_major_locator(MultipleLocator(5))
    box2.xaxis.set_major_formatter('{x:.0f}')
    box2.xaxis.set_minor_locator(MultipleLocator(1))
    plt.xticks(rotation=0)
    box2.errorbar(xx, func1, yerr=ebar, fmt='',ecolor='#C0C0C0',color="#000000")
    box2.plot([0.,60.],[0.,0.], color ="#000000")
    box2.tick_params(axis='both',direction='in',left=True, top=False, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=True)
    box2.tick_params(axis='both', which='minor', direction='in',left=True, top=False, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=False)
    box2.minorticks_on
    box2.yaxis.set_minor_locator(AutoMinorLocator())
    plt.grid()

    box1 = fig.add_subplot(rect_1) #plt.axes(rect_1)
    #plt.ylabel('Flux')
    box1.set_xlim([0, 59])
    box1.set_ylim(ymin=0,ymax=myymax)
    box1.plot(xx, fExp,color ="#000000",)
    box1.plot(xx, fObs,color ="#6495ED") #corn
    box1.tick_params(axis='both', which='major',direction='in', left=True, top=False, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=False)
    box1.tick_params(axis='both', which='minor', direction='in',left=True, top=False, right=True, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=False)
    box1.minorticks_on
    box1.xaxis.set_major_locator(MultipleLocator(5))
    box1.xaxis.set_major_formatter('{x:.0f}')
    box1.xaxis.set_minor_locator(MultipleLocator(1))
    box1.yaxis.set_minor_locator(AutoMinorLocator())
    plt.grid()
    
    #second axis
    ax2 = box1.twiny()
    xtickname=['1100','800','600','500','400','300']
    pollo=[2.01477,9.56506,17.1320,23.4787,35.4242,59.9221]
    ax2.set_xlim(box1.get_xlim())
    ax2.spines['top'].set_color("#6495ED") 
    ax2.set_xticklabels(xtickname,color ="#6495ED")
    ax2.set_xlabel(r"wavelength[nm]",color ="#6495ED")
    ax2.set_xticks(pollo)
    ax2.tick_params(axis='both',  which='major',direction='in', color ="#6495ED",left=False, top=True, right=False, bottom=False, labelleft=False, labeltop=True, labelright=False, labelbottom=False)
    plt.savefig(file,dpi=1000)

##############

pdir = pathlib.Path('.').resolve()
ff1 = 'BP_SPSS-v046a-invWC-019_TEST.csv'
df1 = pd.read_csv(ff1)

###insert loop
for index, row in df1.iterrows():
     dfsel = df1[df1['id'] == row['id'] ]
     nameloop=dfsel['id'].iloc[0]
     if (index == 0 or nameloop != name) :
         name=nameloop
         print(f'Full dataset: {len(df1)}, subset: {len(dfsel)}')

         dfsel['delta']=dfsel['fObs']-dfsel['fExp']
         xx=dfsel[dfsel.columns[1]]
         fObs=dfsel[dfsel.columns[3]]/10000.
         eObs=dfsel[dfsel.columns[4]]/10000.
         fExp=dfsel[dfsel.columns[5]]/10000.
         ebar=eObs/fExp
         func1=(fObs-fExp)/fExp
         aa=((dfsel['id'].iloc[0]))
         fileA="BP_"+str(aa)+"_Apy.pdf"
         show_spectrum(name,xx,fObs,eObs,fExp,ebar,func1,fileA)

         aa=dfsel.loc[dfsel['u'].between(13,50 ),['fObs']].values
         bb=dfsel.loc[dfsel['u'].between(13,50 ),['fExp']].values
         fact=median(aa/bb)
         func2=(fObs/fact-fExp)/fExp
         fObs2=fObs/fact
         aa=((dfsel['id'].iloc[0]))
         fileB="BP_"+str(aa)+"_Bpy.pdf"
         show_spectrum(name,xx,fObs2,eObs,fExp,ebar,func2,fileB)
