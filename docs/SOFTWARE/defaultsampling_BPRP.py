import pyvo
import pandas as pd
from zero_point import zpt
import numpy as np
from numpy import array
from astroquery.gaia import GaiaClass
from gaiaxpy import calibrate
import sys

# Connect to Gaia archive
gaia = GaiaClass(gaia_tap_server='https://gea.esac.esa.int/', 
gaia_data_server='https://gea.esac.esa.int/')
#gaia.login()


datacsv = pd.read_csv('SPSS-result.csv')
data = []
frames = ()
id3=datacsv['source_id']
nlen=len(id3)

for  i in range(0, nlen):
    namedr3 = str(id3[i])
    if namedr3 != 'none' : 
     
      result = gaia.load_data(ids=namedr3, format='csv', data_release='Gaia DR3', data_structure='raw', retrieval_type='XP_CONTINUOUS',avoid_datatype_check=True)
      dict = len(result)
      if dict == 0:
         print ("Dictionary is empty")
      else:
        sampling = np.geomspace(330,1049.9999999999, 361)
        continuous_key = [key for key in result.keys() if 'continuous' in key.lower()][0]
        data = result[continuous_key][0]
        data2=data.to_pandas()
        calibrated_spectra, sampling = calibrate(data2, sampling=sampling)

        wave=sampling.flatten().tolist()  #from nparray to list
        flux=calibrated_spectra['flux'].values #now nparray with ndim=1
        flux=np.hstack(flux).tolist()#from nparray 
        fluxe=calibrated_spectra['flux_error'].values #now nparray with ndim=1
        fluxe=np.hstack(fluxe).tolist()#from nparray 
       
        lala=[wave,flux,fluxe]
        lala=np.asarray(lala) #make a numpy array from a list
        columns = ['wave','flux','flux_error']
        df = pd.DataFrame(lala.T,columns=columns)
        df.to_csv(namedr3+"A.csv", encoding="utf-8")
        
        #sys.exit()    

