import pyvo
import pandas as pd
from zero_point import zpt
import numpy

from astroquery.gaia import GaiaClass

# Connect to Gaia archive
gaia = GaiaClass(gaia_tap_server='https://gea.esac.esa.int/', gaia_data_server='https://gea.esac.esa.int/')
gaia.login()

datacsv = pd.read_csv('../SPSS_DR3_GAIA-BPRP/SPSS-result.csv')
data = []
frames = ()
id3=datacsv['source_id']
nlen=len(id3)

#append them them all in one file
result = gaia.load_data(ids=id3, data_release = 'Gaia DR3', retrieval_type='XP_SAMPLED', data_structure = 'COMBINED', format='csv')
dict = len(result)
if dict == 0:
   print ("Dictionary is empty")
else:
   print ("Dictionary is not empty")
sampled_key = [key for key in result.keys() if 'sampled' in key.lower()][0]
data = result[sampled_key][0]
filename='alltogheter.DTL.csv'
data.write(filename, format='csv',overwrite=True)


#save them them in individual file
for  i in range(0, nlen):
    namedr3 = str(id3[i])

    if namedr3 != 'none' : 
      retrieval_type = 'XP_SAMPLED' # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
      data_structure = 'COMBINED'   # Options are: 'INDIVIDUAL', 'COMBINED', 'RAW'
      data_release   = 'Gaia DR3'     # Options are: 'Gaia DR3' (default), 'Gaia DR2'
      result = gaia.load_data(ids=namedr3, data_release = 'Gaia DR3', retrieval_type='XP_SAMPLED', data_structure = 'INDIVIDUAL', format='csv')
      dict = len(result)
      if dict == 0:
         print ("Dictionary is empty")
      else:
         print ("Dictionary is not empty")
         sampled_key = [key for key in result.keys() if 'sampled' in key.lower()][0]
         data = result[sampled_key][0]
         filename=namedr3+'.DTL.csv'
         data.write(filename, format='csv',overwrite=True)


