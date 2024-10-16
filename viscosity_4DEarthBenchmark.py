###############################################################################

from numba import jit
import numpy as np

#@jit(nopython=True)
def viscosity_4DEarthBenchmark(x,z,R1,R2,crust_eta,lithosphere_eta,\
                               uppermantle_eta,lowermantle_eta,etablob,zblob,Rblob):

    val=lowermantle_eta

    if np.sqrt(x**2+(z-zblob)**2)<Rblob:
       val=etablob
  
    return val       

###############################################################################
