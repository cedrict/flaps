###############################################################################

from numba import jit
import numpy as np

#@jit(nopython=True)
def viscosity_4DEarthBenchmark(x,z,R1,R2,eta_m,etablob,zblob,Rblob):
    val=eta_m

    if np.sqrt(x**2+(z-zblob)**2)<Rblob:
       val=etablob
  
    return val       

###############################################################################
