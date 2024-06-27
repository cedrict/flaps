###############################################################################

from numba import jit
import numpy as np

#@jit(nopython=True)
def density_4DEarthBenchmark(x,z,R1,R2,rho_m,rhoblob,zblob,Rblob):
    val=rho_m

    #-------------------------------------
    #elif exp==3:
    #   r=np.sqrt(x*x+y*y)
    #   theta=np.pi/2-np.arctan2(y,x)
    #   if theta<np.pi/8 and r>R1+3*(R2-R1)/8 and r<R1+5*(R2-R1)/8:
    #      val*=rhoblob
    #else:

    if np.sqrt(x**2+(z-zblob)**2)<Rblob:
       val=rhoblob

    return val

###############################################################################
