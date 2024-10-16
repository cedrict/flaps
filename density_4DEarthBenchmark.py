###############################################################################

from numba import jit
import numpy as np

#@jit(nopython=True)

def density_4DEarthBenchmark(x,z,R1,R2,crust_rho,lithosphere_rho,\
                             uppermantle_rho,lowermantle_rho,blob_rho,blob_z,blob_R):
    val=lowermantle_rho

    #-------------------------------------
    #elif exp==3:
    #   r=np.sqrt(x*x+y*y)
    #   theta=np.pi/2-np.arctan2(y,x)
    #   if theta<np.pi/8 and r>R1+3*(R2-R1)/8 and r<R1+5*(R2-R1)/8:
    #      val*=rhoblob
    #else:

    if np.sqrt(x**2+(z-blob_z)**2)<blob_R:
       val=blob_rho

    return val

###############################################################################
