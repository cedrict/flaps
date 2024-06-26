###############################################################################
#eta_model=0 # constant
#eta_model=1 # steinberger
#eta_model=2 # samuelA
#eta_model=3 # samuelB
###############################################################################

from numba import jit
import numpy as np

#@jit(nopython=True)
def viscosity_Mars(x,z,R1,R2,eta_m,eta_model):
       
    depth=R2-np.sqrt(x**2+z**2)
  
    #--------------------------------------
    if eta_model==0:
       val=eta_m

    #--------------------------------------
    elif eta_model==1: 
       exit('eta_model 1 in viscosity_Mars not implemented')

    #--------------------------------------
    elif eta_model==2: 
       exit('eta_model 2 in viscosity_Mars not implemented')

    #--------------------------------------
    elif eta_model==3: 
       exit('eta_model 3 in viscosity_Mars not implemented')

    #--------------------------------------
    else:
       exit('unknown eta_model')

    return val       

###############################################################################
