###############################################################################

from numba import jit
import numpy as np

#-------------------R2
# eta_c
#-------------------R_c_l
# eta_l
#-------------------R_l_um
# eta_um
#-------------------R_um_lm
# eta_lm
#-------------------R1

#@jit(nopython=True)
def viscosity_BartDisc(x,z,R1,R2,etadisc,R1disc,R2disc,thetadisc,\
                       eta_c,eta_l,eta_um,eta_lm,R_c_l,R_l_um,R_um_lm):
    r=np.sqrt(x**2+z**2)
    theta=np.pi/2-np.arctan2(z,x)

    if r>R_c_l:
       val=eta_c
    elif r>R_l_um:
       val=eta_l
    elif r>R_um_lm:
       val=eta_um
    else:
       val=eta_lm

    if r>R1disc and r<R2disc and theta<thetadisc:
       val=etadisc
  
    return val       

###############################################################################
