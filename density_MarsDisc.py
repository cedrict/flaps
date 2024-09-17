###############################################################################

from numba import jit
import numpy as np

#-------------------R2
# rho_crust
#-------------------R_c_l
# rho_lithosphere
#-------------------R_l_um
# rho_uppermantle
#-------------------R_um_lm
# rho_lowermantle
#-------------------R1

#@jit(nopython=True)
def density_MarsDisc(x,z,R1,R2,rhodisc,R1disc,R2disc,thetadisc,\
                     rho_c,rho_l,rho_um,rho_lm,R_c_l,R_l_um,R_um_lm):

    r=np.sqrt(x**2+z**2)
    theta=np.pi/2-np.arctan2(z,x)

    if r>R_c_l:
       val=rho_c
    elif r>R_l_um:
       val=rho_l
    elif r>R_um_lm:
       val=rho_um
    else:
       val=rho_lm

    if r>R1disc and r<R2disc and theta<thetadisc:
       val=rhodisc

    #val=3550
  
    return val       

###############################################################################
