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
def viscosity_MarsDisc(x,z,R1,R2,\
                       blob_eta,blob_z,blob_R,blob_R1,blob_R2,blob_theta,\
                       crust_eta,crust_depth,\
                       lithosphere_eta,lithosphere_depth,\
                       uppermantle_eta,uppermantle_depth,\
                       lowermantle_eta):

    r=np.sqrt(x**2+z**2)
    theta=np.pi/2-np.arctan2(z,x)

    R_c_l=R2-crust_depth
    R_l_um=R2-lithosphere_depth
    R_um_lm=R2-uppermantle_depth

    if r>R_c_l:
       val=crust_eta
    elif r>R_l_um:
       val=lithosphere_eta
    elif r>R_um_lm:
       val=uppermantle_eta
    else:
       val=lowermantle_eta

    if r>blob_R1 and r<blob_R2 and theta<blob_theta:
       val=blob_eta
  
    return val       

###############################################################################
