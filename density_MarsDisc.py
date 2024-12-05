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
def density_MarsDisc(x,z,R1,R2,\
                     blob_rho,blob_z,blob_R,blob_R1,blob_R2,blob_theta,\
                     crust_rho,crust_depth,\
                     lithosphere_rho,lithosphere_depth,\
                     uppermantle_rho,uppermantle_depth,\
                     lowermantle_rho):

    r=np.sqrt(x**2+z**2)
    theta=np.pi/2-np.arctan2(z,x)

    R_c_l=R2-crust_depth
    R_l_um=R2-lithosphere_depth
    R_um_lm=R2-uppermantle_depth

    if r>R_c_l:
       val=crust_rho
    elif r>R_l_um:
       val=lithosphere_rho
    elif r>R_um_lm:
       val=uppermantle_rho
    else:
       val=lowermantle_rho

    if r>blob_R1 and r<blob_R2 and theta<blob_theta:
       val=blob_rho

    #val=3550
  
    return val       

###############################################################################
