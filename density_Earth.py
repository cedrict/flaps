###############################################################################
#rho_model=0 # constant
#rho_model=1 # PREM
#rho_model=2 # ST105W
#rho_model=3 # AK135F
###############################################################################

from numba import jit
import numpy as np
import math as math

#@jit(nopython=True)
def density_Earth(x,y,R1,R2,rho_m,rho_model,rhoblob,yblob,Rblob):

    match rho_model:

       #------
       case 0:
          val=rho_m

       #------------
       case 1: #PREM
          radius=np.sqrt(x*x+y*y)
          xx=radius/6371.e3
          if radius>6371e3:
             densprem=0
          elif radius<=1221.5e3:
              densprem=13.0885-8.8381*xx**2
          elif radius<=3480e3:
              densprem=12.5815-1.2638*xx-3.6426*xx**2-5.5281*xx**3
          elif radius<=3630.e3:
             densprem=7.9565-6.4761*xx+5.5283*xx**2-3.0807*xx**3
          elif radius<=5600.e3:
             densprem=7.9565-6.4761*xx+5.5283*xx**2-3.0807*xx**3
          elif radius<=5701.e3:
             densprem=7.9565-6.4761*xx+5.5283*xx**2-3.0807*xx**3
          elif radius<=5771.e3:
             densprem=5.3197-1.4836*xx
          elif radius<=5971.e3:
             densprem=11.2494-8.0298*xx
          elif radius<=6151.e3:
             densprem=7.1089-3.8045*xx
          elif radius<=6291.e3:
             densprem=2.6910+0.6924*xx
          elif radius<=6346.e3:
             densprem=2.6910+0.6924*xx
          elif radius<=6356.e3:
             densprem=2.9
          elif radius<=6368.e3:
             densprem=2.6
          else:
             densprem=1.020
          val=densprem*1000

       #------------
       case 2: #ST105W

          depth=R2-np.sqrt(x**2+y**2)
          cell_index=49
          #for kk in range(0,50):
          #    if depth<depth_st105w[kk+1]:
          #       cell_index=kk
          #       break
          #    #end if
          #end for
          #val=(depth-depth_st105w[cell_index])/(depth_st105w[cell_index+1]-depth_st105w[cell_index])\
          #   *(rho_st105w[cell_index+1]-rho_st105w[cell_index])+rho_st105w[cell_index]

       #------------
       case 3: #AK135F

          depth=R2-np.sqrt(x**2+y**2)
          cell_index=144
          #for kk in range(0,145):
          #    if depth<depth_ak135f[kk+1]:
          #       cell_index=kk
          #       break
          #    #end if
          #end for
          #val=(depth-depth_ak135f[cell_index])/(depth_ak135f[cell_index+1]-depth_ak135f[cell_index])\
          #   *(rho_ak135f[cell_index+1]-rho_ak135f[cell_index])+rho_ak135f[cell_index]

       #------------
       case _:
          exit('density_Earth: unknown rho_model')

    return val

###############################################################################
