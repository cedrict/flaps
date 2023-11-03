import numpy as np

from density import *
from viscosity import *

###############################################################################

def compute_rho_eta_fields(mV,NV,nel,xc,yc,iconV,R1,R2,rho_m,eta_m,eta_model,\
                           rho_model,exp,rhoblobstar,zblob,Rblob):

    counter=np.zeros(NV,dtype=np.float64)
    viscosity_elemental=np.zeros(nel,dtype=np.float64)
    density_elemental=np.zeros(nel,dtype=np.float64)
    viscosity_nodal=np.zeros(NV,dtype=np.float64)
    density_nodal=np.zeros(NV,dtype=np.float64)

    for iel in range(0,nel):
           viscosity_elemental[iel]=viscosity(xc[iel],yc[iel],R1,R2,eta_m,eta_model)
           density_elemental[iel]=density(xc[iel],yc[iel],R1,R2,rho_m,rho_model,\
                                          exp,rhoblobstar,zblob,Rblob)
           for i in range(0,mV):
               counter[iconV[i,iel]]+=1
               viscosity_nodal[iconV[i,iel]]+=viscosity_elemental[iel]
               density_nodal[iconV[i,iel]]+=density_elemental[iel]
           #end for
    #end for
    density_nodal/=counter
    viscosity_nodal/=counter

    return  density_elemental,density_nodal,viscosity_elemental,viscosity_nodal

###############################################################################
