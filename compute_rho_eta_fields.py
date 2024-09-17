###############################################################################
import numpy as np
from density_Earth import *
from density_4DEarthBenchmark import *
from density_AnnulusBenchmark import *
from density_MarsDisc import *
from viscosity_Earth import *
from viscosity_4DEarthBenchmark import *
from viscosity_Mars import *
from viscosity_MarsDisc import *

###############################################################################
# this function computes the density and viscosity in the middle of the elt.
# It also projects these values onto the velocity nodes (for visualisation 
# purposes only!)
###############################################################################

def compute_rho_eta_fields(mV,NV,nel,xc,zc,iconV,R1,R2,rho_m,eta_m,eta_model,\
                           rho_model,rhoblob,etablob,zblob,Rblob,\
                           rhodisc,etadisc,R1disc,R2disc,thetadisc,\
                           eta_crust,eta_lithosphere,eta_uppermantle,eta_lowermantle,\
                           rho_crust,rho_lithosphere,rho_uppermantle,rho_lowermantle,\
                           R_c_l,R_l_um,R_um_lm,\
                           planet):

    counter=np.zeros(NV,dtype=np.float64)
    viscosity_elemental=np.zeros(nel,dtype=np.float64)
    density_elemental=np.zeros(nel,dtype=np.float64)
    viscosity_nodal=np.zeros(NV,dtype=np.float64)
    density_nodal=np.zeros(NV,dtype=np.float64)

    for iel in range(0,nel):
        match planet:
           case "Earth":
              viscosity_elemental[iel]=viscosity_Earth(xc[iel],zc[iel],R1,R2,eta_m,eta_model)
              density_elemental[iel]=density_Earth(xc[iel],zc[iel],R1,R2,rho_m,rho_model,\
                                                   rhoblob,zblob,Rblob)
           case "Mars":
              viscosity_elemental[iel]=viscosity_Mars(xc[iel],zc[iel],R1,R2,eta_m,eta_model)
              density_elemental[iel]=density_Mars(xc[iel],zc[iel],R1,R2,rho_m,rho_model,\
                                                  rhoblob,zblob,Rblob)
           case "4DEarthBenchmark":
              viscosity_elemental[iel]=viscosity_4DEarthBenchmark(xc[iel],zc[iel],R1,R2,eta_m,etablob,zblob,Rblob)
              density_elemental[iel]  =density_4DEarthBenchmark  (xc[iel],zc[iel],R1,R2,rho_m,rhoblob,zblob,Rblob)
                                                              
           case "AnnulusBenchmark":
              viscosity_elemental[iel]=1
              density_elemental[iel]=density_AnnulusBenchmark(xc[iel],zc[iel],R1,R2)

           case "AquariumBenchmark":
              viscosity_elemental[iel]=1e22
              density_elemental[iel]=4000

           case "MarsDisc":
              viscosity_elemental[iel]=viscosity_MarsDisc(xc[iel],zc[iel],R1,R2,etadisc,R1disc,R2disc,thetadisc,\
                                                          eta_crust,eta_lithosphere,eta_uppermantle,eta_lowermantle,\
                                                          R_c_l,R_l_um,R_um_lm)

              density_elemental[iel]  =density_MarsDisc(xc[iel],zc[iel],R1,R2,rhodisc,R1disc,R2disc,thetadisc,\
                                                        rho_crust,rho_lithosphere,rho_uppermantle,rho_lowermantle,\
                                                        R_c_l,R_l_um,R_um_lm)
                                                              
           case _:
              exit('pb in compute_rho_eta_fields: unknown planet')

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
