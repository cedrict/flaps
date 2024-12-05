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

def compute_rho_eta_fields(mV,NV,nel,xV,zV,xc,zc,iconV,R1,R2,eta_model,rho_model,\
                           blob_rho,blob_eta,blob_z,blob_R,blob_R1,blob_R2,blob_theta,\
                           crust_rho,crust_eta,crust_depth,\
                           lithosphere_rho,lithosphere_eta,lithosphere_depth,\
                           uppermantle_rho,uppermantle_eta,uppermantle_depth,\
                           lowermantle_rho,lowermantle_eta,\
                           planet):


    ######################
    # elemental quantities

    viscosity_elemental=np.zeros(nel,dtype=np.float64)
    density_elemental=np.zeros(nel,dtype=np.float64)

    for iel in range(0,nel):
        match planet:
           case "Earth":
              viscosity_elemental[iel]=viscosity_Earth(xc[iel],zc[iel],R1,R2,eta_m,eta_model)
              density_elemental[iel]=density_Earth(xc[iel],zc[iel],R1,R2,rho_m,rho_model,\
                                                   blob_rho,blob_z,blob_R)
           case "Mars":
              viscosity_elemental[iel]=viscosity_Mars(xc[iel],zc[iel],R1,R2,eta_m,eta_model)
              density_elemental[iel]=density_Mars(xc[iel],zc[iel],R1,R2,rho_m,rho_model,\
                                                  blob_rho,blob_z,blob_R)
           case "4DEarthBenchmark":
              viscosity_elemental[iel]=viscosity_4DEarthBenchmark(xc[iel],zc[iel],R1,R2,crust_eta,lithosphere_eta,\
                                                                  uppermantle_eta,lowermantle_eta,blob_eta,blob_z,blob_R)
              density_elemental[iel]=density_4DEarthBenchmark(xc[iel],zc[iel],R1,R2,crust_rho,lithosphere_rho,\
                                                              uppermantle_rho,lowermantle_rho,blob_rho,blob_z,blob_R)
                                                              
           case "AnnulusBenchmark":
              viscosity_elemental[iel]=1
              density_elemental[iel]=density_AnnulusBenchmark(xc[iel],zc[iel],R1,R2)

           case "AquariumBenchmark":
              viscosity_elemental[iel]=1e22
              density_elemental[iel]=4000

           case "MarsDisc":
              viscosity_elemental[iel]=viscosity_MarsDisc(xc[iel],zc[iel],R1,R2,\
                                                          blob_eta,blob_z,blob_R,blob_R1,blob_R2,blob_theta,\
                                                          crust_eta,crust_depth,\
                                                          lithosphere_eta,lithosphere_depth,\
                                                          uppermantle_eta,uppermantle_depth,\
                                                          lowermantle_eta)
                                                              


              density_elemental[iel]=density_MarsDisc(xc[iel],zc[iel],R1,R2,\
                                                      blob_rho,blob_z,blob_R,blob_R1,blob_R2,blob_theta,\
                                                      crust_rho,crust_depth,\
                                                      lithosphere_rho,lithosphere_depth,\
                                                      uppermantle_rho,uppermantle_depth,\
                                                      lowermantle_rho)
                                                              
           case _:
              exit('pb1 in compute_rho_eta_fields: unknown planet')

    #end for

    ######################
    # nodal quantities

    viscosity_nodal=np.zeros(NV,dtype=np.float64)
    density_nodal=np.zeros(NV,dtype=np.float64)

    for i in range(0,NV):

        match planet:
           case "4DEarthBenchmark":
              viscosity_nodal[i]=viscosity_4DEarthBenchmark(xV[i],zV[i],R1,R2,crust_eta,lithosphere_eta,\
                                                            uppermantle_eta,lowermantle_eta,blob_eta,blob_z,blob_R)
              density_nodal[i]=density_4DEarthBenchmark(xV[i],zV[i],R1,R2,crust_rho,lithosphere_rho,\
                                                        uppermantle_rho,lowermantle_rho,blob_rho,blob_z,blob_R)

           case "MarsDisc":
              viscosity_nodal[i]=viscosity_MarsDisc(xV[i],zV[i],R1,R2,\
                                                    blob_eta,blob_z,blob_R,blob_R1,blob_R2,blob_theta,\
                                                    crust_eta,crust_depth,\
                                                    lithosphere_eta,lithosphere_depth,\
                                                    uppermantle_eta,uppermantle_depth,\
                                                    lowermantle_eta)

              density_nodal[i]=density_MarsDisc(xV[i],zV[i],R1,R2,\
                                                blob_rho,blob_z,blob_R,blob_R1,blob_R2,blob_theta,\
                                                crust_rho,crust_depth,\
                                                lithosphere_rho,lithosphere_depth,\
                                                uppermantle_rho,uppermantle_depth,\
                                                lowermantle_rho)

           case _:
              exit('pb2 in compute_rho_eta_fields: unknown planet')

    #end for

    return  density_elemental,density_nodal,viscosity_elemental,viscosity_nodal

###############################################################################
