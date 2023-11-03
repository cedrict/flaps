###############################################################################
import numpy as np
from basis_functions import NNN

###############################################################################

def compute_element_center_coords(nel,mapping,xmapping,zmapping):

    xc=np.zeros(nel,dtype=np.float64)
    zc=np.zeros(nel,dtype=np.float64)
    thetac=np.zeros(nel,dtype=np.float64)

    for iel in range(0,nel):
        rq=0
        sq=0
        NNNV=NNN(rq,sq,mapping)
        xc[iel]=np.dot(NNNV[:],xmapping[:,iel])
        zc[iel]=np.dot(NNNV[:],zmapping[:,iel])
        thetac[iel]=np.pi/2-np.arctan2(zc[iel],xc[iel])
    #end for

    return xc,zc,thetac

###############################################################################
