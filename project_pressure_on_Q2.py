###############################################################################
import numpy as np
from basis_functions import *

###############################################################################

def project_pressure_on_Q2(NV,nel,mV,mP,p,iconP,iconV):

    count=np.zeros(NV,dtype=np.int32)
    q=np.zeros(NV,dtype=np.float64)

    rVnodes=[-1,1,1,-1,0,1,0,-1,0]
    sVnodes=[-1,-1,1,1,-1,0,1,0,0]

    for iel in range(0,nel):
        for i in range(0,mV):
            inode=iconV[i,iel]
            rq=rVnodes[i]
            sq=sVnodes[i]
            NNNP=NNN(rq,sq,'Q1')
            q[inode]+=np.dot(p[iconP[0:mP,iel]],NNNP[0:mP])
            count[inode]+=1
        #end for
    #end for
    q/=count

    return q

###############################################################################
