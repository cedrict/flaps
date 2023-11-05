###############################################################################
import numpy as np
from basis_functions import *
from analytical_solution import *

###############################################################################

def compute_errors(nel,nqel,mapping,xmapping,zmapping,qcoords_r,qcoords_s,qweights,\
                   R1,R2,rho_m,g0,exp,u,v,p,axisymmetric,iconV,iconP,total_volume):

    jcb=np.zeros((2,2),dtype=np.float64)
    vrms=0.
    errv=0.
    errp=0.
    for iel in range (0,nel):

        for kq in range(0,nqel):
            rq=qcoords_r[kq]
            sq=qcoords_s[kq]

            #compute coords of quadrature points
            NNNV=NNN(rq,sq,mapping)
            xq=np.dot(NNNV[:],xmapping[:,iel])
            zq=np.dot(NNNV[:],zmapping[:,iel])

            #compute jacobian matrix
            dNNNVdr=dNNNdr(rq,sq,mapping)
            dNNNVds=dNNNds(rq,sq,mapping)
            jcb[0,0]=np.dot(dNNNVdr[:],xmapping[:,iel])
            jcb[0,1]=np.dot(dNNNVdr[:],zmapping[:,iel])
            jcb[1,0]=np.dot(dNNNVds[:],xmapping[:,iel])
            jcb[1,1]=np.dot(dNNNVds[:],zmapping[:,iel])
            jcob=np.linalg.det(jcb)
            JxW=jcob*qweights[kq]
            if axisymmetric: JxW*=2*np.pi*xq

            NNNV=NNN(rq,sq,'Q2')
            dNNNVdr=dNNNdr(rq,sq,'Q2')
            dNNNVds=dNNNds(rq,sq,'Q2')
            NNNP=NNN(rq,sq,'Q1')

            uq=np.dot(NNNV[:],u[iconV[:,iel]])
            vq=np.dot(NNNV[:],v[iconV[:,iel]])

            errv+=((uq-velocity_x(xq,zq,R1,R2,exp))**2+\
                   (vq-velocity_y(xq,zq,R1,R2,exp))**2)*JxW

            vrms+=(uq**2+vq**2)*JxW

            pq=np.dot(NNNP[:],p[iconP[:,iel]])
            errp+=(pq-pressure(xq,zq,R1,R2,rho_m,g0,exp))**2*JxW

        # end for kq
    # end for iel

    errv=np.sqrt(errv)
    errp=np.sqrt(errp)

    vrms=np.sqrt(vrms/total_volume)

    return errv,errp,vrms

###############################################################################
