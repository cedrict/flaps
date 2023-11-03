###############################################################################
import numpy as np
from basis_functions import *
from analytical_solution import *

###############################################################################

def compute_errors(nel,nqel,mapping,xmapping,zmapping,qcoords_r,qcoords_s,qweights,\
                   R1,R2,rho_m,g0,exp,u,v,p,axisymmetric,iconV,iconP):

    jcb=np.zeros((2,2),dtype=np.float64)
    vrms=0.
    errv=0.
    errp=0.
    #errexx1=0. ; erreyy1=0. ; errexy1=0.
    #errexx2=0. ; erreyy2=0. ; errexy2=0.
    #errexx3=0. ; erreyy3=0. ; errexy3=0.
    for iel in range (0,nel):

        for kq in range(0,nqel):
            rq=qcoords_r[kq]
            sq=qcoords_s[kq]

            #compute coords of quadrature points
            NNNV=NNN(rq,sq,mapping)
            xq=np.dot(NNNV[:],xmapping[:,iel])
            yq=np.dot(NNNV[:],zmapping[:,iel])

            #compute jacobian matrix
            dNNNVdr=dNNNdr(rq,sq,mapping)
            dNNNVds=dNNNds(rq,sq,mapping)
            jcb[0,0]=np.dot(dNNNVdr[:],xmapping[:,iel])
            jcb[0,1]=np.dot(dNNNVdr[:],zmapping[:,iel])
            jcb[1,0]=np.dot(dNNNVds[:],xmapping[:,iel])
            jcb[1,1]=np.dot(dNNNVds[:],zmapping[:,iel])
            jcob=np.linalg.det(jcb)
            JxW=jcob*qweights[kq]

            #basis functions
            NNNV=NNN(rq,sq,'Q2')
            dNNNVdr=dNNNdr(rq,sq,'Q2')
            dNNNVds=dNNNds(rq,sq,'Q2')
            NNNP=NNN(rq,sq,'Q1')

            uq=np.dot(NNNV[:],u[iconV[:,iel]])
            vq=np.dot(NNNV[:],v[iconV[:,iel]])

            errv+=((uq-velocity_x(xq,yq,R1,R2,exp))**2+\
                   (vq-velocity_y(xq,yq,R1,R2,exp))**2)*JxW

            #exx1q=np.dot(NNNV[:],exx1[iconV[:,iel]])
            #eyy1q=np.dot(NNNV[:],eyy1[iconV[:,iel]])
            #exy1q=np.dot(NNNV[:],exy1[iconV[:,iel]])
            #exx2q=np.dot(NNNV[:],exx2[iconV[:,iel]])
            #eyy2q=np.dot(NNNV[:],eyy2[iconV[:,iel]])
            #exy2q=np.dot(NNNV[:],exy2[iconV[:,iel]])
            #exx3q=np.dot(NNNV[:],exx3[iconV[:,iel]])
            #eyy3q=np.dot(NNNV[:],eyy3[iconV[:,iel]])
            #exy3q=np.dot(NNNV[:],exy3[iconV[:,iel]])
            #errexx1+=(exx1q-sr_xx(xq,yq,R1,R2,kk))**2*JxW
            #erreyy1+=(eyy1q-sr_yy(xq,yq,R1,R2,kk))**2*JxW
            #errexy1+=(exy1q-sr_xy(xq,yq,R1,R2,kk))**2*JxW
            #errexx2+=(exx2q-sr_xx(xq,yq,R1,R2,kk))**2*JxW
            #erreyy2+=(eyy2q-sr_yy(xq,yq,R1,R2,kk))**2*JxW
            #errexy2+=(exy2q-sr_xy(xq,yq,R1,R2,kk))**2*JxW
            #errexx3+=(exx3q-sr_xx(xq,yq,R1,R2,kk))**2*JxW
            #erreyy3+=(eyy3q-sr_yy(xq,yq,R1,R2,kk))**2*JxW
            #errexy3+=(exy3q-sr_xy(xq,yq,R1,R2,kk))**2*JxW

            if axisymmetric: 
               JxW*=2*np.pi*xq

            vrms+=(uq**2+vq**2)*JxW

            pq=np.dot(NNNP[:],p[iconP[:,iel]])
            errp+=(pq-pressure(xq,yq,R1,R2,rho_m,g0,exp))**2*JxW

        # end for kq
    # end for iel

    errv=np.sqrt(errv)
    errp=np.sqrt(errp)

    if axisymmetric :
       vrms=np.sqrt(vrms/ (4/3*np.pi*(R2**3-R1**3)) )
    else:
       vrms=np.sqrt(vrms/ (np.pi*(R2**2-R1**2)) )

    return errv,errp,vrms

###############################################################################
