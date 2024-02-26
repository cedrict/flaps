###############################################################################
from basis_functions import *
import numpy as np

###############################################################################

def compute_strain_rate1(nel,mV,NV,iconV,mapping,xmapping,zmapping,u,v):

    exxc=np.zeros(nel,dtype=np.float64)
    eyyc=np.zeros(nel,dtype=np.float64)
    exyc=np.zeros(nel,dtype=np.float64)
    exx1=np.zeros(NV,dtype=np.float64)  
    exy1=np.zeros(NV,dtype=np.float64)  
    eyy1=np.zeros(NV,dtype=np.float64)  
    count=np.zeros(NV,dtype=np.int32)  
    jcb=np.zeros((2,2),dtype=np.float64)
    dNNNVdx=np.zeros(mV,dtype=np.float64)
    dNNNVdy=np.zeros(mV,dtype=np.float64)

    for iel in range(0,nel):

           rq=0
           sq=0

           #compute jacobian matrix
           dNNNVdr=dNNNdr(rq,sq,mapping)
           dNNNVds=dNNNds(rq,sq,mapping)
           jcb[0,0]=np.dot(dNNNVdr[:],xmapping[:,iel])
           jcb[0,1]=np.dot(dNNNVdr[:],zmapping[:,iel])
           jcb[1,0]=np.dot(dNNNVds[:],xmapping[:,iel])
           jcb[1,1]=np.dot(dNNNVds[:],zmapping[:,iel])
           jcbi=np.linalg.inv(jcb)

           #basis functions
           dNNNVdr=dNNNdr(rq,sq,'Q2')
           dNNNVds=dNNNds(rq,sq,'Q2')

           for k in range(0,mV):
               dNNNVdx[k]=jcbi[0,0]*dNNNVdr[k]+jcbi[0,1]*dNNNVds[k]
               dNNNVdy[k]=jcbi[1,0]*dNNNVdr[k]+jcbi[1,1]*dNNNVds[k]
           #end for

           e_xx=np.dot(dNNNVdx[:],u[iconV[:,iel]])
           e_xy=np.dot(dNNNVdx[:],v[iconV[:,iel]])*0.5\
               +np.dot(dNNNVdy[:],u[iconV[:,iel]])*0.5
           e_yy=np.dot(dNNNVdy[:],v[iconV[:,iel]])

           exxc[iel]=e_xx
           eyyc[iel]=e_yy
           exyc[iel]=e_xy

           for i in range(0,mV):
               inode=iconV[i,iel]
               exx1[inode]+=e_xx
               exy1[inode]+=e_xy
               eyy1[inode]+=e_yy
               count[inode]+=1
           #end for

    #end for
    exx1/=count
    exy1/=count
    eyy1/=count

    return exx1,eyy1,exy1,exxc,eyyc,exyc

###############################################################################
