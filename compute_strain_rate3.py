###############################################################################
from basis_functions import *
import numpy as np

###############################################################################

def compute_strain_rate3(nel,mV,NV,nqel,iconV,mapping,xmapping,zmapping,u,v):

    exx3=np.zeros(NV,dtype=np.float64)  
    eyy3=np.zeros(NV,dtype=np.float64)  
    exy3=np.zeros(NV,dtype=np.float64)  

    M_mat=lil_matrix((NV,NV),dtype=np.float64)
    rhsLxx=np.zeros(NV,dtype=np.float64)
    rhsLyy=np.zeros(NV,dtype=np.float64)
    rhsLxy=np.zeros(NV,dtype=np.float64)
    rhsLyx=np.zeros(NV,dtype=np.float64)

    for iel in range(0,nel):

           M_el =np.zeros((mV,mV),dtype=np.float64)
           fLxx_el=np.zeros(mV,dtype=np.float64)
           fLyy_el=np.zeros(mV,dtype=np.float64)
           fLxy_el=np.zeros(mV,dtype=np.float64)
           fLyx_el=np.zeros(mV,dtype=np.float64)
           NNNV1 =np.zeros((mV,1),dtype=np.float64) 

           for kq in range(0,nqel):
               rq=qcoords_r[kq]
               sq=qcoords_s[kq]
               weightq=qweights[kq]

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
               jcbi=np.linalg.inv(jcb)

               #basis functions
               NNNV1[:,0]=NNN(rq,sq,'Q2')
               dNNNVdr=dNNNdr(rq,sq,'Q2')
               dNNNVds=dNNNds(rq,sq,'Q2')

               # compute dNdx & dNdy
               Lxxq=0.
               Lyyq=0.
               Lxyq=0.
               Lyxq=0.
               for k in range(0,mV):
                   dNNNVdx[k]=jcbi[0,0]*dNNNVdr[k]+jcbi[0,1]*dNNNVds[k]
                   dNNNVdy[k]=jcbi[1,0]*dNNNVdr[k]+jcbi[1,1]*dNNNVds[k]
                   Lxxq+=dNNNVdx[k]*u[iconV[k,iel]]
                   Lyyq+=dNNNVdy[k]*v[iconV[k,iel]]
                   Lxyq+=dNNNVdx[k]*v[iconV[k,iel]]
                   Lyxq+=dNNNVdy[k]*u[iconV[k,iel]]
               #end for 

               M_el +=NNNV1.dot(NNNV1.T)*weightq*jcob

               fLxx_el[:]+=NNNV1[:,0]*Lxxq*jcob*weightq
               fLyy_el[:]+=NNNV1[:,0]*Lyyq*jcob*weightq
               fLxy_el[:]+=NNNV1[:,0]*Lxyq*jcob*weightq
               fLyx_el[:]+=NNNV1[:,0]*Lyxq*jcob*weightq

           #end for kq

           for k1 in range(0,mV):
               m1=iconV[k1,iel]
               for k2 in range(0,mV):
                   m2=iconV[k2,iel]
                   M_mat[m1,m2]+=M_el[k1,k2]
               #end for
               rhsLxx[m1]+=fLxx_el[k1]
               rhsLyy[m1]+=fLyy_el[k1]
               rhsLxy[m1]+=fLxy_el[k1]
               rhsLyx[m1]+=fLyx_el[k1]
           #end for

    #end for
    #sparse_matrix=A_sparse.tocsr()

    Lxx3 = sps.linalg.spsolve(sps.csr_matrix(M_mat),rhsLxx)
    Lyy3 = sps.linalg.spsolve(sps.csr_matrix(M_mat),rhsLyy)
    Lxy3 = sps.linalg.spsolve(sps.csr_matrix(M_mat),rhsLxy)
    Lyx3 = sps.linalg.spsolve(sps.csr_matrix(M_mat),rhsLyx)

    exx3[:]=Lxx3[:]
    eyy3[:]=Lyy3[:]
    exy3[:]=0.5*(Lxy3[:]+Lyx3[:])

    return exx3,eyy3,exy3
###############################################################################
