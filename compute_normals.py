import numpy as np

def compute_normals(normal_type,NV,xV,yV,mV,nqel,surfaceV,hull,\
                    R1,R2,eps,iconV,axisymmetric,rad,theta,\
                    qcoords_r,qcoords_s,qweights,):

   nx=np.zeros(NV,dtype=np.float64) 
   ny=np.zeros(NV,dtype=np.float64) 

   if normal_type==1:
      nx1=np.zeros(NV,dtype=np.float64) 
      ny1=np.zeros(NV,dtype=np.float64) 
      for i in range(0,NV):
          if axisymmetric and xV[i]/R2<eps:
             nx1[i]=-1
             ny1[i]=0
          if rad[i]/R1<1+eps:
             nx1[i]=-np.sin(theta[i])
             ny1[i]=-np.cos(theta[i])
          if rad[i]/R2>1-eps:
             nx1[i]=np.sin(theta[i])
             ny1[i]=np.cos(theta[i])
      #end for

      if axisymmetric:
         for i in range(0,NV):
             if xV[i]/R1<eps and yV[i]/R2>1-eps:
                nx1[i]=0
                ny1[i]=0
             if xV[i]/R1<eps and yV[i]/R2<-1+eps:
                nx1[i]=0
                ny1[i]=0
         #end for
      #end if
      nx[:]=nx1[:]
      ny[:]=ny1[:]


   if normal_type==2:
      exit('this normal type calulation requires xy-mapping')
      nx2=np.zeros(NV,dtype=np.float64) 
      ny2=np.zeros(NV,dtype=np.float64) 
      dNNNVdx=np.zeros(mV,dtype=np.float64) 
      dNNNVdy=np.zeros(mV,dtype=np.float64) 
      jcb=np.zeros((2,2),dtype=np.float64)
      for iel in range(0,nel):
          if True: #hull[iconV[0,iel]] or hull[iconV[2,iel]]: 
             for kq in range(0,nqel):
                 rq=qcoords_r[kq]
                 sq=qcoords_s[kq]
                 weightq=qweights[kq]
                 #compute jacobian matrix
                 dNNNVdr=dNNNdr(rq,sq,mapping)
                 dNNNVds=dNNNds(rq,sq,mapping)
                 jcb[0,0]=np.dot(dNNNVdr[:],xmapping[:,iel])
                 jcb[0,1]=np.dot(dNNNVdr[:],ymapping[:,iel])
                 jcb[1,0]=np.dot(dNNNVds[:],xmapping[:,iel])
                 jcb[1,1]=np.dot(dNNNVds[:],ymapping[:,iel])
                 jcob=np.linalg.det(jcb)
                 jcbi=np.linalg.inv(jcb)
                 #basis functions
                 dNNNVdr=dNNNdr(rq,sq,'Q2')
                 dNNNVds=dNNNds(rq,sq,'Q2')
                 # compute dNdx & dNdy
                 for k in range(0,mV):
                     dNNNVdx[k]=jcbi[0,0]*dNNNVdr[k]+jcbi[0,1]*dNNNVds[k]
                     dNNNVdy[k]=jcbi[1,0]*dNNNVdr[k]+jcbi[1,1]*dNNNVds[k]
                 #end for 
                 nx2[iconV[0,iel]]+=dNNNVdx[0]*jcob*weightq
                 ny2[iconV[0,iel]]+=dNNNVdy[0]*jcob*weightq
                 nx2[iconV[1,iel]]+=dNNNVdx[1]*jcob*weightq
                 ny2[iconV[1,iel]]+=dNNNVdy[1]*jcob*weightq
                 nx2[iconV[2,iel]]+=dNNNVdx[2]*jcob*weightq
                 ny2[iconV[2,iel]]+=dNNNVdy[2]*jcob*weightq
                 nx2[iconV[3,iel]]+=dNNNVdx[3]*jcob*weightq
                 ny2[iconV[3,iel]]+=dNNNVdy[3]*jcob*weightq
                 nx2[iconV[4,iel]]+=dNNNVdx[4]*jcob*weightq
                 ny2[iconV[4,iel]]+=dNNNVdy[4]*jcob*weightq
                 nx2[iconV[5,iel]]+=dNNNVdx[5]*jcob*weightq
                 ny2[iconV[5,iel]]+=dNNNVdy[5]*jcob*weightq
                 nx2[iconV[6,iel]]+=dNNNVdx[6]*jcob*weightq
                 ny2[iconV[6,iel]]+=dNNNVdy[6]*jcob*weightq
                 nx2[iconV[7,iel]]+=dNNNVdx[7]*jcob*weightq
                 ny2[iconV[7,iel]]+=dNNNVdy[7]*jcob*weightq
             #end for
          #end if
      #end for

      for i in range(0,NV):
          if hull[i]:
             norm=np.sqrt(nx2[i]**2+ny2[i]**2)
             nx2[i]/=norm
             ny2[i]/=norm


      if axisymmetric:
         for i in range(0,NV):
             if xV[i]/R1<eps and yV[i]/R2>1-eps:
                nx2[i]=0
                ny2[i]=0
             if xV[i]/R1<eps and yV[i]/R2<-1+eps:
                nx2[i]=0
                ny2[i]=0
   
      nx[:]=nx2[:]
      ny[:]=ny2[:]

   return nx,ny




