###############################################################################
import numpy as np

###############################################################################
# Q1 nodes for mapping are corners of Q2 basis functions
# Q2 nodes for mapping are same as Q2 basis functions
# Q3,Q4,Q5,Q6 nodes for mapping are built
# note that python uses row-major storage for 2S arrays and since we need 
# to do dot products with xmapping and zmapping it makes more sense 
# to use column-major, i.e. F_CONTIGUOUS in python jargon.
###############################################################################
    
def define_mapping(mapping,mmapping,xV,yV,iconV,nel,axisymmetric,rad,theta):

    xmapping=np.zeros((mmapping,nel),dtype=np.float64,order='F')
    zmapping=np.zeros((mmapping,nel),dtype=np.float64,order='F')

    if mapping=='Q1':
       for iel in range(0,nel):
           xmapping[0,iel]=xV[iconV[0,iel]] ; zmapping[0,iel]=yV[iconV[0,iel]]
           xmapping[1,iel]=xV[iconV[1,iel]] ; zmapping[1,iel]=yV[iconV[1,iel]]
           xmapping[2,iel]=xV[iconV[2,iel]] ; zmapping[2,iel]=yV[iconV[2,iel]]
           xmapping[3,iel]=xV[iconV[3,iel]] ; zmapping[3,iel]=yV[iconV[3,iel]]

    if mapping=='Q2':
       for iel in range(0,nel):
           xmapping[0,iel]=xV[iconV[0,iel]] ; zmapping[0,iel]=yV[iconV[0,iel]]
           xmapping[1,iel]=xV[iconV[1,iel]] ; zmapping[1,iel]=yV[iconV[1,iel]]
           xmapping[2,iel]=xV[iconV[2,iel]] ; zmapping[2,iel]=yV[iconV[2,iel]]
           xmapping[3,iel]=xV[iconV[3,iel]] ; zmapping[3,iel]=yV[iconV[3,iel]]
           xmapping[4,iel]=xV[iconV[4,iel]] ; zmapping[4,iel]=yV[iconV[4,iel]]
           xmapping[5,iel]=xV[iconV[5,iel]] ; zmapping[5,iel]=yV[iconV[5,iel]]
           xmapping[6,iel]=xV[iconV[6,iel]] ; zmapping[6,iel]=yV[iconV[6,iel]]
           xmapping[7,iel]=xV[iconV[7,iel]] ; zmapping[7,iel]=yV[iconV[7,iel]]
           xmapping[8,iel]=xV[iconV[8,iel]] ; zmapping[8,iel]=yV[iconV[8,iel]]

    if mapping=='Q3':
       if not axisymmetric:
          dtheta=2*np.pi/nelt/3
       else:
          dtheta=np.pi/nelt/3
       for iel in range(0,nel):
           thetamin=theta[iconV[0,iel]]
           rmin=rad[iconV[0,iel]]
           rmax=rad[iconV[2,iel]]
           counter=0
           for j in range(0,4):
               for i in range(0,4):
                   ttt=thetamin+i*dtheta
                   rrr=rmin+j*(rmax-rmin)/3
                   xmapping[counter,iel]=math.sin(ttt)*rrr
                   zmapping[counter,iel]=math.cos(ttt)*rrr
                   counter+=1

    if mapping=='Q4':
       if not axisymmetric:
          dtheta=2*np.pi/nelt/4
       else:
          dtheta=np.pi/nelt/4
       for iel in range(0,nel):
           thetamin=theta[iconV[0,iel]]
           rmin=rad[iconV[0,iel]]
           rmax=rad[iconV[2,iel]]
           counter=0
           for j in range(0,5):
               for i in range(0,5):
                   ttt=thetamin+i*dtheta
                   rrr=rmin+j*(rmax-rmin)/4
                   xmapping[counter,iel]=math.sin(ttt)*rrr
                   zmapping[counter,iel]=math.cos(ttt)*rrr
                   counter+=1

    if mapping=='Q5':
       if not axisymmetric:
          dtheta=2*np.pi/nelt/5
       else:
          dtheta=np.pi/nelt/5
       for iel in range(0,nel):
           thetamin=theta[iconV[0,iel]]
           rmin=rad[iconV[0,iel]]
           rmax=rad[iconV[2,iel]]
           counter=0
           for j in range(0,6):
               for i in range(0,6):
                   ttt=thetamin+i*dtheta
                   rrr=rmin+j*(rmax-rmin)/5
                   xmapping[counter,iel]=math.sin(ttt)*rrr
                   zmapping[counter,iel]=math.cos(ttt)*rrr
                   counter+=1

    if mapping=='Q6':
       if not axisymmetric:
          dtheta=2*np.pi/nelt/6
       else:
          dtheta=np.pi/nelt/6
       for iel in range(0,nel):
           thetamin=theta[iconV[0,iel]]
           rmin=rad[iconV[0,iel]]
           rmax=rad[iconV[2,iel]]
           counter=0
           for j in range(0,7):
               for i in range(0,7):
                   ttt=thetamin+i*dtheta
                   rrr=rmin+j*(rmax-rmin)/6
                   xmapping[counter,iel]=math.sin(ttt)*rrr
                   zmapping[counter,iel]=math.cos(ttt)*rrr
                   counter+=1

    return xmapping,zmapping

###############################################################################
