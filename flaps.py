import numpy as np
import math as math
import sys as sys
import scipy
import scipy.sparse as sps
from scipy.sparse.linalg import *
import time as timing
import matplotlib.pyplot as plt
from scipy import sparse
from numba import jit

###############################################################################
# list of other python files containing the code
###############################################################################

from basis_functions import *
from density import *
from viscosity import *
from analytical_solution import *
from gravity_vector import *
from compute_gravity_at_point import *
from compute_strain_rate2 import *
from export_to_vtu import *
from quadrature_setup import *
from compute_normals import *
from define_mapping import *
from compute_element_center_coords import *
from compute_rho_eta_fields import *
from compute_errors import *
from project_pressure_on_Q2 import *

###############################################################################

axisymmetric=True

surface_free_slip=True

bottom_free_slip=True # not yet implemented

solve_stokes=True

use_elemental_density=False
use_elemental_viscosity=False

self_gravitation=False

###############################################################################
# list of available Earth density and viscosity profiles
###############################################################################

planet_is_Earth=False

rho_model=0 # constant
#rho_model=1 # PREM
#rho_model=2 # STW105
#rho_model=3 # AK135F

eta_model=0 # constant
#eta_model=1 # civs12 a
#eta_model=2 # civs12 b
#eta_model=3 # stho08
#eta_model=4 # yohk01

###############################################################################
# list of available Mars density and viscosity profiles
###############################################################################

planet_is_Mars=False

use_eta_steinberger=False
use_eta_samuelA=False
use_eta_samuelB=False

###############################################################################

planet_is_4DEarth=False

###############################################################################

spacing="                                                 "

year=365.25*3600*24
Ggrav=6.67430e-11
mGal=1e-5

print("+=================================================")
print("+=================================================")
print("+---------------------FL.A.P.S--------------------")
print("+=================================================")
print("+=================================================")

ndim=2   # number of dimensions
mV=9     # number of nodes making up an element
mP=4     # number of nodes making up an element
ndofV=2  # number of velocity degrees of freedom per node
ndofP=1  # number of pressure degrees of freedom 

if int(len(sys.argv) == 14):
   exp         = int(sys.argv[1])
   nelr        = int(sys.argv[2])
   visu        = int(sys.argv[3])
   nqperdim    = int(sys.argv[4])
   mapping     = int(sys.argv[5])
   xi          = int(sys.argv[6])
   etablobstar = float(sys.argv[7])
   rhoblobstar = float(sys.argv[8])
   zblob       = float(sys.argv[9])
   Rblob       = float(sys.argv[10])
   etalithstar = float(sys.argv[11])
   compute_g   = int(sys.argv[12])
   nel_phi     = int(sys.argv[13])
   if mapping==1: mapping='Q1'
   if mapping==2: mapping='Q2'
   if mapping==3: mapping='Q3'
   if mapping==4: mapping='Q4'
   if mapping==5: mapping='Q5'
   if mapping==6: mapping='Q6'
   compute_gravity=(compute_g==1)
   print(sys.argv)

else:
   # exp=0: cyl benchmark
   # exp=1: aquarium
   # exp=2: blob
   exp         = 1
   nelr        = 128 # Q1 cells!
   visu        = 1
   nqperdim    = 3
   mapping     = 'Q2' 
   xi          = 6
   etablobstar = 1
   rhoblobstar = 0.9
   zblob       = 4900e3
   Rblob       = 400e3
   etalithstar = 1
   compute_gravity=False
   nel_phi     = 100

gravity_model=2
rho_c=6000

dt=50*year
np_grav=48 # nb of satellites positions
height=250e3

normal_type=1

if exp==0:
   R1=1.
   R2=2.
   rho_m=0.
   g0=1.
   eta_ref=1
   eta_m=1
   vel_unit=1
   velunit=' '
   surface_free_slip=False
   compute_gravity=False
   axisymmetric=False
   nstep=1
   rhoc=0
else:
   R1=3400e3
   R2=6400e3
   g0=10.
   eta_ref=1e21
   eta_m=1e21
   rho_m=4000.
   vel_unit=0.01/year
   velunit='cm/year'
   nstep=1
   rhoc=6000

rhoblob=rhoblobstar*rho_m

eps=1.e-10

debug=False

compute_sr1=False
compute_sr3=False
    
if mapping=='Q1': mmapping=4
if mapping=='Q2': mmapping=9
if mapping=='Q3': mmapping=16
if mapping=='Q4': mmapping=25
if mapping=='Q5': mmapping=36
if mapping=='Q6': mmapping=49

###############################################################################
###############################################################################
start = timing.time()

if planet_is_Earth:
   R1=3480e3
   R2=6371e3
   g0=9.81

   print('max depth',R2-R1)

   #--------------------------------------
   if eta_model==1 or eta_model==2: #civs12a,b
      # reading data from civs12
      # file is 153 lines long 
      # first 51 lines are viscA, then 51 lines are viscB 
      # and last 51 lines are depths, from 0 to 2900km 
      # I have removed all ",&"

      viscA_civs12 = np.empty(51,dtype=np.float64)
      viscB_civs12 = np.empty(51,dtype=np.float64)
      depths_civs12 = np.empty(51,dtype=np.float64)

      f = open('DATA/EARTH/eta_civs12.ascii','r')
      counter=0
      for line in f:
          line=line.strip()
          columns=line.split()
          if counter<51:
             viscA_civs12[counter]=columns[0]
          elif counter<102:
             viscB_civs12[counter-51]=columns[0]
          else:
             depths_civs12[counter-102]=columns[0]
          counter+=1

      depths_civs12[:]=np.flip(depths_civs12)
      viscA_civs12[:]=np.flip(viscA_civs12)
      viscB_civs12[:]=np.flip(viscB_civs12)

      np.savetxt('civs12.ascii',np.array([depths_civs12,viscA_civs12,viscB_civs12]).T)

      print('     -> read eta_civs12.ascii ok') 

   #-----------------
   if eta_model==3: #stho08
      # reading data from  Steinberger & Holmes 2008
      # file counts 22 lines
      # first column is number between 0 and 1 (normalised radii)
      # second column is viscosity
      # I have added a last line R=1, eta=1e24

      depths_stho08 = np.empty(23,dtype=np.float64)
      visc_stho08 = np.empty(23,dtype=np.float64)
      f = open('DATA/EARTH/eta_stho08.ascii','r')
      counter=0
      for line in f:
          line=line.strip()
          columns=line.split()
          depths_stho08[counter]=columns[0]
          visc_stho08[counter]=columns[1]
          counter+=1
      #end for
      depths_stho08[:]=6371e3*(1-depths_stho08[:])
      depths_stho08[:]=np.flip(depths_stho08)
      visc_stho08[:]=np.flip(visc_stho08)

      np.savetxt('stho08.ascii',np.array([depths_stho08,visc_stho08]).T)

      print('     -> read eta_stho08.ascii ok')

   #-----------------
   if rho_model==3: #AK135F:
      depth_ak135f = np.empty(145,dtype=np.float64)
      rho_ak135f = np.empty(145,dtype=np.float64)
      f = open('DATA/EARTH/rho_ak135f.ascii','r')
      counter=0
      for line in f:
          line=line.strip()
          columns=line.split()
          depth_ak135f[counter]=columns[0]
          rho_ak135f[counter]=columns[1]
          counter+=1
      #end for
      depth_ak135f[:]*=1000
      rho_ak135f[:]*=1000

      np.savetxt('density_ak135f.ascii',np.array([depth_ak135f,rho_ak135f]).T)

      print('     -> read rho_ak135f.ascii ok')

   #-----------------
   if rho_model==2: #STW105
      depth_stw105 = np.empty(750,dtype=np.float64)
      rho_stw105 = np.empty(750,dtype=np.float64)
      f = open('DATA/EARTH/rho_stw105.ascii','r')
      counter=0
      for line in f:
          line=line.strip()
          columns=line.split()
          depth_stw105[counter]=columns[0]
          rho_stw105[counter]=columns[1]
          counter+=1
      #end for

      depth_stw105[:]=6371e3-depth_stw105[:]
      depth_stw105[:]=np.flip(depth_stw105)
      rho_stw105[:]=np.flip(rho_stw105)

      np.savetxt('density_stw105.ascii',np.array([depth_stw105,rho_stw105]).T)

      print('     -> read rho_stw105.ascii ok')


   #--------------------------------------
   if grav_model==1: #prem
      depth_prem= np.empty(44,dtype=np.float64)
      g_prem = np.empty(44,dtype=np.float64)

print("read EARTH data...........................(%.3fs)" % (timing.time() - start))

###############################################################################
###############################################################################
start = timing.time()

if planet_is_Mars:
   R1=1700e3
   R2=3390e3
   g0=3.71 #https://en.wikipedia.org/wiki/Mars   

   if use_steinberger:
      R1=R2-1967e3
      R_disc1 = R2-49.5e3
      R_disc2 = R2-1111.5e3
      R_disc3 = R2-1160e3
      R_disc4 = R2-1951.5e3
      eta_max=1e25

   if use_samuelA:
      R1=1839.5976879540331e3
      R_disc1 = 3317.7417442558781e3
      R_disc2 = 2836.6008937146739e3
      R_disc3 = 2350.4998282194360e3
      R_disc4 = 1918.9611272185618e3
      eta_max=1e25

   if use_samuelB:
      R1=1624.2975658322634e3
      R_disc1 = 3324.3388640802909e3
      R_disc2 = 3090.3851276356227e3
      R_disc3 = 2313.0549710614014e3
      R_disc4 = 1822.5068139999998e3
      eta_max=1e25

print("read MARS data............................(%.3fs)" % (timing.time() - start))

###############################################################################
###############################################################################
start = timing.time()

if planet_is_4DEarth:

   R1=3400e3
   R2=6400e3


###############################################################################
# quadrature points and associated weights
###############################################################################
start = timing.time()

nqel,qcoords_r,qcoords_s,qweights=quadrature_setup(nqperdim)

print("quadrature setup..........................(%.3fs)" % (timing.time() - start))

###############################################################################
# grid point setup
###############################################################################
start = timing.time()

if not axisymmetric:

   nelt=xi*nelr 
   nel=nelr*nelt  
   nnr=nelr+1
   nnt=nelt
   NV=nnr*nnt  # number of V nodes

   xV=np.zeros(NV,dtype=np.float64) 
   yV=np.zeros(NV,dtype=np.float64) 
   rad=np.zeros(NV,dtype=np.float64)  
   theta=np.zeros(NV,dtype=np.float64) 

   Louter=2.*math.pi*R2
   Lr=R2-R1
   sx=Louter/float(nelt)
   sz=Lr/float(nelr)

   counter=0
   for j in range(0,nnr):
       for i in range(0,nelt):
           xV[counter]=i*sx
           yV[counter]=j*sz
           counter += 1

   counter=0
   for j in range(0,nnr):
       for i in range(0,nnt):
           xi=xV[counter]
           yi=yV[counter]
           t=xi/Louter*2.*math.pi    
           xV[counter]=math.cos(t)*(R1+yi)
           yV[counter]=math.sin(t)*(R1+yi)
           rad[counter]=R1+yi
           #theta[counter]=np.arctan2(yV[counter],xV[counter])
           theta[counter]=np.pi/2-np.arctan2(yV[counter],xV[counter])
           if theta[counter]<0.:
              theta[counter]+=2.*math.pi
           counter+=1

else:

   nelt=xi*nelr 
   nel=nelr*nelt  
   nnr=nelr+1
   nnt=nelt+1
   NV=nnr*nnt  # number of V nodes

   xV=np.zeros(NV,dtype=np.float64) 
   yV=np.zeros(NV,dtype=np.float64) 
   rad=np.zeros(NV,dtype=np.float64)  
   theta=np.zeros(NV,dtype=np.float64) 

   Louter=math.pi*R2
   Lr=R2-R1
   sx=Louter/float(nelt)
   sz=Lr/float(nelr)

   counter=0
   for j in range(0,nnr):
       for i in range(0,nnt):
           xV[counter]=i*sx
           yV[counter]=j*sz
           counter += 1

   counter=0
   for j in range(0,nnr):
       for i in range(0,nnt):
           x_i=xV[counter]
           yi=yV[counter]
           t=math.pi/2-x_i/Louter*math.pi 
           xV[counter]=math.cos(t)*(R1+yi)
           yV[counter]=math.sin(t)*(R1+yi)
           rad[counter]=R1+yi
           theta[counter]=np.pi/2-np.arctan2(yV[counter],xV[counter])
           if i==0:
              theta[counter]=0
              xV[counter]=0
           if i==nnt-1:
              theta[counter]=np.pi
              xV[counter]=0
           counter+=1

   if debug:
      np.savetxt('grid.ascii',np.array([xV,yV,theta]).T,header='# x,y')

print("coordinate arrays.........................(%.3fs)" % (timing.time() - start))

###############################################################################
# build iconQ1 array 
###############################################################################
start = timing.time()

iconQ1 =np.zeros((4,nel),dtype=np.int32)

if not axisymmetric:

   counter = 0
   for j in range(0, nelr):
       for i in range(0, nelt):
           icon1=counter
           icon2=counter+1
           icon3=i+(j+1)*nelt+1
           icon4=i+(j+1)*nelt
           if i==nelt-1:
              icon2-=nelt
              icon3-=nelt
           iconQ1[0,counter] = icon2 
           iconQ1[1,counter] = icon1
           iconQ1[2,counter] = icon4
           iconQ1[3,counter] = icon3
           counter += 1
       #end for
   #end for

else:

   counter = 0
   for j in range(0,nelr):
       for i in range(0,nelt):
           iconQ1[0,counter] = i + j * (nelt + 1)
           iconQ1[1,counter] = i + 1 + j * (nelt + 1)
           iconQ1[2,counter] = i + 1 + (j + 1) * (nelt + 1)
           iconQ1[3,counter] = i + (j + 1) * (nelt + 1)
           counter += 1

if debug: export_Q1_mesh_to_vtu(NV,nel,xV,yV,iconQ1)

print("build iconQ1..............................(%.3fs)" % (timing.time() - start))

###############################################################################
# now that the grid has been built as if it was a Q1 grid, 
# we can simply use these same points to arrive at a Q2 
# connectivity array with 4 times less elements.
###############################################################################
start = timing.time()

nelr=nelr//2
nelt=nelt//2
nel=nel//4

if not axisymmetric:
   NP=nelt*(nelr+1)
else:
   NP=(nelt+1)*(nelr+1)

NfemV=NV*ndofV   # Total number of degrees of V freedom 
NfemP=NP*ndofP   # Total number of degrees of P freedom
Nfem=NfemV+NfemP # total number of dofs

h_r=(R2-R1)/nelr # radial resolution

print('  -------------------')
print('  nelr=',nelr,' ; nelt=',nelt,' ; nel=',nel)
print('  NfemV=',NfemV,' ; NfemP=',NfemP, '; Nfem=',Nfem)
print('  nqel=',nqel)
print('  mapping=',mapping)
print('  axisymmetric=',axisymmetric)
print('  exp=',exp)
print('  xi=',xi)
print('  eta_model=',eta_model)
print('  rho_model=',rho_model)
print('  rhoblob=',rhoblob)
print('  Rblob=',Rblob)
print('  zblob=',zblob)
print('  compute_gravity=',compute_gravity)
print('  nel_phi=',nel_phi)
print('  h_r=',h_r)
print('  -------------------')

print("compute nelr,nelt,nel.....................(%.3fs)" % (timing.time() - start))

###############################################################################
# connectivity
###############################################################################
start = timing.time()

iconV =np.zeros((mV,nel),dtype=np.int32,order='F')
iconP =np.zeros((mP,nel),dtype=np.int32,order='F')

if not axisymmetric:

   counter = 0
   for j in range(0, nelr):
       for i in range(0, nelt):
           iconV[0,counter]=2*counter+2 +2*j*nelt
           iconV[1,counter]=2*counter   +2*j*nelt
           iconV[2,counter]=iconV[1,counter]+4*nelt
           iconV[3,counter]=iconV[1,counter]+4*nelt+2
           iconV[4,counter]=iconV[0,counter]-1
           iconV[5,counter]=iconV[1,counter]+2*nelt
           iconV[6,counter]=iconV[2,counter]+1
           iconV[7,counter]=iconV[5,counter]+2
           iconV[8,counter]=iconV[5,counter]+1
           if i==nelt-1:
              iconV[0,counter]-=2*nelt
              iconV[7,counter]-=2*nelt
              iconV[3,counter]-=2*nelt
           #print(j,i,counter,'|',iconV[0:mV,counter])
           counter += 1
       #end for
   #end for

   counter = 0
   for j in range(0, nelr):
       for i in range(0, nelt):
           icon1=counter
           icon2=counter+1
           icon3=i+(j+1)*nelt+1
           icon4=i+(j+1)*nelt
           if i==nelt-1:
              icon2-=nelt
              icon3-=nelt
           iconP[0,counter] = icon2 
           iconP[1,counter] = icon1
           iconP[2,counter] = icon4
           iconP[3,counter] = icon3
           counter += 1
       #end for
   #end for

else:

   counter = 0
   for j in range(0,nelr):
       for i in range(0,nelt):
           iconV[0,counter]=(i)*2+1+(j)*2*nnt -1
           iconV[1,counter]=(i)*2+3+(j)*2*nnt -1
           iconV[2,counter]=(i)*2+3+(j)*2*nnt+nnt*2 -1
           iconV[3,counter]=(i)*2+1+(j)*2*nnt+nnt*2 -1
           iconV[4,counter]=(i)*2+2+(j)*2*nnt -1
           iconV[5,counter]=(i)*2+3+(j)*2*nnt+nnt -1
           iconV[6,counter]=(i)*2+2+(j)*2*nnt+nnt*2 -1
           iconV[7,counter]=(i)*2+1+(j)*2*nnt+nnt -1
           iconV[8,counter]=(i)*2+2+(j)*2*nnt+nnt -1
           counter += 1
       #end for
   #end for

   counter = 0
   for j in range(0,nelr):
       for i in range(0,nelt):
           iconP[0,counter]=i+j*(nelt+1)
           iconP[1,counter]=i+1+j*(nelt+1)
           iconP[2,counter]=i+1+(j+1)*(nelt+1)
           iconP[3,counter]=i+(j+1)*(nelt+1)
           counter += 1
       #end for
   #end for

print("connectivity array........................(%.3fs)" % (timing.time() - start))

###############################################################################

#for iel in range(0,nel):
#    r0=rad[iconV[0,iel]]
#    r3=rad[iconV[3,iel]]
#    r8=rad[iconV[8,iel]]
    #center 0
    #xV[iconV[8,iel]]=0.25*(xV[iconV[0,iel]]+xV[iconV[1,iel]]+xV[iconV[2,iel]]+xV[iconV[3,iel]])
    #yV[iconV[8,iel]]=0.25*(yV[iconV[0,iel]]+yV[iconV[1,iel]]+yV[iconV[2,iel]]+yV[iconV[3,iel]])
    #center 1
    #xV[iconV[8,iel]]=0.125*(xV[iconV[0,iel]]+xV[iconV[1,iel]]+xV[iconV[2,iel]]+xV[iconV[3,iel]]+\
    #                        xV[iconV[4,iel]]+xV[iconV[5,iel]]+xV[iconV[6,iel]]+xV[iconV[7,iel]])
    #yV[iconV[8,iel]]=0.125*(yV[iconV[0,iel]]+yV[iconV[1,iel]]+yV[iconV[2,iel]]+yV[iconV[3,iel]]+\
    #                        yV[iconV[4,iel]]+yV[iconV[5,iel]]+yV[iconV[6,iel]]+yV[iconV[7,iel]])
    #center 2
    #xV[iconV[8,iel]]=0.0625*(xV[iconV[0,iel]]+xV[iconV[1,iel]]+xV[iconV[2,iel]]+xV[iconV[3,iel]]+\
    #                         3*xV[iconV[4,iel]]+3*xV[iconV[5,iel]]+3*xV[iconV[6,iel]]+3*xV[iconV[7,iel]])
    #yV[iconV[8,iel]]=0.0625*(yV[iconV[0,iel]]+yV[iconV[1,iel]]+yV[iconV[2,iel]]+yV[iconV[3,iel]]+\
    #                         3*yV[iconV[4,iel]]+3*yV[iconV[5,iel]]+3*yV[iconV[6,iel]]+3*yV[iconV[7,iel]])

###############################################################################
#now that I have both connectivity arrays I can easily build xP,yP (Q1 space)
###############################################################################
start = timing.time()

xP=np.empty(NP,dtype=np.float64)  # x coordinates
yP=np.empty(NP,dtype=np.float64)  # y coordinates

for iel in range(0,nel):
    xP[iconP[0,iel]]=xV[iconV[0,iel]]
    xP[iconP[1,iel]]=xV[iconV[1,iel]]
    xP[iconP[2,iel]]=xV[iconV[2,iel]]
    xP[iconP[3,iel]]=xV[iconV[3,iel]]
    yP[iconP[0,iel]]=yV[iconV[0,iel]]
    yP[iconP[1,iel]]=yV[iconV[1,iel]]
    yP[iconP[2,iel]]=yV[iconV[2,iel]]
    yP[iconP[3,iel]]=yV[iconV[3,iel]]

print("compute xP,yP.............................(%.3fs)" % (timing.time() - start))

###############################################################################
# find out nodes on hull
###############################################################################
start = timing.time()

hull=np.zeros(NV,dtype=bool) 

for i in range(0,NV):
    if rad[i]/R1<1+eps:
       hull[i]=True
    if rad[i]/R2>1-eps:
       hull[i]=True
    if axisymmetric and xV[i]/R2<eps:
       hull[i]=True

print("flag nodes on hull........................(%.3fs)" % (timing.time() - start))

###############################################################################
# flag surface nodes and elements
###############################################################################
start = timing.time()

surfaceV=np.zeros(NV,dtype=bool) 
surface_element=np.zeros(nel,dtype=bool) 
cmbV=np.zeros(NV,dtype=bool) 
cmb_element=np.zeros(nel,dtype=bool) 

for i in range(0,NV):
    if rad[i]/R2>1-eps:
       surfaceV[i]=True
    if rad[i]/R1<1+eps:
       cmbV[i]=True

for iel in range(0,nel):
    if surfaceV[iconV[2,iel]]:
       surface_element[iel]=True 
    if cmbV[iconV[0,iel]]:
       cmb_element[iel]=True 

print("flag surf and cmb nodes+elts..............(%.3fs)" % (timing.time() - start))

###############################################################################
# compute normal vectors
###############################################################################
start = timing.time()

nx,ny=compute_normals(normal_type,NV,xV,yV,mV,nqel,surfaceV,hull,\
                      R1,R2,eps,iconV,axisymmetric,rad,theta,\
                      qcoords_r,qcoords_s,qweights,)

if debug:
   np.savetxt('normals.ascii',np.array([xV[surfaceV],yV[surfaceV],\
                                        nx[surfaceV],ny[surfaceV],theta[surfaceV]]).T)

print("compute surface normals...................(%.3fs)" % (timing.time() - start))

###############################################################################
# define boundary conditions
###############################################################################
start = timing.time()

bc_fix = np.zeros(Nfem,dtype=bool)  
bc_val = np.zeros(Nfem,dtype=np.float64) 

for i in range(0,NV):
    if rad[i]/R1<1+eps:
       bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = velocity_x(xV[i],yV[i],R1,R2,exp)
       bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = velocity_y(xV[i],yV[i],R1,R2,exp)
    if rad[i]/R2>(1-eps) and not surface_free_slip:
       bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = velocity_x(xV[i],yV[i],R1,R2,exp)
       bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = velocity_y(xV[i],yV[i],R1,R2,exp)

    #vertical wall x=0
    if axisymmetric and xV[i]/R1<eps:
       bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = 0 #u=0
       if rad[i]/R2>1-eps:
          bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0 #also v=0 for 2 pts
          #bc_fix[i*ndofV] = False 

print("defining boundary conditions..............(%.3fs)" % (timing.time() - start))

###############################################################################
# compute array for assembly
###############################################################################
start = timing.time()

ndofV_el=mV*ndofV
local_to_globalV=np.zeros((ndofV_el,nel),dtype=np.int32)

for iel in range(0,nel):
    for k1 in range(0,mV):
        for i1 in range(0,ndofV):
            ikk=ndofV*k1          +i1
            m1 =ndofV*iconV[k1,iel]+i1
            local_to_globalV[ikk,iel]=m1
    
print("compute local_to_global...................(%.3fs)" % (timing.time() - start))

###############################################################################
# fill I,J arrays
###############################################################################
start = timing.time()

bignb=nel*( (mV*ndofV)**2 + 2*(mV*ndofV*mP) )

I=np.zeros(bignb,dtype=np.int32)
J=np.zeros(bignb,dtype=np.int32)
V=np.zeros(bignb,dtype=np.float64)

counter=0
for iel in range(0,nel):
    for ikk in range(ndofV_el):
        m1=local_to_globalV[ikk,iel]
        for jkk in range(ndofV_el):
            m2=local_to_globalV[jkk,iel]
            I[counter]=m1
            J[counter]=m2
            counter+=1
        for jkk in range(0,mP):
            m2 =iconP[jkk,iel]+NfemV
            I[counter]=m1
            J[counter]=m2
            counter+=1
            I[counter]=m2
            J[counter]=m1
            counter+=1

print("fill I,J arrays...........................(%.3fs)" % (timing.time() - start))

##############################################################################################
##############################################################################################
##############################################################################################
# time stepping
##############################################################################################
##############################################################################################
##############################################################################################

for istep in range(0,nstep):

    print('**************************************************')
    print('******************** istep='+str(istep)+' *********************')
    print('**************************************************')

    ###############################################################################
    ###############################################################################
    start = timing.time()

    xmapping,zmapping=define_mapping(mapping,mmapping,xV,yV,iconV,nel,axisymmetric,rad,theta)

    print(spacing+" -> xmapping (m,M) %.2e %.2e " %(np.min(xmapping),np.max(xmapping)))
    print(spacing+" -> zmapping (m,M) %.2e %.2e " %(np.min(zmapping),np.max(zmapping)))

    if debug:
       np.savetxt('xzmapping'+mapping+'.ascii',np.array([xmapping[0,:],zmapping[0,:]]).T)
       export_mapping_points_to_vtu(mapping,mmapping,xmapping,zmapping)
       export_quadrature_points_to_vtu(nqperdim,nqel,qcoords_r,qcoords_s,mapping,\
                                       xmapping,zmapping)

    print("define mapping............................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute element center coordinates
    ###############################################################################
    start = timing.time()

    xc,zc,thetac=compute_element_center_coords(nel,mapping,xmapping,zmapping)

    print(spacing+" -> xc (m,M) %.2e %.2e " %(np.min(xc),np.max(xc)))
    print(spacing+" -> zc (m,M) %.2e %.2e " %(np.min(zc),np.max(zc)))

    print("compute center coords.....................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute elemental and nodal viscosity+density fields for 1st timestep only
    ###############################################################################
    start = timing.time()

    if istep==0:

       density_elemental,density_nodal,viscosity_elemental,viscosity_nodal=\
       compute_rho_eta_fields(mV,NV,nel,xc,zc,iconV,R1,R2,rho_m,eta_m,eta_model,\
                              rho_model,exp,rhoblobstar,zblob,Rblob)

       print(spacing+" -> eta_e (m,M) %.2e %.2e " %(np.min(viscosity_elemental),np.max(viscosity_elemental)))
       print(spacing+" -> rho_e (m,M) %.2e %.2e " %(np.min(density_elemental),np.max(density_elemental)))
       print(spacing+" -> eta_n (m,M) %.2e %.2e " %(np.min(viscosity_nodal),np.max(viscosity_nodal)))
       print(spacing+" -> rho_n (m,M) %.2e %.2e " %(np.min(density_nodal),np.max(density_nodal)))

    if debug:
       np.savetxt('xycenter'+mapping+'.ascii',np.array([xc,yc,density_elemental,viscosity_elemental]).T)

    print("compute rho+eta fields....................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute area and volume of elements
    # note that for the volume calculation we only consider the elements x>0
    ###############################################################################
    start = timing.time()

    area=np.zeros(nel,dtype=np.float64) 
    vol=np.zeros(nel,dtype=np.float64) 
    massq  = np.zeros(nel*nqel,dtype=np.float64) 
    radq  = np.zeros(nel*nqel,dtype=np.float64) 
    thetaq  = np.zeros(nel*nqel,dtype=np.float64) 
    zq  = np.zeros(nel*nqel,dtype=np.float64) 
    rhoq = np.zeros(nel*nqel,dtype=np.float64) 
    etaq = np.zeros(nel*nqel,dtype=np.float64) 

    counterq=0
    jcb=np.zeros((2,2),dtype=np.float64)
    for iel in range(0,nel):
        for kq in range(0,nqel):
            rq=qcoords_r[kq]
            sq=qcoords_s[kq]
            weightq=qweights[kq]
            NNNV=NNN(rq,sq,mapping)
            dNNNVdr=dNNNdr(rq,sq,mapping)
            dNNNVds=dNNNds(rq,sq,mapping)
            xq=np.dot(NNNV[:],xmapping[:,iel])
            yq=np.dot(NNNV[:],zmapping[:,iel])
            zq[counterq]=yq
            radq[counterq]=np.sqrt(xq**2+yq**2)
            thetaq[counterq]=np.pi/2-np.arctan2(yq,xq)
            jcb[0,0]=np.dot(dNNNVdr[:],xmapping[:,iel])
            jcb[0,1]=np.dot(dNNNVdr[:],zmapping[:,iel])
            jcb[1,0]=np.dot(dNNNVds[:],xmapping[:,iel])
            jcb[1,1]=np.dot(dNNNVds[:],zmapping[:,iel])
            jcob = np.linalg.det(jcb)
            area[iel]+=jcob*weightq
            massq[counterq]=density_elemental[iel]*weightq*jcob*2*np.pi*xq
            if xq>0: vol[iel]+=jcob*weightq*2*np.pi*xq

            if use_elemental_density:
               rhoq[counterq]=density_elemental[iel]
            else:
               rhoq[counterq]=density(xq,yq,R1,R2,rho_m,rho_model,\
                                      exp,rhoblobstar,zblob,Rblob)

            if use_elemental_viscosity:
               etaq[counterq]=viscosity_elemental[iel]
            else:
               etaq[counterq]=viscosity(xq,yq,R1,R2,eta_m,eta_model)

            counterq+=1
       #end for
    #end for

    print(spacing+" -> area (m,M) %.6e %.6e " %(np.min(area),np.max(area)))
    print(spacing+" -> total area (meas) %.12e | nel= %d" %(area.sum(),nel))
    if not axisymmetric:
       print(spacing+" -> total area (anal) %e " %(np.pi*(R2**2-R1**2)))
    else:
       print(spacing+" -> total area (anal) %e " %(np.pi*(R2**2-R1**2)/2))

    print(spacing+" -> total volume (meas) %.12e | nel= %d" %(vol.sum(),nel))
    print(spacing+" -> total volume (anal) %.12e" %(4*np.pi/3*(R2**3-R1**3)))

    print(spacing+" -> rhoq (m,M) %.3e %.3e " %(np.min(rhoq),np.max(rhoq)))
    print(spacing+" -> etaq (m,M) %.3e %.3e " %(np.min(etaq),np.max(etaq)))
    print(spacing+" -> massq (m,M) %.3e %.3e " %(np.min(massq),np.max(massq)))
    print(spacing+" -> total mass (meas) %.12e | nel= %d" %(np.sum(density_elemental*vol),nel))
    
    print("sanity check..............................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # build FE matrix
    ###############################################################################
    start = timing.time()

    f_rhs = np.zeros(NfemV,dtype=np.float64) 
    h_rhs = np.zeros(NfemP,dtype=np.float64) 
    dNNNVdx  = np.zeros(mV,dtype=np.float64) 
    dNNNVdy  = np.zeros(mV,dtype=np.float64) 
    jcb=np.zeros((2,2),dtype=np.float64)

    if axisymmetric:
       #c_mat=np.array([[2,0,0,0],[0,2,0,0],[0,0,2,0],[0,0,0,1]],dtype=np.float64) 
       c_mat=np.array([[4/3,-2/3,-2/3,0],[-2/3,4/3,-2/3,0],[-2/3,-2/3,4/3,0],[0,0,0,1]],dtype=np.float64) 
       b_mat= np.zeros((4,ndofV*mV),dtype=np.float64) # gradient matrix B 
       N_mat= np.zeros((4,ndofP*mP),dtype=np.float64) # matrix  
    else:
       c_mat=np.array([[2,0,0],[0,2,0],[0,0,1]],dtype=np.float64) 
       b_mat= np.zeros((3,ndofV*mV),dtype=np.float64) # gradient matrix B 
       N_mat= np.zeros((3,ndofP*mP),dtype=np.float64) # matrix  

    counter=0
    counterq=0
    for iel in range(0,nel):

        if not solve_stokes: break 

        f_el =np.zeros((mV*ndofV),dtype=np.float64)
        K_el =np.zeros((mV*ndofV,mV*ndofV),dtype=np.float64)
        G_el=np.zeros((mV*ndofV,mP*ndofP),dtype=np.float64)
        h_el=np.zeros((mP*ndofP),dtype=np.float64)

        for kq in range(0,nqel):
            rq=qcoords_r[kq]
            sq=qcoords_s[kq]
            weightq=qweights[kq]

            #compute coords of quadrature points
            NNNV=NNN(rq,sq,mapping)
            xq=np.dot(NNNV[:],xmapping[:,iel])
            yq=np.dot(NNNV[:],zmapping[:,iel])
            #zq[counterq]=yq
            #radq[counterq]=np.sqrt(xq**2+yq**2)
            #thetaq[counterq]=np.pi/2-np.arctan2(yq,xq)

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
            NNNV=NNN(rq,sq,'Q2')
            dNNNVdr=dNNNdr(rq,sq,'Q2')
            dNNNVds=dNNNds(rq,sq,'Q2')
            NNNP=NNN(rq,sq,'Q1')

            # compute dNdx & dNdy
            for k in range(0,mV):
                dNNNVdx[k]=jcbi[0,0]*dNNNVdr[k]+jcbi[0,1]*dNNNVds[k]
                dNNNVdy[k]=jcbi[1,0]*dNNNVdr[k]+jcbi[1,1]*dNNNVds[k]
            #end for 

            if axisymmetric:

               coeffq=weightq*jcob*2*np.pi*xq

               for i in range(0,mV):
                   b_mat[0:4,2*i:2*i+2] = [[dNNNVdx[i],0.       ],
                                           [NNNV[i]/xq,0.       ],
                                           [0.        ,dNNNVdy[i]],
                                           [dNNNVdy[i],dNNNVdx[i]]]
               for i in range(0,mP):
                   N_mat[0,i]=NNNP[i]
                   N_mat[1,i]=NNNP[i]
                   N_mat[2,i]=NNNP[i]
            else:

               coeffq=weightq*jcob

               for i in range(0,mV):
                   b_mat[0:3,2*i:2*i+2] = [[dNNNVdx[i],0.     ],
                                           [0.        ,dNNNVdy[i]],
                                           [dNNNVdy[i],dNNNVdx[i]]]
               for i in range(0,mP):
                   N_mat[0,i]=NNNP[i]
                   N_mat[1,i]=NNNP[i]

            K_el+=b_mat.T.dot(c_mat.dot(b_mat))*etaq[counterq]*coeffq

            G_el-=b_mat.T.dot(N_mat)*coeffq

            g_x,g_y=gravity_acceleration(xq,yq,R1,R2,gravity_model,g0,rho_m,\
                                         rhoc,rhoblob,Rblob,zblob)
            for i in range(0,mV):
                f_el[ndofV*i  ]+=NNNV[i]*coeffq*g_x*rhoq[counterq]
                f_el[ndofV*i+1]+=NNNV[i]*coeffq*g_y*rhoq[counterq]
            #end for 

            counterq+=1

        #end for kq

        G_el*=eta_ref/h_r

        #Kmax=np.max(abs(K_el))
        #Gmax=np.max(abs(G_el))
        #fmax=np.max(abs(f_el))
        #K_el/=Kmax
        #G_el/=Gmax
        #f_el/=fmax
        #for i in range(0,mV): #counter rotation by angle alpha
        #    RotMat[2*i  ,2*i]= np.cos(alpha) ; RotMat[2*i  ,2*i+1]=np.sin(alpha)
        #    RotMat[2*i+1,2*i]=-np.sin(alpha) ; RotMat[2*i+1,2*i+1]=np.cos(alpha)

        if surface_free_slip and surface_element[iel]==1:
           for k in range(0,mV):
               inode=iconV[k,iel]
               #if surfaceV[inode] and xV[inode]>0 and (not bc_fix[inode*ndofV]):
               if surfaceV[inode] and (not bc_fix[inode*ndofV]):

                  #print(xV[inode],yV[inode],np.arctan2(nx[inode],ny[inode]),theta[inode])
                  #alpha=-np.arctan2(nx[inode],ny[inode])

                  #if theta[inode]<np.pi/4:
                  #alpha=-theta[inode]
                  #o=1 #y-component set to 0
                  #else:
                  alpha=np.pi/2-theta[inode]
                  o=0 #x-component set to 0

                  RotMat=np.zeros((mV*ndofV,mV*ndofV),dtype=np.float64)
                  for i in range(0,mV*ndofV):
                      RotMat[i,i]=1.
                  RotMat[2*k  ,2*k]= np.cos(alpha) ; RotMat[2*k  ,2*k+1]=np.sin(alpha)
                  RotMat[2*k+1,2*k]=-np.sin(alpha) ; RotMat[2*k+1,2*k+1]=np.cos(alpha)

                  # apply counter rotation
                  K_el=RotMat.dot(K_el.dot(RotMat.T))
                  f_el=RotMat.dot(f_el)
                  G_el=RotMat.dot(G_el)

                  # apply boundary conditions,
                  ikk=ndofV*k+o
                  K_ref=K_el[ikk,ikk]
                  for jkk in range(0,mV*ndofV):
                      K_el[ikk,jkk]=0
                      K_el[jkk,ikk]=0
                  K_el[ikk,ikk]=K_ref
                  f_el[ikk]=0#K_ref*bc_val[m1]
                  #h_el[:]-=G_el[ikk,:]*bc_val[m1]
                  G_el[ikk,:]=0

                  # rotate back
                  K_el=RotMat.T.dot(K_el.dot(RotMat))
                  f_el=RotMat.T.dot(f_el)
                  G_el=RotMat.T.dot(G_el)

               #end if
           #end for

        # impose Dirichlet b.c. 
        for k1 in range(0,mV):
            for i1 in range(0,ndofV):
                ikk=ndofV*k1+i1
                m1 =ndofV*iconV[k1,iel]+i1
                if bc_fix[m1]:
                   K_ref=K_el[ikk,ikk] 
                   for jkk in range(0,mV*ndofV):
                       f_el[jkk]-=K_el[jkk,ikk]*bc_val[m1]
                       K_el[ikk,jkk]=0
                       K_el[jkk,ikk]=0
                   #end for 
                   K_el[ikk,ikk]=K_ref
                   f_el[ikk]=K_ref*bc_val[m1]
                   h_el[:]-=G_el[ikk,:]*bc_val[m1]
                   G_el[ikk,:]=0
                #end if 
            #end for 
        #end for 

        # assemble matrix K_mat and right hand side rhs
        for ikk in range(ndofV_el):
            m1=local_to_globalV[ikk,iel]
            for jkk in range(ndofV_el):
                V[counter]=K_el[ikk,jkk]
                counter+=1
            for jkk in range(0,mP):
                V[counter]=G_el[ikk,jkk]
                counter+=1
                V[counter]=G_el[ikk,jkk]
                counter+=1
            f_rhs[m1]+=f_el[ikk]
        for k2 in range(0,mP):
            m2=iconP[k2,iel]
            h_rhs[m2]+=h_el[k2]

    #end for iel


    print("build FE matrixs & rhs....................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # solve system
    ###############################################################################
    start = timing.time()

    if solve_stokes:

       rhs=np.zeros(Nfem,dtype=np.float64)
       rhs[0:NfemV]=f_rhs
       rhs[NfemV:NfemV+NfemP]=h_rhs
       sparse_matrix = sparse.coo_matrix((V,(I,J)),shape=(Nfem,Nfem)).tocsr()
       sol=sps.linalg.spsolve(sparse_matrix,rhs)

    else:

       sol=np.zeros(Nfem,dtype=np.float64)  

    print("solving system............................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # put solution into separate x,y velocity arrays
    ###############################################################################
    start = timing.time()

    u,v=np.reshape(sol[0:NfemV],(NV,2)).T
    p=sol[NfemV:NfemV+NfemP]*eta_ref/h_r


    if debug:
       np.savetxt('velocity.ascii',np.array([xV,yV,u/vel_unit,v/vel_unit]).T,header='# x,y,u,v')
 
    vr= np.sin(theta)*u+np.cos(theta)*v
    vt= np.cos(theta)*u-np.sin(theta)*v
   
    vel=np.zeros(NV,dtype=np.float64)  
    vel[:]=np.sqrt(u[:]**2+v[:]**2)

    print(spacing+" -> nelr= %d | u   (m,M) %.7e %.7e " %(nelr,np.min(u)/vel_unit,np.max(u)/vel_unit),velunit)
    print(spacing+" -> nelr= %d | v   (m,M) %.7e %.7e " %(nelr,np.min(v)/vel_unit,np.max(v)/vel_unit),velunit)
    print(spacing+" -> nelr= %d | v_r (m,M) %.7e %.7e " %(nelr,np.min(vr)/vel_unit,np.max(vr)/vel_unit),velunit)
    print(spacing+" -> nelr= %d | v_t (m,M) %.7e %.7e " %(nelr,np.min(vt)/vel_unit,np.max(vt)/vel_unit),velunit)
    print(spacing+" -> nelr= %d | vel (m,M) %.7e %.7e " %(nelr,np.min(vel)/vel_unit,np.max(vel)/vel_unit),velunit)

    print("reshape solution..........................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute elemental pressure
    ###############################################################################
    start = timing.time()

    pc=np.zeros(nel,dtype=np.float64)

    for iel in range(0,nel):
        pc[iel]=np.sum(p[iconP[0:4,iel]])/4

    print("compute elemental pressure................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute strain rate - center to nodes - method 1
    ###############################################################################
    start = timing.time()
   
    exxc=np.zeros(nel,dtype=np.float64)
    eyyc=np.zeros(nel,dtype=np.float64)
    exyc=np.zeros(nel,dtype=np.float64)

    exx1=np.zeros(NV,dtype=np.float64)  
    exy1=np.zeros(NV,dtype=np.float64)  
    eyy1=np.zeros(NV,dtype=np.float64)  

    if compute_sr1:

       count=np.zeros(NV,dtype=np.int32)  

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

       print(spacing+"     -> exx1 (m,M) %e %e " %(np.min(exx1),np.max(exx1)))
       print(spacing+"     -> eyy1 (m,M) %e %e " %(np.min(eyy1),np.max(eyy1)))
       print(spacing+"     -> exy1 (m,M) %e %e " %(np.min(exy1),np.max(exy1)))

    #end if

    print("compute strain rate meth-1................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute strain rate - corners to nodes - method 2
    ###############################################################################
    start = timing.time()

    exx2,eyy2,exy2=compute_strain_rate2(nel,mV,NV,iconV,mapping,xmapping,zmapping,u,v)

    print(spacing+" -> exx2 (m,M) %e %e " %(np.min(exx2),np.max(exx2)))
    print(spacing+" -> eyy2 (m,M) %e %e " %(np.min(eyy2),np.max(eyy2)))
    print(spacing+" -> exy2 (m,M) %e %e " %(np.min(exy2),np.max(exy2)))

    sr2=np.zeros(NV,dtype=np.float64)  
    sr2[:]=np.sqrt(0.5*(exx2[:]**2+eyy2[:]**2)+exy2[:]**2)

    if debug:
       np.savetxt('strainrate'+str(istep)+'.ascii',np.array([xV,yV,exx2,eyy2,exy2]).T)

    print("compute strain rate meth-2................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # project pressure onto Q2 mesh
    ###############################################################################
    start = timing.time()

    q=project_pressure_on_Q2(NV,nel,mV,mP,p,iconP,iconV)

    print(spacing+" -> q (m,M) %e %e " %(np.min(q),np.max(q)))

    if debug:
       np.savetxt('pressure'+str(istep)+'.ascii',np.array([xV,yV,q]).T)

    print("project pressure onto Q2 mesh.............(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # convert strain rate tensor to spherical coordinates
    # Note that there is a pb when we run the model in plane strain. In that 
    # case I set the err,ert,ett to zero for x<0
    ###############################################################################
    start = timing.time()

    e_rr2=np.zeros(NV,dtype=np.float64)  
    e_tt2=np.zeros(NV,dtype=np.float64)  
    e_rt2=np.zeros(NV,dtype=np.float64)  

    for i in range(0,NV):
        if xV[i]>=0:
           e_rr2[i]=exx2[i]*np.sin(theta[i])**2+\
                   2*exy2[i]*np.sin(theta[i])*np.cos(theta[i])+\
                   eyy2[i]*np.cos(theta[i])**2
           e_tt2[i]=exx2[i]*np.cos(theta[i])**2-\
                   2*exy2[i]*np.sin(theta[i])*np.cos(theta[i])+\
                   eyy2[i]*np.sin(theta[i])**2
           e_rt2[i]=(exx2[i]-eyy2[i])*np.sin(theta[i])*np.cos(theta[i])+\
                   exy2[i]*(-np.sin(theta[i])**2+\
                   np.cos(theta[i])**2)

    e_rrc=np.zeros(nel,dtype=np.float64)  
    for iel in range(0,nel):
        if xc[iel]>=0:
           e_rrc[iel]=exxc[iel]*np.sin(thetac[iel])**2+\
                      2*exyc[iel]*np.sin(thetac[iel])*np.cos(thetac[iel])+\
                      eyyc[iel]*np.cos(thetac[iel])**2

    print(spacing+" -> e_rr (m,M) %e %e | nel= %d" %(np.min(e_rr2),np.max(e_rr2),nel))
    print(spacing+" -> e_tt (m,M) %e %e | nel= %d" %(np.min(e_tt2),np.max(e_tt2),nel))
    print(spacing+" -> e_rt (m,M) %e %e | nel= %d" %(np.min(e_rt2),np.max(e_rt2),nel))

    print("compute strain rate in sph. coords........(%.3fs)" % (timing.time() - start))

    ###############################################################################
    start = timing.time()

    exx3=np.zeros(NV,dtype=np.float64)  
    eyy3=np.zeros(NV,dtype=np.float64)  
    exy3=np.zeros(NV,dtype=np.float64)  

    if compute_sr3:

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

       print(spacing+"     -> Lxx3 (m,M) %e %e " %(np.min(Lxx3),np.max(Lxx3)))
       print(spacing+"     -> Lyy3 (m,M) %e %e " %(np.min(Lyy3),np.max(Lyy3)))
       print(spacing+"     -> Lxy3 (m,M) %e %e " %(np.min(Lxy3),np.max(Lxy3)))
       print(spacing+"     -> Lxy3 (m,M) %e %e " %(np.min(Lyx3),np.max(Lyx3)))

       exx3[:]=Lxx3[:]
       eyy3[:]=Lyy3[:]
       exy3[:]=0.5*(Lxy3[:]+Lyx3[:])

    print("compute strain rate meth-3................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # normalise pressure
    # note that the integration could be improved (I currently only sample the 
    # pressure in 1 point in the middle of the edge)
    ###############################################################################
    start = timing.time()

    if not axisymmetric:
       if exp==0:
          poffset=np.sum(q[0:2*nelt])/(2*nelt)
       else:
          poffset=0
          for iel in range(0,nel):
              if surface_element[iel]:
                 dtheta=2*np.pi/nelt 
                 pmean=0.5*(p[iconP[2,iel]]+p[iconP[3,iel]])
                 poffset+=dtheta*pmean
          poffset/=2*np.pi

       q-=poffset
       p-=poffset

    else: #zero average pressure on surface

       poffset=0
       for iel in range(0,nel):
           if surface_element[iel]:
              dtheta=theta[iconV[2,iel]]-theta[iconV[3,iel]]
              pmean=0.5*(p[iconP[2,iel]]+p[iconP[3,iel]])
              poffset+=np.sin((theta[iconV[2,iel]]+theta[iconV[3,iel]])/2)*dtheta\
                       *2*np.pi*R2**2 * pmean
       poffset/=4*np.pi*R2**2
       q-=poffset
       p-=poffset

       poffset=0
       for iel in range(0,nel):
           if surface_element[iel]:
              dtheta=theta[iconV[2,iel]]-theta[iconV[3,iel]]
              poffset+=np.sin((theta[iconV[2,iel]]+theta[iconV[3,iel]])/2)*dtheta * 2*np.pi*R2**2 * pc[iel]
       poffset/=4*np.pi*R2**2
       pc-=poffset

    print(spacing+" -> p (m,M) %e %e " %(np.min(p),np.max(p)))
    print(spacing+" -> q (m,M) %e %e " %(np.min(q),np.max(q)))

    if debug: np.savetxt('pressure.ascii',np.array([xP,yP,p]).T,header='# x,y,p')

    print("normalise pressure........................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute nodal error fields (for plotting)
    ###############################################################################
    start = timing.time()

    u_err = np.zeros(NV,dtype=np.float64) 
    v_err = np.zeros(NV,dtype=np.float64)    
    p_err = np.zeros(NP,dtype=np.float64)    

    if exp==0:
       for i in range(0,NV):
           u_err[i]=u[i]-velocity_x(xV[i],yV[i],R1,R2,exp)
           v_err[i]=v[i]-velocity_y(xV[i],yV[i],R1,R2,exp)

       for i in range(0,NP):
           p_err[i]=p[i]-pressure(xP[i],yP[i],R1,R2,rho_m,g0,exp)

       print(spacing+" -> u_err (m,M) %.10e %.10e | nelr= %d" %(np.min(u_err),np.max(u_err),nelr))
       print(spacing+" -> v_err (m,M) %.10e %.10e | nelr= %d" %(np.min(v_err),np.max(v_err),nelr))
       print(spacing+" -> p_err (m,M) %.10e %.10e | nelr= %d" %(np.min(p_err),np.max(p_err),nelr))

    print("compute error fields......................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # since I have zeroed the the spherical components of the strain rate for x<0
    # then I also zero the dynamic topography
    # note that this is only valid for constant viscosity!
    ###############################################################################
    start = timing.time()

    dyn_topo_nodal=np.zeros(NV,dtype=np.float64)
    if exp>0:
       for i in range(0,NV):
           if surfaceV[i] and xV[i]>=0:
              dyn_topo_nodal[i]= -(2*viscosity_nodal[i]*e_rr2[i]-q[i])/(rho_m*g0) 

    dyn_topo_eltal=np.zeros(nel,dtype=np.float64)
    if exp>0:
       for iel in range(0,nel):
           if surface_element[iel] and xc[iel]>=0:
              dyn_topo_eltal[iel]= -(2*viscosity_elemental[iel]*e_rrc[iel]-pc[iel])/(rho_m*g0) 

    print("compute dynamic topography................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # export fields at both surfaces
    ###############################################################################
    start = timing.time()

    surfaceQ1=np.zeros(NV,dtype=bool) 
    surfaceQ2=np.zeros(NV,dtype=bool) 
    for iel in range(0,nel):
        if surface_element[iel]:
           surfaceQ1[iconV[2,iel]]=True
           surfaceQ1[iconV[3,iel]]=True
           surfaceQ2[iconV[2,iel]]=True
           surfaceQ2[iconV[3,iel]]=True
           surfaceQ2[iconV[6,iel]]=True

    leftV=np.zeros(NV,dtype=bool) 
    for i in range(0,NV):
        if xV[i]<eps and yV[i]>0: leftV[i]=True

    #src=np.zeros(nel,dtype=np.float64)  
    #src[:]=np.sqrt(0.5*(exxc[:]**2+eyyc[:]**2)+exyc[:]**2)
    #np.savetxt('pc_R2.ascii',np.array([thetac[surface_element],pc[surface_element]]).T)
    #np.savetxt('src_R2.ascii',np.array([thetac[surface_element],src[surface_element]]).T)
    #np.savetxt('errc_R2.ascii',np.array([thetac[surface_element],e_rrc[surface_element]]).T)
    #np.savetxt('d_tc_R2.ascii',np.array([thetac[surface_element],dyn_topo_eltal[surface_element]]).T)

    #np.savetxt('qqq_R1.ascii',np.array([xV[0:nnt],yV[0:nnt],q[0:nnt],theta[0:nnt]]).T)
    #np.savetxt('sr2_R1.ascii',np.array([xV[0:nnt],yV[0:nnt],sr2[0:nnt],theta[0:nnt]]).T)
    #np.savetxt('vel_R1.ascii',np.array([xV[0:nnt],yV[0:nnt],vel[0:nnt],theta[0:nnt]]).T)

    np.savetxt('qqq_R2_'+str(istep)+'.ascii',np.array([theta[surfaceQ1],q[surfaceQ1]]).T)
    np.savetxt('sr2_R2_'+str(istep)+'.ascii',np.array([theta[surfaceQ2],sr2[surfaceQ2]]).T)
    np.savetxt('vel_R2_'+str(istep)+'.ascii',np.array([theta[surfaceQ2],vel[surfaceQ2],vr[surfaceQ2],vt[surfaceQ2]]).T)
    np.savetxt('v_r_R2_'+str(istep)+'.ascii',np.array([theta[surfaceQ2],vr[surfaceQ2]]).T)
    np.savetxt('v_t_R2_'+str(istep)+'.ascii',np.array([theta[surfaceQ2],vt[surfaceQ2]]).T)
    np.savetxt('err_R2_'+str(istep)+'.ascii',np.array([theta[surfaceQ2],e_rr2[surfaceQ2]]).T)
    np.savetxt('d_t_R2_'+str(istep)+'.ascii',np.array([theta[surfaceQ2],dyn_topo_nodal[surfaceQ2]]).T)
    np.savetxt('sr2_left_'+str(istep)+'.ascii',np.array([xV[leftV],yV[leftV],sr2[leftV],rad[leftV]]).T)
    np.savetxt('vel_left_'+str(istep)+'.ascii',np.array([xV[leftV],yV[leftV],vel[leftV],rad[leftV]]).T)

    if compute_sr1:
       sr1=np.zeros(NV,dtype=np.float64)  
       sr1[:]=np.sqrt(0.5*(exx1[:]**2+eyy1[:]**2)+exy1[:]**2)
       np.savetxt('sr1_R1.ascii',np.array([xV[0:nnt],yV[0:nnt],sr1[0:nnt],theta[0:nnt]]).T)
       np.savetxt('sr1_R2.ascii',np.array([xV[surfaceQ2],yV[surfaceQ2],sr1[surfaceQ2],theta[surfaceQ2]]).T)

    if compute_sr3:
       sr3=np.zeros(NV,dtype=np.float64)  
       sr3[:]=np.sqrt(0.5*(exx3[:]**2+eyy3[:]**2)+exy3[:]**2)
       np.savetxt('sr3_R1.ascii',np.array([xV[0:nnt],yV[0:nnt],sr3[0:nnt],theta[0:nnt]]).T)
       np.savetxt('sr3_R2.ascii',np.array([xV[surfaceQ2],yV[surfaceQ2],sr3[surfaceQ2],theta[surfaceQ2]]).T)

    print("export p&q on R1,R2.......................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute error
    ###############################################################################
    start = timing.time()

    errv,errp,vrms=compute_errors(nel,nqel,mapping,xmapping,zmapping,qcoords_r,\
                                  qcoords_s,qweights,R1,R2,rho_m,g0,exp,u,v,p,\
                                  axisymmetric,iconV,iconP)

    print(spacing+' -> nelr= %d ; vrms= %.12e' %(nelr,vrms/vel_unit))
    if exp==0:
       print(spacing+" -> nelr= %d ; errv= %.12e ; errp= %.14e " %(nelr,errv,errp))

    print("compute errors............................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute gravity acceleration on pressure nodes
    ###############################################################################
    start = timing.time()

    if self_gravitation:

       g_x=np.zeros(NP,dtype=np.float64)   
       g_y=np.zeros(NP,dtype=np.float64)   
       g_z=np.zeros(NP,dtype=np.float64)   
       gnorm=np.zeros(NP,dtype=np.float64)   

       for i in range(0,NP):
           if i%100==0: print('i=',i,'/',NP)
           g_x[i],g_y[i],g_z[i]=\
           compute_gravity_at_point(xP[i],yP[i],nel,nqel,zq,radq,thetaq,massq,nel_phi)
           #print(xV[i],yV[i],g_x[i],g_y[i],g_z[i])

       gnorm[:]=np.sqrt(g_x[:]**2+g_y[:]**2+g_z[:]**2)

       np.savetxt('gravity_self.ascii',np.array([xP,yP,g_x,g_y,g_z,gnorm]).T)

       print(spacing+"     -> g_x (m,M) %.4e %.4e " %(np.min(g_x),np.max(g_x)))
       print(spacing+"     -> g_y (m,M) %.4e %.4e " %(np.min(g_y),np.max(g_y)))
       print(spacing+"     -> g_z (m,M) %.4e %.4e " %(np.min(g_z),np.max(g_z)))

       print("self gravity (%.3fs)" % (timing.time() - start))

    ###############################################################################
    # plot of solution
    ###############################################################################
    start = timing.time()

    if visu==1:

       export_solutionQ2_to_vtu(istep,NV,nel,xV,yV,vel_unit,u,v,viscosity_nodal,nx,ny,\
                                hull,surfaceV,cmbV,bc_fix,area,viscosity_elemental,\
                                density_elemental,surface_element,cmb_element,iconV)
       export_solutionQ1_to_vtu(istep,exp,NV,nel,xV,yV,vel_unit,u,v,vr,vt,rad,theta,\
                                exx2,eyy2,exy2,sr2,e_rr2,e_tt2,e_rt2,q,viscosity_nodal,\
                                density_nodal,iconQ1,g0,R1,R2,rho_m,gravity_model,\
                                rhoblob,Rblob,zblob,rho_c)

    print("export to vtu file........................(%.3fs)" % (timing.time() - start))

    ###############################################################################

    if exp>0:
       print(spacing+" -> EARTH4D | nelr= %d v_r %.5f %.5f " %\
       (nelr,np.min(vr)/vel_unit,np.max(vr)/vel_unit))
       print(spacing+" -> EARTH4D | nelr= %d v_t %.5f %.5f " %\
       (nelr,np.min(vt)/vel_unit,np.max(vt)/vel_unit))
       print(spacing+" -> EARTH4D | nelr= %d v_rms %.5f " %\
       (nelr,vrms/vel_unit))

    ###############################################################################
    # M is in (x,z) y=0 plane, so yM=0
    # massq contains factor 2 pi too much, so we must divide by 2*pi here
    ###############################################################################

    if compute_gravity and axisymmetric:

       start=timing.time()

       print("------------------------------")
       print(" compute gravity ")
       print("------------------------------")

       print('np_grav=',np_grav)
       print('height=',height/1e3,'km')

       xM=np.zeros(np_grav,dtype=np.float64)     
       zM=np.zeros(np_grav,dtype=np.float64)     
       gvect_x=np.zeros(np_grav,dtype=np.float64)   
       gvect_y=np.zeros(np_grav,dtype=np.float64)   
       gvect_z=np.zeros(np_grav,dtype=np.float64)   
       angleM=np.zeros(np_grav,dtype=np.float64)   
       gnorm=np.zeros(np_grav,dtype=np.float64)   

       for i in range(0,np_grav):

           start2=timing.time()

           angleM[i]=np.pi/2/(np_grav-1)*i
           xM[i]=(R2+height)*np.sin(angleM[i])
           zM[i]=(R2+height)*np.cos(angleM[i])

           gvect_x[i],gvect_y[i],gvect_z[i]=\
           compute_gravity_at_point(xM[i],zM[i],nel,nqel,zq,radq,thetaq,massq,nel_phi)

           tend2=timing.time()
           print("point ",i,"/",np_grav,\
                 "-> time: %.3fs, time per elt per slice: %.9fs" %\
                   (tend2-start2,(tend2-start2)/nel/nel_phi))
       #end for i

       gnorm[:]=np.sqrt(gvect_x[:]**2+gvect_y[:]**2+gvect_z[:]**2)
       gnorm_avrg=np.sum(gnorm)/np_grav
       rM=np.sqrt(xM**2+zM**2)

       print('     -> gvect_x (m,M):',min(gvect_x),max(gvect_x))
       print('     -> gvect_y (m,M):',min(gvect_y),max(gvect_y))
       print('     -> gvect_z (m,M):',min(gvect_z),max(gvect_z))
       print('     -> gnorm   (m,M):',min(gnorm),  max(gnorm  ))

       if exp==1:
          gnorm_th=Ggrav*4*np.pi/3*(R2**3-R1**3)*rho_m/(R2+height)**2
          print('     -> grav (th)',gnorm_th)
       else:
          gnorm_th=1e50
   
       print("compute gravity (%.3fs)" % (timing.time() - start))

       gravfile1='gravity_{:04d}.ascii'.format(istep)
       np.savetxt(gravfile1,np.array([xM,zM,rM,angleM,gvect_x,gvect_y,gvect_z,\
                  gnorm,abs(gnorm-gnorm_th)/gnorm_th]).T,fmt='%.10e',\
                  header='1:xM, 2:zM, 3:rM, 4:angleM, 5:gvect_x, 6:gvect_y,\
                           7:gvect_z, 8:gvect, 9:error')

       if istep==0: 
          gnorm0=np.zeros(np_grav,dtype=np.float64)   
          gnorm0[:]=gnorm[:]
       else:
          gnorm_rate=(gnorm-gnorm0)/dt/mGal*year
          np.savetxt('gravity_rate.ascii',np.array([angleM,gnorm_rate]).T,\
                  header='1:angleM, 2:gnorm_rate')

          print("gnorm rate (m/M)",min(gnorm_rate),max(gnorm_rate),' mGal/year')

       #--------------------------------------------

       export_gravity_to_vtu(istep,np_grav,xM,zM,gvect_x,gvect_z)

    #end if compute gravity

    ###############################################################################
    # export the nel_phi slices as vtu file for educational purposes
    ###############################################################################

    #export_slices(nel_phi,NV,nel,rad,theta,iconV,density_elemental)

    ###############################################################################
    # mesh advection
    ###############################################################################
    start = timing.time()

    if nstep>1:

       if debug:
          np.savetxt('grid_before_advection_{:04d}.ascii'.format(istep),np.array([xV,yV]).T,header='# x,y')

       for i in range(0,NV): 
           if not surfaceV[i]:
              xV[i]+=u[i]*dt
              yV[i]+=v[i]*dt

       for iel in range(0,nel):
           xP[iconP[0,iel]]=xV[iconV[0,iel]]
           xP[iconP[1,iel]]=xV[iconV[1,iel]]
           xP[iconP[2,iel]]=xV[iconV[2,iel]]
           xP[iconP[3,iel]]=xV[iconV[3,iel]]
           yP[iconP[0,iel]]=yV[iconV[0,iel]]
           yP[iconP[1,iel]]=yV[iconV[1,iel]]
           yP[iconP[2,iel]]=yV[iconV[2,iel]]
           yP[iconP[3,iel]]=yV[iconV[3,iel]]

       if debug:
          np.savetxt('grid_after_advection_{:04d}.ascii'.format(istep),np.array([xV,yV]).T,header='# x,y')
          np.savetxt('grid_distortion_{:04d}.ascii'.format(istep),np.array([xV,yV,u*dt,v*dt]).T,header='# x,y')

    print("advection step............................(%.3fs)" % (timing.time() - start))

#end for istep

print("+=================================================")
print("----------------the end--------------------------")
print("+=================================================")

###############################################################################
