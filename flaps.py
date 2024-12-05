###############################################################################
#  
#  FFFFF L     AAAAA PPPPP SSSSS
#  F     L     A   A P   P S
#  FFF   L     AAAAA PPPPP SSSSS
#  F     L     A   A P         S
#  F     LLLLL A   A P     SSSSS
#  
###############################################################################

import numpy as np
import sys as sys
import scipy
import scipy.sparse as sps
from scipy.sparse.linalg import *
import time as timing
import matplotlib.pyplot as plt
from scipy import sparse
from numba import jit

###############################################################################

from basis_functions import *
from density_Earth import *
from density_4DEarthBenchmark import *
from density_AnnulusBenchmark import *
from viscosity_Earth import *
from analytical_solution import *
from gravity_vector import *
from compute_gravity_at_point import *
from compute_strain_rate1 import *
from compute_strain_rate2 import *
from compute_strain_rate3 import *
from export_to_vtu import *
from quadrature_setup import *
from compute_normals import *
from define_mapping import *
from compute_element_center_coords import *
from compute_rho_eta_fields import *
from compute_errors import *
from project_pressure_on_Q2 import *
from constants import *

###############################################################################

#planet='Earth'
#planet='Mars'
#planet='4DEarthBenchmark'
#planet='AnnulusBenchmark'
#planet='AquariumBenchmark'
planet='MarsDisc'

###############################################################################
# list of available Earth density and viscosity profiles
###############################################################################

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

eta_model=0 # constant
eta_model=1 # steinberger
eta_model=2 # samuelA
eta_model=3 # samuelB

###############################################################################
# carry out tests on inputs 
###############################################################################

print("+=================================================")
print("+=================================================")
print("   xxxxx  x      xxxxx  xxxxx  xxxxx ")
print("   x      x      x   x  x   x  x     ")
print("   xxx    x      xxxxx  xxxxx  xxxxx ")
print("   x      x      x   x  x          x ")
print("   x      xxxxx  x   x  x      xxxxx ")
print("+=================================================")
print("+=================================================")

ndim=2   # number of dimensions
mV=9     # number of nodes making up an element
mP=4     # number of nodes making up an element
ndofV=2  # number of velocity degrees of freedom per node
ndofP=1  # number of pressure degrees of freedom 

###############################################################################
# reading parameters from command line, or not
###############################################################################

if int(len(sys.argv) == 24):
   nelr              = int(sys.argv[1])
   visu              = int(sys.argv[2])
   nqperdim          = int(sys.argv[3])
   mapping           = int(sys.argv[4])
   xi                = int(sys.argv[5])
   blob_eta          = float(sys.argv[6])
   blob_rho          = float(sys.argv[7])
   blob_z            = float(sys.argv[8])
   blob_R            = float(sys.argv[9])
   blob_theta        = float(sys.argv[10])
   blob_R1           = float(sys.argv[11])
   blob_R2           = float(sys.argv[12])
   crust_rho         = float(sys.argv[13])
   crust_eta         = float(sys.argv[14])
   crust_depth       = float(sys.argv[15])
   lithosphere_rho   = float(sys.argv[16])
   lithosphere_eta   = float(sys.argv[17])
   lithosphere_depth = float(sys.argv[18])
   uppermantle_rho   = float(sys.argv[19])
   uppermantle_eta   = float(sys.argv[20])
   uppermantle_depth = float(sys.argv[21])
   lowermantle_rho   = float(sys.argv[22])
   lowermantle_eta   = float(sys.argv[23])

   if mapping==1: mapping='Q1'
   if mapping==2: mapping='Q2'
   if mapping==3: mapping='Q3'
   if mapping==4: mapping='Q4'
   if mapping==5: mapping='Q5'
   if mapping==6: mapping='Q6'

   blob_eta        = 10**blob_eta
   crust_eta       = 10**crust_eta
   lithosphere_eta = 10**lithosphere_eta
   uppermantle_eta = 10**uppermantle_eta
   lowermantle_eta = 10**lowermantle_eta

   print(sys.argv)

else:
   nelr              = 128 # Q1 cells!
   visu              = 1
   nqperdim          = 3
   mapping           = 'Q2' 
   xi                = 6
   blob_eta          = 1e23 
   blob_rho          = 3960
   blob_z            = 4900e3
   blob_R            = 400e3
   blob_theta        = 0
   blob_R1           = 0
   blob_R2           = 0
   crust_rho         = 4000 
   crust_eta         = 1e21
   crust_depth       = 50e3
   lithosphere_rho   = 4000
   lithosphere_eta   = 1e21
   lithosphere_depth = 100e3
   uppermantle_rho   = 4000
   uppermantle_eta   = 1e21
   uppermantle_depth = 700e3
   lowermantle_rho   = 4000
   lowermantle_eta   = 1e21

DJ=False

###########################################################

if mapping=='Q1': mmapping=4
if mapping=='Q2': mmapping=9
if mapping=='Q3': mmapping=16
if mapping=='Q4': mmapping=25
if mapping=='Q5': mmapping=36
if mapping=='Q6': mmapping=49

###########################################################

#NEW
# projection is n,e,q
density_projection='q'
viscosity_projection='e'

match planet:
   case "Earth":
      nstep=1
      dt=50*year
      solve_stokes=True
      axisymmetric=True
      eta_ref=1e21
      vel_unit=0.01/year
      velunit='cm/year'
      surface_free_slip=True
      bottom_free_slip=False 
      compute_gravity=False
      gravity_model=0
      rho_core=6000
      nel_phi = 200
   
   case "Mars":
      nstep=1
      dt=50*year
      solve_stokes=True
      axisymmetric=True
      eta_ref=1e21
      vel_unit=0.01/year
      velunit='cm/year'
      surface_free_slip=True
      bottom_free_slip=False 
      compute_gravity=False
      gravity_model=0
      nel_phi = 200
         
   case "AnnulusBenchmark":
      nstep=1
      dt=0
      solve_stokes=True
      axisymmetric=False
      eta_ref=1
      vel_unit=1
      velunit=' '
      surface_free_slip=False
      bottom_free_slip=False
      compute_gravity=False
      gravity_model=0
      rho_core=0
      nel_phi = 200
        
   case "AquariumBenchmark":
      nstep=1
      dt=50*year
      solve_stokes=True
      axisymmetric=True
      eta_ref=1e21
      vel_unit=0.01/year
      velunit='cm/year'
      surface_free_slip=True
      bottom_free_slip=False 
      compute_gravity=False
      gravity_model=0
      rho_core=6000
      nel_phi = 200

   case "4DEarthBenchmark":
      nstep=1
      dt=50*year
      solve_stokes=True
      axisymmetric=True
      eta_ref=1e21
      vel_unit=0.01/year
      velunit='cm/year'
      surface_free_slip=True
      bottom_free_slip=False 
      compute_gravity=False
      gravity_model=0
      rho_core=6000
      nel_phi = 100

   case "MarsDisc":

      nstep=1
      dt=50*year
      solve_stokes=True
      axisymmetric=True
      eta_ref=1e21
      vel_unit=0.01/year
      velunit='cm/year'
      surface_free_slip=True
      bottom_free_slip=False 
      compute_gravity=False
      gravity_model=0
      rho_core=8050
      nel_phi = 50

      crust_rho         = 3050 
      crust_eta         = 1e23
      crust_depth       = 60e3

      lithosphere_rho   = 3550
      lithosphere_eta   = 1e23
      lithosphere_depth = 83e3

      uppermantle_rho   = lithosphere_rho
      uppermantle_eta   = 6e20
      uppermantle_depth = 700e3 

      lowermantle_rho   = 3550 
      lowermantle_eta   = 1e21

      blob_R1           = 2296e3-100e3
      blob_R2           = blob_R1+235e3
      blob_theta        = 1700./(3396-1100+100) #np.pi/8
      blob_rho          = lowermantle_rho-70 
      blob_eta          = 1e21

   case _:
      exit('pb1 in flaps setup: unknown planet')

###########################################################

debug=False
compute_sr1=False
compute_sr3=False

###############################################################################
# read data for planets
###############################################################################

match planet:
   ############################################################################
   case 'Earth':
   ############################################################################
      start = timing.time()
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

         np.savetxt('civs12.ascii',np.array([depths_civs12,viscA_civs12,viscB_civs12]).T,fmt='%1.4e')

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

         np.savetxt('stho08.ascii',np.array([depths_stho08,visc_stho08]).T,fmt='%1.4e')

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

         np.savetxt('density_ak135f.ascii',np.array([depth_ak135f,rho_ak135f]).T,fmt='%1.4e')

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

         np.savetxt('density_stw105.ascii',np.array([depth_stw105,rho_stw105]).T,fmt='%1.4e')

         print('     -> read rho_stw105.ascii ok')

      print("read EARTH data...........................(%.3fs)" % (timing.time() - start))

   ############################################################################
   case 'Mars':
   ############################################################################
      start = timing.time()

      R1=1700e3
      R2=3390e3
      g0=3.71 #https://en.wikipedia.org/wiki/Mars   

      if eta_model==1: #steinberger
         R1=R2-1967e3
         R_disc1 = R2-49.5e3
         R_disc2 = R2-1111.5e3
         R_disc3 = R2-1160e3
         R_disc4 = R2-1951.5e3
         eta_max=1e25

      if eta_model==2: #samuelA:
         R1=1839.5976879540331e3
         R_disc1 = 3317.7417442558781e3
         R_disc2 = 2836.6008937146739e3
         R_disc3 = 2350.4998282194360e3
         R_disc4 = 1918.9611272185618e3
         eta_max=1e25

      if eta_model==3: #samuelB:
         R1=1624.2975658322634e3
         R_disc1 = 3324.3388640802909e3
         R_disc2 = 3090.3851276356227e3
         R_disc3 = 2313.0549710614014e3
         R_disc4 = 1822.5068139999998e3
         eta_max=1e25

      print("read MARS data............................(%.3fs)" % (timing.time() - start))

   ############################################################################
   case '4DEarthBenchmark':
   ############################################################################
      R1=3400e3
      R2=6400e3
      g0=10  

   ############################################################################
   case 'MarsDisc':
   ############################################################################
      R1=1835e3
      R2=3396e3
      g0=3.72 #https://en.wikipedia.org/wiki/Mars


   ############################################################################
   case 'AnnulusBenchmark':
   ############################################################################
      R1=1
      R2=2
      g0=1

   ############################################################################
   case 'AquariumBenchmark':
   ############################################################################
      R1=3400
      R2=6400
      g0=10

   ############################################################################
   case _:
      exit('unknown planet in flaps')

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

#xV,zV,rad,theta=mesh_nodes(nelr,nelt,axisymmetric,R1,R2) TODOOOOOO

if not axisymmetric:

   nelt=xi*nelr 
   nel=nelr*nelt  
   nnr=nelr+1
   nnt=nelt
   NV=nnr*nnt  # number of V nodes

   xV=np.zeros(NV,dtype=np.float64) 
   zV=np.zeros(NV,dtype=np.float64) 
   rad=np.zeros(NV,dtype=np.float64)  
   theta=np.zeros(NV,dtype=np.float64) 

   Louter=2.*np.pi*R2
   Lr=R2-R1
   sx=Louter/float(nelt)
   sz=Lr/float(nelr)

   counter=0
   for j in range(0,nnr):
       for i in range(0,nelt):
           xV[counter]=i*sx
           zV[counter]=j*sz
           counter += 1

   counter=0
   for j in range(0,nnr):
       for i in range(0,nnt):
           x_i=xV[counter]
           z_i=zV[counter]
           t=x_i/Louter*2.*np.pi    
           xV[counter]=np.cos(t)*(R1+z_i)
           zV[counter]=np.sin(t)*(R1+z_i)
           rad[counter]=R1+z_i
           theta[counter]=np.pi/2-np.arctan2(zV[counter],xV[counter])
           if theta[counter]<0.:
              theta[counter]+=2.*np.pi
           counter+=1

else:

   nelt=xi*nelr 
   nel=nelr*nelt  
   nnr=nelr+1
   nnt=nelt+1
   NV=nnr*nnt  # number of V nodes

   xV=np.zeros(NV,dtype=np.float64) 
   zV=np.zeros(NV,dtype=np.float64) 
   rad=np.zeros(NV,dtype=np.float64)  
   theta=np.zeros(NV,dtype=np.float64) 

   Louter=np.pi*R2
   Lr=R2-R1
   sx=Louter/float(nelt)
   sz=Lr/float(nelr)

   counter=0
   for j in range(0,nnr):
       for i in range(0,nnt):
           xV[counter]=i*sx
           zV[counter]=j*sz
           counter += 1

   counter=0
   for j in range(0,nnr):
       for i in range(0,nnt):
           x_i=xV[counter]
           z_i=zV[counter]
           t=np.pi/2-x_i/Louter*np.pi 
           xV[counter]=np.cos(t)*(R1+z_i)
           zV[counter]=np.sin(t)*(R1+z_i)
           rad[counter]=R1+z_i
           theta[counter]=np.pi/2-np.arctan2(zV[counter],xV[counter])
           if i==0:
              theta[counter]=0
              xV[counter]=0
           if i==nnt-1:
              theta[counter]=np.pi
              xV[counter]=0
           counter+=1

   if debug:
      np.savetxt('grid.ascii',np.array([xV,zV,theta]).T,header='# x,y',fmt='%1.4e')

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

if debug: export_Q1_mesh_to_vtu(NV,nel,xV,zV,iconQ1)

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

if axisymmetric:
   total_volume=4*np.pi/3*(R2**3-R1**3)
else:
   total_volume=np.pi*(R2**2-R1**2)

print('  -------------------')
print('  planet=',planet)
print('  nelr=',nelr,' ; nelt=',nelt,' ; nel=',nel)
print('  NfemV=',NfemV,' ; NfemP=',NfemP, '; Nfem=',Nfem)
print('  nqel=',nqel)
print('  mapping=',mapping)
print('  axisymmetric=',axisymmetric)
print('  xi=',xi)
print('  h_r=',h_r)
print('  eta_model=',eta_model)
print('  rho_model=',rho_model)
print('  compute_gravity=',compute_gravity)
print('  blob_rho=',blob_rho)
print('  blob_eta=',blob_eta)
print('  blob_R=',blob_R)
print('  blob_z=',blob_z)
print('  blob_R1=',blob_R1)
print('  blob_R2=',blob_R2)
print('  blob_theta=',blob_theta)
print('  nel_phi=',nel_phi)
print('  solve_stokes=',solve_stokes)
print('  density_projection=',density_projection)
print('  viscosity_projection=',viscosity_projection)
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
#now that I have both connectivity arrays I can easily build xP,zP (Q1 space)
###############################################################################
start = timing.time()

xP=np.empty(NP,dtype=np.float64)  # x coordinates
zP=np.empty(NP,dtype=np.float64)  # y coordinates

for iel in range(0,nel):
    xP[iconP[0,iel]]=xV[iconV[0,iel]]
    xP[iconP[1,iel]]=xV[iconV[1,iel]]
    xP[iconP[2,iel]]=xV[iconV[2,iel]]
    xP[iconP[3,iel]]=xV[iconV[3,iel]]
    zP[iconP[0,iel]]=zV[iconV[0,iel]]
    zP[iconP[1,iel]]=zV[iconV[1,iel]]
    zP[iconP[2,iel]]=zV[iconV[2,iel]]
    zP[iconP[3,iel]]=zV[iconV[3,iel]]

print("compute xP,zP.............................(%.3fs)" % (timing.time() - start))

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
outer_element=np.zeros(nel,dtype=bool) 
cmbV=np.zeros(NV,dtype=bool) 
inner_element=np.zeros(nel,dtype=bool) 

for i in range(0,NV):
    if rad[i]/R2>1-eps:
       surfaceV[i]=True
    if rad[i]/R1<1+eps:
       cmbV[i]=True

for iel in range(0,nel):
    if surfaceV[iconV[2,iel]]:
       outer_element[iel]=True 
    if cmbV[iconV[0,iel]]:
       inner_element[iel]=True 

outerQ2=np.zeros(NV,dtype=bool) 
innerQ2=np.zeros(NV,dtype=bool) 
for iel in range(0,nel):
    if outer_element[iel]:
       outerQ2[iconV[2,iel]]=True
       outerQ2[iconV[3,iel]]=True
       outerQ2[iconV[6,iel]]=True
    if inner_element[iel]:
       innerQ2[iconV[0,iel]]=True
       innerQ2[iconV[1,iel]]=True
       innerQ2[iconV[4,iel]]=True


print("flag surf and cmb nodes+elts..............(%.3fs)" % (timing.time() - start))

###############################################################################
# compute normal vectors
###############################################################################
start = timing.time()

normal_type=1

nx,nz=compute_normals(normal_type,NV,xV,zV,mV,nqel,surfaceV,hull,\
                      R1,R2,eps,iconV,axisymmetric,rad,theta,\
                      qcoords_r,qcoords_s,qweights,)

if debug:
   np.savetxt('normals.ascii',np.array([xV[surfaceV],zV[surfaceV],\
                                        nx[surfaceV],nz[surfaceV],theta[surfaceV]]).T,fmt='%1.4e')

print("compute surface normals...................(%.3fs)" % (timing.time() - start))

###############################################################################
# define boundary conditions
###############################################################################
start = timing.time()

bc_fix = np.zeros(Nfem,dtype=bool)  
bc_val = np.zeros(Nfem,dtype=np.float64) 

match planet:
   #case 'Earth':
   #case 'Mars':

   #------------------------
   case 'AquariumBenchmark':
      for i in range(0,NV):
          #vertical wall x=0
          if axisymmetric and xV[i]/R1<eps:
             bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = 0 #u=0
             if rad[i]/R2>1-eps:
                bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0 #also v=0 for 2 pts
          #surface
          if not surface_free_slip and rad[i]/R2>(1-eps):
             bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = 0 
             bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0
          if rad[i]/R1<1+eps:
             bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = 0 
             bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0

   #------------------------
   case '4DEarthBenchmark':
      for i in range(0,NV):
          #vertical wall x=0
          if axisymmetric and xV[i]/R1<eps:
             bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = 0 #u=0
             if rad[i]/R2>1-eps:
                bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0 #also v=0 for 2 pts
          #surface
          if not surface_free_slip and rad[i]/R2>(1-eps):
             bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = 0 
             bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0
          if rad[i]/R1<1+eps:
             bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = 0 
             bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0

   #------------------------
   case 'MarsDisc':
      for i in range(0,NV):
          #vertical wall x=0
          if axisymmetric and xV[i]/R1<eps:
             bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = 0 #u=0
             if rad[i]/R2>1-eps:
                bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0 #also v=0 for 2 pts
          #surface
          if not surface_free_slip and rad[i]/R2>(1-eps):
             bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = 0 
             bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0
          if rad[i]/R1<1+eps:
             bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = 0 
             bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0

   #------------------------
   case 'AnnulusBenchmark':
      for i in range(0,NV):
         if rad[i]/R1<1+eps:
            bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = velocity_x(xV[i],zV[i],R1,R2,planet)
            bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = velocity_y(xV[i],zV[i],R1,R2,planet)
         if rad[i]/R2>(1-eps) and not surface_free_slip:
            bc_fix[i*ndofV]   = True ; bc_val[i*ndofV]   = velocity_x(xV[i],zV[i],R1,R2,planet)
            bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = velocity_y(xV[i],zV[i],R1,R2,planet)
   case _:
      exit('unknown planet in flaps boundary conditions')

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

    xmapping,zmapping=define_mapping(mapping,mmapping,xV,zV,iconV,nel,\
                                     axisymmetric,rad,theta,nelt)

    print(spacing+" -> xmapping (m,M) %.2e %.2e " %(np.min(xmapping),np.max(xmapping)))
    print(spacing+" -> zmapping (m,M) %.2e %.2e " %(np.min(zmapping),np.max(zmapping)))

    if debug:
       np.savetxt('xzmapping'+mapping+'.ascii',np.array([xmapping[0,:],zmapping[0,:]]).T,fmt='%1.4e')
       export_mapping_points_to_vtu(mapping,mmapping,xmapping,zmapping)
       export_elt_quadrature_points_to_vtu(nqperdim,nqel,qcoords_r,qcoords_s,mapping,\
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
       compute_rho_eta_fields(mV,NV,nel,xV,zV,xc,zc,iconV,R1,R2,eta_model,rho_model,\
                              blob_rho,blob_eta,blob_z,blob_R,blob_R1,blob_R2,blob_theta,\
                              crust_rho,crust_eta,crust_depth,\
                              lithosphere_rho,lithosphere_eta,lithosphere_depth,\
                              uppermantle_rho,uppermantle_eta,uppermantle_depth,\
                              lowermantle_rho,lowermantle_eta,\
                              planet)

       print(spacing+" -> eta_e (m,M) %.2e %.2e " %(np.min(viscosity_elemental),np.max(viscosity_elemental)))
       print(spacing+" -> rho_e (m,M) %.2e %.2e " %(np.min(density_elemental),np.max(density_elemental)))
       print(spacing+" -> eta_n (m,M) %.2e %.2e " %(np.min(viscosity_nodal),np.max(viscosity_nodal)))
       print(spacing+" -> rho_n (m,M) %.2e %.2e " %(np.min(density_nodal),np.max(density_nodal)))

    if debug:
       np.savetxt('xycenter'+mapping+'.ascii',np.array([xc,yc,density_elemental,viscosity_elemental]).T,fmt='%1.4e')

    print("compute rho+eta fields....................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute area and volume of elements ('sanity check')
    # note that for the volume calculation we only consider the elements x>0
    ###############################################################################
    start = timing.time()

    area=np.zeros(nel,dtype=np.float64) 
    vol_elt=np.zeros(nel,dtype=np.float64) 
    mass_elt=np.zeros(nel,dtype=np.float64) 
    radq=np.zeros(nel*nqel,dtype=np.float64) 
    thetaq=np.zeros(nel*nqel,dtype=np.float64) 
    rhoq=np.zeros(nel*nqel,dtype=np.float64) 
    massq=np.zeros(nel*nqel,dtype=np.float64) 
    etaq=np.zeros(nel*nqel,dtype=np.float64) 
    coords_xq=np.zeros(nel*nqel,dtype=np.float64) 
    coords_zq=np.zeros(nel*nqel,dtype=np.float64) 

    counterq=0
    jcb=np.zeros((2,2),dtype=np.float64)
    for iel in range(0,nel):
        for kq in range(0,nqel):
            rq=qcoords_r[kq]
            sq=qcoords_s[kq]
            weightq=qweights[kq]

            if DJ:
               NNNV=NNN(rq,sq,'Q2')
               #radq[counterq]=np.dot(NNNV[:],rad[:,iel])
               #thetaq[counterq]=np.dot(NNNV[:],theta[:,iel])
               #xq=radq[counterq]*np.sin(thetaq[counterq])
               #zq=radq[counterq]*np.cos(thetaq[counterq])
               #dtheta2=(theta[icon[2,iel]]-theta[icon[0,iel]])/2
               #dr2=(rad[icon[2,iel]]-rad[icon[0,iel]])/2
	       #jcob=dtheta2*dr2
            else:
               NNNV=NNN(rq,sq,mapping)
               dNNNVdr=dNNNdr(rq,sq,mapping)
               dNNNVds=dNNNds(rq,sq,mapping)
               xq=np.dot(NNNV[:],xmapping[:,iel])
               zq=np.dot(NNNV[:],zmapping[:,iel])
               coords_xq[counterq]=xq
               coords_zq[counterq]=zq
               radq[counterq]=np.sqrt(xq**2+zq**2)
               thetaq[counterq]=np.pi/2-np.arctan2(zq,xq)
               jcb[0,0]=np.dot(dNNNVdr[:],xmapping[:,iel])
               jcb[0,1]=np.dot(dNNNVdr[:],zmapping[:,iel])
               jcb[1,0]=np.dot(dNNNVds[:],xmapping[:,iel])
               jcb[1,1]=np.dot(dNNNVds[:],zmapping[:,iel])
               jcob = np.linalg.det(jcb)

            area[iel]+=jcob*weightq
            if xq>0: vol_elt[iel]+=jcob*weightq*2*np.pi*xq


            #NEW
            ######################################
            #assign density to quadrature points
            match density_projection:
               case "e":
                  rhoq[counterq]=density_elemental[iel]
               case "q":
                  match planet:
                     case "Earth":
                        rhoq[counterq]=density_Earth(xq,zq,R1,R2,rho_m,rho_model,\
                                                     blob_rho,blob_z,blob_R)
                     case "Mars":
                        rhoq[counterq]=density_Mars(xq,zq,R1,R2,rho_m,rho_model,\
                                                    blob_rho,blob_z,blob_R)
                     case "4DEarthBenchmark":
                        rhoq[counterq]=density_4DEarthBenchmark(xq,zq,R1,R2,crust_rho,lithosphere_rho,\
                                                                uppermantle_rho,lowermantle_rho,blob_rho,blob_z,blob_R)
                     case "AnnulusBenchmark":
                        rhoq[counterq]=density_AnnulusBenchmark(xq,zq,R1,R2)

                     case "MarsDisc":
                        rhoq[counterq]=density_MarsDisc(xq,zq,R1,R2,\
                                                        blob_rho,blob_z,blob_R,blob_R1,blob_R2,blob_theta,\
                                                        crust_rho,crust_depth,\
                                                        lithosphere_rho,lithosphere_depth,\
                                                        uppermantle_rho,uppermantle_depth,\
                                                        lowermantle_rho)


                     case _:
                        exit('pb2 in flaps rhoq: unknown planet')

               case "n":

                  NNNV=NNN(rq,sq,'Q2')
                  rhoq[counterq]=np.dot(NNNV[:],density_nodal[iconV[:,iel]])

               case  _ :
                  exit('pb3 in flaps: unknown density_projection')

            #NEW
            ######################################
            #assign viscosity to quadrature points
            match viscosity_projection:
               case "e":
                  etaq[counterq]=viscosity_elemental[iel]
               case "q":
                  match planet:
                     case "Earth":
                        etaq[counterq]=viscosity_Earth(xq,zq,R1,R2,eta_m,eta_model)
                     case "Mars":
                        etaq[counterq]=viscosity_Mars(xq,zq,R1,R2,eta_m,eta_model)
                     case "4DEarthBenchmark":
                        etaq[counterq]=viscosity_4DEarthBenchmark(xq,zq,R1,R2,crust_eta,lithosphere_eta,\
                                                                  uppermantle_eta,lowermantle_eta,blob_eta,blob_z,blob_R)
                     case "AnnulusBenchmark":
                        etaq[counterq]=1
                     case _:
                        exit('pb4 in flaps etaq: unknown planet')

               case "n":
                  NNNV=NNN(rq,sq,'Q2')
                  etaq[counterq]=np.dot(NNNV[:],viscosity_nodal[iconV[:,iel]])

               case  _  :
                  exit('pb5: unknown viscosity_projection')

            massq[counterq]=rhoq[counterq]*weightq*jcob*2*np.pi*xq
            if xq>0: mass_elt[iel]+=massq[counterq]

            counterq+=1
       #end for
    #end for


    if axisymmetric:
       print(spacing+" -> total volume (meas) %.12e | nel= %d" %(vol_elt.sum(),nel))
       print(spacing+" -> total volume (anal) %.12e" %(total_volume))
    else:
       print(spacing+" -> area (m,M) %.6e %.6e " %(np.min(area),np.max(area)))
       print(spacing+" -> total area (meas) %.12e | nel= %d" %(area.sum(),nel))
       print(spacing+" -> total area (anal) %e " %(total_volume))

    print(spacing+" -> rhoq (m,M) %.3e %.3e " %(np.min(rhoq),np.max(rhoq)))
    print(spacing+" -> etaq (m,M) %.3e %.3e " %(np.min(etaq),np.max(etaq)))
    print(spacing+" -> massq (m,M) %.3e %.3e " %(np.min(massq),np.max(massq)))

    print(spacing+" -> total mass (meas) %.12e | nel= %d" %(np.sum(mass_elt),nelr))
    
    print("sanity check.............................. %.3fs  | %d" % (timing.time()-start,Nfem))

    ###############################################################################

    if debug:
       np.savetxt('qpoints.ascii',np.array([coords_xq,coords_zq,rhoq,etaq]).T,header='# x,y,rho,eta',fmt='%1.4e')

    export_quadrature_points_to_vtu(nel*nqel,coords_xq,coords_zq,rhoq,etaq)

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
            zq=np.dot(NNNV[:],zmapping[:,iel])

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

            g_x,g_y=gravity_acceleration(xq,zq,R1,R2,gravity_model,g0,lowermantle_rho,\
                                         rho_core,blob_rho,blob_R,blob_z)
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

        if surface_free_slip and outer_element[iel]==1:
           for k in range(0,mV):
               inode=iconV[k,iel]
               #if surfaceV[inode] and xV[inode]>0 and (not bc_fix[inode*ndofV]):
               if surfaceV[inode] and (not bc_fix[inode*ndofV]):

                  #print(xV[inode],zV[inode],np.arctan2(nx[inode],nz[inode]),theta[inode])
                  #alpha=-np.arctan2(nx[inode],nz[inode])

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

    print("build FE system........................... %.3fs  | %d" % (timing.time()-start,Nfem))

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

       #sol = scipy.sparse.linalg.gmres(sparse_matrix, rhs, restart=2000,tol=1e-8)[0]
       #sol = scipy.sparse.linalg.lgmres(sparse_matrix, rhs,atol=1e-16,tol=1e-5)[0]

    else:
       sol=np.zeros(Nfem,dtype=np.float64)  

    print("solving system............................ %.3fs  | %d" % (timing.time()-start,Nfem))

    ###############################################################################
    # put solution into separate x,y velocity arrays
    ###############################################################################
    start = timing.time()

    u,v=np.reshape(sol[0:NfemV],(NV,2)).T
    p=sol[NfemV:NfemV+NfemP]*eta_ref/h_r

    if debug:
       np.savetxt('velocity.ascii',np.array([xV,zV,u/vel_unit,v/vel_unit]).T,header='# x,y,u,v',fmt='%1.4e')
       np.savetxt('pressure.ascii',np.array([xP,zP,p]).T,header='# x,y,p',fmt='%1.4e')
 
    vr= np.sin(theta)*u+np.cos(theta)*v
    vt= np.cos(theta)*u-np.sin(theta)*v
   
    vel=np.zeros(NV,dtype=np.float64)  
    vel[:]=np.sqrt(u[:]**2+v[:]**2)

    print(spacing+" -> nelr= %d | u   (m,M) %.7e %.7e " %(nelr,np.min(u)/vel_unit,np.max(u)/vel_unit),velunit)
    print(spacing+" -> nelr= %d | v   (m,M) %.7e %.7e " %(nelr,np.min(v)/vel_unit,np.max(v)/vel_unit),velunit)
    print(spacing+" -> nelr= %d | v_r (m,M) %.7e %.7e " %(nelr,np.min(vr)/vel_unit,np.max(vr)/vel_unit),velunit)
    print(spacing+" -> nelr= %d | v_t (m,M) %.7e %.7e " %(nelr,np.min(vt)/vel_unit,np.max(vt)/vel_unit),velunit)
    print(spacing+" -> nelr= %d | vel (m,M) %.7e %.7e  "%(nelr,np.min(vel)/vel_unit,np.max(vel)/vel_unit),velunit)

    print("reshape solution..........................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute strain rate - center to nodes - method 1
    ###############################################################################
    start = timing.time()

    sr1=np.zeros(NV,dtype=np.float64)  

    if compute_sr1:

       exx1,ezz1,exz1,exxc,ezzc,exzc=compute_strain_rate1(nel,mV,NV,iconV,mapping,xmapping,zmapping,u,v)

       print(spacing+" -> exxc (m,M) %e %e " %(np.min(exxc),np.max(exxc)))
       print(spacing+" -> ezzc (m,M) %e %e " %(np.min(ezzc),np.max(ezzc)))
       print(spacing+" -> exzc (m,M) %e %e " %(np.min(exzc),np.max(exzc)))
       print(spacing+" -> exx1 (m,M) %e %e " %(np.min(exx1),np.max(exx1)))
       print(spacing+" -> ezz1 (m,M) %e %e " %(np.min(ezz1),np.max(ezz1)))
       print(spacing+" -> exz1 (m,M) %e %e " %(np.min(exz1),np.max(exz1)))

       sr1=np.zeros(NV,dtype=np.float64)  
       sr1[:]=np.sqrt(0.5*(exx1[:]**2+ezz1[:]**2)+exz1[:]**2)
       np.savetxt('sr1_R1.ascii',np.array([xV[innerQ2],zV[innerQ2],sr1[innerQ2],theta[innerQ2]]).T,fmt='%1.4e')
       np.savetxt('sr1_R2.ascii',np.array([xV[outerQ2],zV[outerQ2],sr1[outerQ2],theta[outerQ2]]).T,fmt='%1.4e')

    print("compute strain rate meth-1................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute strain rate - corners to nodes - method 2
    ###############################################################################
    start = timing.time()

    if solve_stokes:
       exx2,ezz2,exz2=compute_strain_rate2(nel,mV,NV,iconV,mapping,xmapping,zmapping,u,v)
    else:
       exx2=np.zeros(NV,dtype=np.float64)  
       ezz2=np.zeros(NV,dtype=np.float64)  
       exz2=np.zeros(NV,dtype=np.float64)  

    print(spacing+" -> exx2 (m,M) %e %e " %(np.min(exx2),np.max(exx2)))
    print(spacing+" -> ezz2 (m,M) %e %e " %(np.min(ezz2),np.max(ezz2)))
    print(spacing+" -> exz2 (m,M) %e %e " %(np.min(exz2),np.max(exz2)))

    sr2=np.zeros(NV,dtype=np.float64)  
    sr2[:]=np.sqrt(0.5*(exx2[:]**2+ezz2[:]**2)+exz2[:]**2)

    if debug:
       np.savetxt('strainrate'+str(istep)+'.ascii',np.array([xV,zV,exx2,ezz2,exz2]).T,fmt='%1.4e')

    print("compute strain rate meth-2................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute strain rate - via mass matrix - method 3
    ###############################################################################
    start = timing.time()

    sr3=np.zeros(NV,dtype=np.float64)  

    if compute_sr3:

       exx3,ezz3,exz3=compute_strain_rate3(nel,mV,NV,nqel,iconV,mapping,xmapping,zmapping,u,v,\
                                           qcoords_r,qcoords_s,qweights)

       print(spacing+" -> exx3 (m,M) %e %e " %(np.min(exx3),np.max(exx3)))
       print(spacing+" -> ezz3 (m,M) %e %e " %(np.min(ezz3),np.max(ezz3)))
       print(spacing+" -> exz3 (m,M) %e %e " %(np.min(exz3),np.max(exz3)))

       sr3[:]=np.sqrt(0.5*(exx3[:]**2+ezz3[:]**2)+exz3[:]**2)
       np.savetxt('sr3_R1.ascii',np.array([xV[0:nnt],zV[0:nnt],sr3[0:nnt],theta[0:nnt]]).T,fmt='%1.4e')
       np.savetxt('sr3_R2.ascii',np.array([xV[outerQ2],zV[outerQ2],sr3[outerQ2],theta[outerQ2]]).T,fmt='%1.4e')

    print("compute strain rate meth-3................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # project pressure onto Q2 mesh
    ###############################################################################
    start = timing.time()

    if solve_stokes:
       q=project_pressure_on_Q2(NV,nel,mV,mP,p,iconP,iconV)
    else:
       q=np.zeros(NV,dtype=np.float64)  

    print(spacing+" -> q (m,M) %e %e " %(np.min(q),np.max(q)))

    if debug:
       np.savetxt('pressure'+str(istep)+'.ascii',np.array([xV,zV,q]).T,fmt='%1.4e')

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

    if solve_stokes:
       for i in range(0,NV):
           if xV[i]>=0:
              e_rr2[i]=exx2[i]*np.sin(theta[i])**2+\
                      2*exz2[i]*np.sin(theta[i])*np.cos(theta[i])+\
                      ezz2[i]*np.cos(theta[i])**2
              e_tt2[i]=exx2[i]*np.cos(theta[i])**2-\
                      2*exz2[i]*np.sin(theta[i])*np.cos(theta[i])+\
                      ezz2[i]*np.sin(theta[i])**2
              e_rt2[i]=(exx2[i]-ezz2[i])*np.sin(theta[i])*np.cos(theta[i])+\
                      exz2[i]*(-np.sin(theta[i])**2+\
                      np.cos(theta[i])**2)

       print(spacing+" -> e_rr (m,M) %e %e | nel= %d" %(np.min(e_rr2),np.max(e_rr2),nel))
       print(spacing+" -> e_tt (m,M) %e %e | nel= %d" %(np.min(e_tt2),np.max(e_tt2),nel))
       print(spacing+" -> e_rt (m,M) %e %e | nel= %d" %(np.min(e_rt2),np.max(e_rt2),nel))

    print("compute strain rate in sph. coords........(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # normalise pressure: zero average pressure on surface
    # note that the integration could be improved (I currently only sample the 
    # pressure in 1 point in the middle of the edge)
    ###############################################################################
    start = timing.time()

    if not axisymmetric:
       #if planet=='AnnulusBenchmark':
       #   poffset=np.sum(q[0:2*nelt])/(2*nelt)
       #else:
       poffset=0
       for iel in range(0,nel):
           if outer_element[iel]:
              dtheta=2*np.pi/nelt 
              pmean=0.5*(p[iconP[2,iel]]+p[iconP[3,iel]])
              poffset+=dtheta*pmean
       poffset/=2*np.pi

       q-=poffset
       p-=poffset

    else: 

       poffset=0
       for iel in range(0,nel):
           if outer_element[iel]:
              dtheta=theta[iconV[2,iel]]-theta[iconV[3,iel]]
              pmean=0.5*(p[iconP[2,iel]]+p[iconP[3,iel]])
              poffset+=np.sin((theta[iconV[2,iel]]+theta[iconV[3,iel]])/2)*dtheta\
                       *2*np.pi*R2**2 * pmean
       poffset/=4*np.pi*R2**2
       q-=poffset
       p-=poffset

    print(spacing+" -> p (m,M) %e %e " %(np.min(p),np.max(p)))
    print(spacing+" -> q (m,M) %e %e " %(np.min(q),np.max(q)))

    if debug: 
       np.savetxt('pressure.ascii',np.array([xP,zP,p]).T,header='# x,z,p',fmt='%1.4e')

    print("normalise pressure........................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute average pressure at bottom (needed for dyn topo)
    ###############################################################################

    if not axisymmetric:

       exit('not done')

    else: 

       bottom_p_avrg=0
       for iel in range(0,nel):
           if inner_element[iel]:
              dtheta=theta[iconV[1,iel]]-theta[iconV[0,iel]]
              pmean=0.5*(p[iconP[0,iel]]+p[iconP[1,iel]])
              bottom_p_avrg+=np.sin((theta[iconV[0,iel]]+theta[iconV[1,iel]])/2)*dtheta\
                            *2*np.pi*R1**2*pmean
       bottom_p_avrg/=4*np.pi*R1**2

    print(spacing+" -> bottom <p> (m,M) %e " %(bottom_p_avrg))

    ###############################################################################

    p_th=np.zeros(NP,dtype=np.float64)
    q_th=np.zeros(NV,dtype=np.float64)
    for i in range(0,NP):
        p_th[i]=pressure(xP[i],zP[i],R1,R2,lowermantle_rho,g0,planet)
    for i in range(0,NV):
        q_th[i]=pressure(xV[i],zV[i],R1,R2,lowermantle_rho,g0,planet)
    
    if debug:    
       np.savetxt('p_th.ascii',np.array([xP,zP,p_th]).T,header='# x,z,p',fmt='%1.4e')
       np.savetxt('q_th.ascii',np.array([xV,zV,q_th]).T,header='# x,z,p',fmt='%1.4e')
       np.savetxt('q_th_R2.ascii',np.array([theta[outerQ2],q_th[outerQ2]]).T,fmt='%1.4e')

    ###############################################################################
    # since I have zeroed the the spherical components of the strain rate for x<0
    # then I also zero the dynamic topography
    ###############################################################################
    start = timing.time()

    dyn_topo_nodal=np.zeros(NV,dtype=np.float64)
    for i in range(0,NV):
        if surfaceV[i] and xV[i]>=0:
           dyn_topo_nodal[i]= -(2*viscosity_nodal[i]*e_rr2[i]-q[i])/(crust_rho*g0) 
        if cmbV[i] and xV[i]>=0:
           dyn_topo_nodal[i]= -(2*viscosity_nodal[i]*e_rr2[i]-(q[i]-bottom_p_avrg))/((lowermantle_rho-rho_core)*g0) 

    print("compute dynamic topography................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # export fields at both surfaces
    ###############################################################################
    start = timing.time()

    np.savetxt('qqq_R1_'+str(istep)+'.ascii',np.array([theta[innerQ2],q[innerQ2]]).T,fmt='%1.6e')
    np.savetxt('qqq_R2_'+str(istep)+'.ascii',np.array([theta[outerQ2],q[outerQ2]]).T,fmt='%1.6e')
    np.savetxt('sr2_R2_'+str(istep)+'.ascii',np.array([theta[outerQ2],sr2[outerQ2]]).T,fmt='%1.4e')
    np.savetxt('vel_R2_'+str(istep)+'.ascii',np.array([theta[outerQ2],vel[outerQ2]/vel_unit,vr[outerQ2]/vel_unit,vt[outerQ2]/vel_unit]).T,fmt='%1.4e')
    np.savetxt('err_R2_'+str(istep)+'.ascii',np.array([theta[outerQ2],e_rr2[outerQ2]]).T,fmt='%1.4e')
    np.savetxt('d_t_R1_'+str(istep)+'.ascii',np.array([theta[innerQ2],dyn_topo_nodal[innerQ2]]).T,fmt='%1.4e')
    np.savetxt('d_t_R2_'+str(istep)+'.ascii',np.array([theta[outerQ2],dyn_topo_nodal[outerQ2]]).T,fmt='%1.4e')
    #if axisymmetric:
    #   leftV=np.zeros(NV,dtype=bool) 
    #   for i in range(0,NV):
    #       if xV[i]<eps and zV[i]>0: leftV[i]=True
    #   np.savetxt('sr2_left_'+str(istep)+'.ascii',np.array([xV[leftV],zV[leftV],sr2[leftV],rad[leftV]]).T,fmt='%1.4e')
    #   np.savetxt('vel_left_'+str(istep)+'.ascii',np.array([xV[leftV],zV[leftV],vel[leftV],rad[leftV]]).T,fmt='%1.4e')

    print("export p&q on R1,R2.......................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute errors and vrms
    ###############################################################################
    start = timing.time()

    if solve_stokes:
       errv,errp,vrms=compute_errors(nel,nqel,mapping,xmapping,zmapping,qcoords_r,\
                                     qcoords_s,qweights,R1,R2,lowermantle_rho,g0,u,v,p,\
                                     axisymmetric,iconV,iconP,total_volume,planet)

       print(spacing+' -> nelr= %d ; vrms= %.12e' %(nelr,vrms/vel_unit))
       print(spacing+" -> nelr= %d ; errv= %.12e ; errp= %.14e " %(nelr,errv,errp))

    else:
       errv=errp=vrms=0

    print("compute errors............................(%.3fs)" % (timing.time() - start))

    ###############################################################################
    # compute gravity acceleration on pressure nodes
    ###############################################################################
    start = timing.time()

    self_gravitation=False

    if self_gravitation:

       g_x=np.zeros(NP,dtype=np.float64)   
       g_y=np.zeros(NP,dtype=np.float64)   
       g_z=np.zeros(NP,dtype=np.float64)   
       gnorm=np.zeros(NP,dtype=np.float64)   

       for i in range(0,NP):
           if i%100==0: print('i=',i,'/',NP)
           g_x[i],g_y[i],g_z[i]=\
           compute_gravity_at_point(xP[i],zP[i],nel,nqel,zq,radq,thetaq,massq,nel_phi)
           #print(xV[i],zV[i],g_x[i],g_y[i],g_z[i])

       gnorm[:]=np.sqrt(g_x[:]**2+g_y[:]**2+g_z[:]**2)

       np.savetxt('gravity_self.ascii',np.array([xP,zP,g_x,g_y,g_z,gnorm]).T,fmt='%1.4e')

       print(spacing+"     -> g_x (m,M) %.4e %.4e " %(np.min(g_x),np.max(g_x)))
       print(spacing+"     -> g_y (m,M) %.4e %.4e " %(np.min(g_y),np.max(g_y)))
       print(spacing+"     -> g_z (m,M) %.4e %.4e " %(np.min(g_z),np.max(g_z)))

       print("self gravity (%.3fs)" % (timing.time() - start))

    ###############################################################################
    # plot of solution
    ###############################################################################
    start = timing.time()

    if visu==1:
       export_solution_to_vtu(istep,NV,nel,xV,zV,iconV,u,v,vr,vt,q,vel_unit,rad,\
                              theta,nx,nz,sr1,sr2,sr3,density_nodal,density_elemental,\
                              viscosity_nodal,viscosity_elemental,R1,R2,lowermantle_rho,\
                              gravity_model,g0,rho_core,blob_rho,blob_R,blob_z,hull,\
                              inner_element,outer_element,innerQ2,outerQ2,bc_fix,\
                              e_rr2,e_tt2,e_rt2,vol_elt,mass_elt,planet)

    print("export to vtu file........................(%.3fs)" % (timing.time() - start))

    ###############################################################################

    if planet=='4DEarthBenchmark':
       print(spacing+" -> DTGbenchmark | nelr= %d v_r %.5f %.5f " %\
       (nelr,np.min(vr)/vel_unit,np.max(vr)/vel_unit))
       print(spacing+" -> DTGbenchmark | nelr= %d v_t %.5f %.5f " %\
       (nelr,np.min(vt)/vel_unit,np.max(vt)/vel_unit))

       if density_projection=='e':
          print(spacing+" -> DTGbenchmark | mass_blob",nelr,np.sum(density_elemental[np.where(density_elemental<4000)]*vol_elt[np.where(density_elemental<4000)]))

    ###############################################################################
    # M is in (x,z) y=0 plane, so yM=0
    # massq contains factor 2 pi too much, so we must divide by 2*pi here
    ###############################################################################

    if compute_gravity and axisymmetric:

       start=timing.time()

       print("------------------------------")
       print(" compute gravity ")
       print("------------------------------")

       np_grav=100 # nb of satellites positions
       height=1e3
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

           angleM[i]=np.pi/(np_grav-1)*i
           xM[i]=(R2+height)*np.sin(angleM[i])
           zM[i]=(R2+height)*np.cos(angleM[i])

           gvect_x[i],gvect_y[i],gvect_z[i]=\
           compute_gravity_at_point(xM[i],zM[i],nel,nqel,coords_zq,radq,thetaq,massq,nel_phi)

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

       gnorm_th=1 #Ggrav*4*np.pi/3*(R2**3-R1**3)*rho_m/(R2+height)**2
       print('     -> grav (th)',gnorm_th)
   
       print("compute gravity (%.3fs)" % (timing.time() - start))

       gravfile1='gravity_{:04d}.ascii'.format(istep)
       np.savetxt(gravfile1,np.array([xM,zM,rM,angleM,gvect_x,gvect_y,gvect_z,\
                  gnorm,abs(gnorm-gnorm_th)/gnorm_th]).T,\
                  header='1:xM, 2:zM, 3:rM, 4:angleM, 5:gvect_x, 6:gvect_y,\
                           7:gvect_z, 8:gvect, 9:error',fmt='%1.4e')

       if istep==0: 
          gnorm0=np.zeros(np_grav,dtype=np.float64)   
          gnorm0[:]=gnorm[:]
       else:
          gnorm_rate=(gnorm-gnorm0)/dt/mGal*year
          np.savetxt('gravity_rate.ascii',np.array([angleM,gnorm_rate]).T,\
                  header='1:angleM, 2:gnorm_rate',fmt='%1.4e')

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
          np.savetxt('grid_before_advection_{:04d}.ascii'.format(istep),np.array([xV,zV]).T,header='# x,y',fmt='%1.4e')

       for i in range(0,NV): 
           if not surfaceV[i]:
              xV[i]+=u[i]*dt
              zV[i]+=v[i]*dt

       for iel in range(0,nel):
           xP[iconP[0,iel]]=xV[iconV[0,iel]]
           xP[iconP[1,iel]]=xV[iconV[1,iel]]
           xP[iconP[2,iel]]=xV[iconV[2,iel]]
           xP[iconP[3,iel]]=xV[iconV[3,iel]]
           zP[iconP[0,iel]]=zV[iconV[0,iel]]
           zP[iconP[1,iel]]=zV[iconV[1,iel]]
           zP[iconP[2,iel]]=zV[iconV[2,iel]]
           zP[iconP[3,iel]]=zV[iconV[3,iel]]

       if debug:
          np.savetxt('grid_after_advection_{:04d}.ascii'.format(istep),np.array([xV,zV]).T,header='# x,y',fmt='%1.4e')
          np.savetxt('grid_distortion_{:04d}.ascii'.format(istep),np.array([xV,zV,u*dt,v*dt]).T,header='# x,y',fmt='%1.4e')

    print("advection step............................(%.3fs)" % (timing.time() - start))

#end for istep

print("+=================================================")
print("----------------the end--------------------------")
print("+=================================================")

###############################################################################
