###############################################################################

from numba import jit
import numpy as np
import math as math # <-- remove ?

###############################################################################

@jit(nopython=True)
def velocity_x(x,y,R1,R2,planet):
    match planet:
       case 'AnnulusBenchmark':
          k=4
          r=np.sqrt(x*x+y*y)
          theta=math.atan2(y,x)
          A=2.*(math.log(R1)-math.log(R2))/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          B=(R2**2-R1**2)/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          fr=A*r+B/r
          fpr=A-B/r**2
          gr=A/2.*r + B/r*math.log(r) - 1./r
          hr=(2*gr-fr)/r
          vr=k *gr * math.sin (k * theta)
          vtheta = fr *math.cos(k* theta)
          val=vr*math.cos(theta)-vtheta*math.sin(theta)
       case 'AquariumBenchmark':
          val=0
       case _:
          val=0
    return val

###############################################################################

@jit(nopython=True)
def velocity_y(x,y,R1,R2,planet):
    match planet:
       case 'AnnulusBenchmark':
          k=4
          r=np.sqrt(x*x+y*y)
          theta=math.atan2(y,x)
          A=2.*(math.log(R1)-math.log(R2))/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          B=(R2**2-R1**2)/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          fr=A*r+B/r
          fpr=A-B/r**2
          gr=A/2.*r + B/r*math.log(r) - 1./r
          hr=(2*gr-fr)/r
          vr=k *gr * math.sin (k * theta)
          vtheta = fr *math.cos(k* theta)
          val=vr*math.sin(theta)+vtheta*math.cos(theta)
       case 'AquariumBenchmark':
          val=0
       case _:
          val=0
    return val

###############################################################################

@jit(nopython=True)
def pressure(x,y,R1,R2,rho0,g0,planet):
    r=np.sqrt(x*x+y*y)
    match planet:
       case 'AnnulusBenchmark':
          k=4
          theta=math.atan2(y,x)
          A=2.*(math.log(R1)-math.log(R2))/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          B=(R2**2-R1**2)/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          fr=A*r+B/r
          gr=A/2.*r + B/r*math.log(r) - 1./r
          hr=(2*gr-fr)/r
          val=k*hr*math.sin(k*theta) #+ rho0*g0*(r-R2)
       case 'AquariumBenchmark':
          val=rho0*g0*(R2-r)
       case _:
          val=0
    return val

###############################################################################

@jit(nopython=True)
def sr_xx(x,y,R1,R2,planet):
    match planet:
       case 'AnnulusBenchmark':
          k=4
          r=np.sqrt(x*x+y*y)
          theta=math.atan2(y,x)
          A=2.*(math.log(R1)-math.log(R2))/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          B=(R2**2-R1**2)/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          gr=A/2.*r + B/r*math.log(r) - 1./r
          gpr=A/2 + B*((1-math.log(r)) / r**2 ) +1./r**2
          fr=A*r+B/r
          fpr=A-B/r**2
          err=gpr*k*math.sin(k*theta)
          ert=0.5*(k**2/r*gr+fpr-fr/r)*math.cos(k*theta)
          ett=(gr-fr)/r*k*math.sin(k*theta)
          val=err*(math.cos(theta))**2\
             +ett*(math.sin(theta))**2\
             -2*ert*math.sin(theta)*math.cos(theta)
       case 'AquariumBenchmark':
          val=0
       case _:
          val=0
    return val

###############################################################################

@jit(nopython=True)
def sr_yy(x,y,R1,R2,planet):
    match planet:
       case 'AnnulusBenchmark':
          k=4
          r=np.sqrt(x*x+y*y)
          theta=math.atan2(y,x)
          A=2.*(math.log(R1)-math.log(R2))/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          B=(R2**2-R1**2)/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          gr=A/2.*r + B/r*math.log(r) - 1./r
          gpr=A/2 + B*((1-math.log(r)) / r**2 ) +1./r**2
          fr=A*r+B/r
          fpr=A-B/r**2
          err=gpr*k*math.sin(k*theta)
          ert=0.5*(k**2/r*gr+fpr-fr/r)*math.cos(k*theta)
          ett=(gr-fr)/r*k*math.sin(k*theta)
          val=err*(math.sin(theta))**2\
             +ett*(math.cos(theta))**2\
             +2*ert*math.sin(theta)*math.cos(theta)
       case 'AquariumBenchmark':
          val=0
       case _:
          val=0
    return val

###############################################################################

@jit(nopython=True)
def sr_xy(x,y,R1,R2,planet):
    match planet:
       case 'AnnulusBenchmark':
          k=4
          r=np.sqrt(x*x+y*y)
          theta=math.atan2(y,x)
          A=2.*(math.log(R1)-math.log(R2))/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          B=(R2**2-R1**2)/(R2**2*math.log(R1)-R1**2*math.log(R2) )
          gr=A/2.*r + B/r*math.log(r) - 1./r
          gpr=A/2 + B*((1-math.log(r)) / r**2 ) +1./r**2
          fr=A*r+B/r
          fpr=A-B/r**2
          err=gpr*k*math.sin(k*theta)
          ert=0.5*(k**2/r*gr+fpr-fr/r)*math.cos(k*theta)
          ett=(gr-fr)/r*k*math.sin(k*theta)
          val=ert*(math.cos(theta)**2-math.sin(theta)**2)\
             +(err-ett)*math.cos(theta)*math.sin(theta)
       case 'AquariumBenchmark':
          val=0
       case _:
          val=0
    return val

###############################################################################
