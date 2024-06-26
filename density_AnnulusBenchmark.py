from numba import jit
import numpy as np
import math as math # <-- remove ?

#@jit(nopython=True)
def density_AnnulusBenchmark(x,y,R1,R2):
    k=4
    r=np.sqrt(x*x+y*y)
    theta=math.atan2(y,x)
    A=2.*(math.log(R1)-math.log(R2))/(R2**2*math.log(R1)-R1**2*math.log(R2) )
    B=(R2**2-R1**2)/(R2**2*math.log(R1)-R1**2*math.log(R2) )
    fr=A*r+B/r
    fpr=A-B/r**2
    gr=A/2.*r + B/r*math.log(r) - 1./r
    gpr=A/2.+B/r**2*(1.-math.log(r))+1./r**2
    gppr=-B/r**3*(3.-2.*math.log(r))-2./r**3
    alephr=gppr - gpr/r -gr/r**2*(k**2-1.) +fr/r**2  +fpr/r
    val=k*math.sin(k*theta)*alephr #+ rho0 
    return val
