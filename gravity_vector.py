from numba import jit
import numpy as np

@jit(nopython=True)
def gx(x,y,g0):
    val= -x/np.sqrt(x*x+y*y)*g0
    return val

@jit(nopython=True)
def gy(x,y,g0):
    val= -y/np.sqrt(x*x+y*y)*g0
    return val

def gravity_acceleration(x,z,R1,R2,gravity_model,g0,rhom,rhoc,rhoblob,Rblob,zblob):
    r=np.sqrt(x**2+z**2)

    #-------------------
    if gravity_model==0:
       gx=-x/r*g0
       gz=-z/r*g0

    #---------------------------------
    elif gravity_model==1: # aquarium

       C=R1**3/3*(rhoc-rhom)
       E=rhom/3*(R2**3-R1**3)+R1**3*rhoc/3
       if r<R1:
          gr=rhoc*r/3
       elif r<R2:
          gr=rhom/3*r+C/r**2
       else:
          gr=E/r**2
       gx=-x/r*gr
       gz=-z/r*gr

    #---------------------------------
    elif gravity_model==2: # aquarium+blob
       C=R1**3/3*(rhoc-rhom)
       E=rhom/3*(R2**3-R1**3)+R1**3*rhoc/3
       if r<R1:
          gr=rhoc*r/3
       elif r<R2:
          gr=rhom/3*r+C/r**2
       else:
          gr=E/r**2
       coeff=4*np.pi*Ggrav
       gx=-x/r*gr*coeff
       gz=-z/r*gr*coeff

       Mblob=4*np.pi/3*rhoblob*Rblob**3
       dist=np.sqrt(x**2+(z-zblob)**2)
       if dist<Rblob:
          gr=coeff*rhoblob*dist/3
       else:
          gr=Ggrav*Mblob/dist**2

       gx+=-x/dist*gr
       gz+=-(z-zblob)/dist*gr

    #---------------------------------
    elif gravity_model==3: # prem

       exit('abs')

    else:
       exit('gravity_model unknown')
    
    return gx,gz


###############################################################################
