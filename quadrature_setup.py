import numpy as np

def quadrature_setup(nqperdim):


   if nqperdim==2:
      coords=[-1/np.sqrt(3),1/np.sqrt(3)]
      weights=[1,1]

   if nqperdim==3:
      coords=[-np.sqrt(3./5.),0.,np.sqrt(3./5.)]
      weights=[5./9.,8./9.,5./9.]

   if nqperdim==4:
      qc4a=np.sqrt(3./7.+2./7.*np.sqrt(6./5.))
      qc4b=np.sqrt(3./7.-2./7.*np.sqrt(6./5.))
      qw4a=(18-np.sqrt(30.))/36.
      qw4b=(18+np.sqrt(30.))/36.
      coords=[-qc4a,-qc4b,qc4b,qc4a]
      weights=[qw4a,qw4b,qw4b,qw4a]

   if nqperdim==5:
      qc5a=np.sqrt(5.+2.*np.sqrt(10./7.))/3.
      qc5b=np.sqrt(5.-2.*np.sqrt(10./7.))/3.
      qc5c=0.
      qw5a=(322.-13.*np.sqrt(70.))/900.
      qw5b=(322.+13.*np.sqrt(70.))/900.
      qw5c=128./225.
      coords=[-qc5a,-qc5b,qc5c,qc5b,qc5a]
      weights=[qw5a,qw5b,qw5c,qw5b,qw5a]

   if nqperdim==6:
      coords=[-0.932469514203152,-0.661209386466265,\
              -0.238619186083197,+0.238619186083197,\
              +0.661209386466265,+0.932469514203152]
      weights=[0.171324492379170,0.360761573048139,\
               0.467913934572691,0.467913934572691,\
               0.360761573048139,0.171324492379170]
   if nqperdim==7:
      coords=[-0.949107912342759,-0.741531185599394,\
              -0.405845151377397,0.000000000000000,\
               0.405845151377397,0.741531185599394,\
               0.949107912342759]
      weights=[0.129484966168870,0.279705391489277,\
               0.381830050505119,0.417959183673469,\
               0.381830050505119,0.279705391489277,\
               0.129484966168870]

   nqel=nqperdim**2

   qcoords_r=np.empty(nqel,dtype=np.float64)
   qcoords_s=np.empty(nqel,dtype=np.float64)
   qweights=np.empty(nqel,dtype=np.float64)

   counterq=0
   for iq in range(0,nqperdim):
       for jq in range(0,nqperdim):
           qcoords_r[counterq]=coords[iq]
           qcoords_s[counterq]=coords[jq]
           qweights[counterq]=weights[iq]*weights[jq]
           counterq+=1
       #end for
   #end for

   return nqel,qcoords_r,qcoords_s,qweights






