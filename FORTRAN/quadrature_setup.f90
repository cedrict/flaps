!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine quadrature_setup

use iso_fortran_env
use module_core
use module_parameters,only: nel,iel,debug,iq,g0
use module_mesh 
use module_timing
use module_arrays,only: NVq,dNVqdr,dNVqds,NPq

implicit none

real(real64), parameter:: qc2a=1.d0/sqrt(3d0)
real(real64), parameter:: qw2a=1.d0
real(real64), parameter:: qc3a=sqrt(3d0/5d0)
real(real64), parameter:: qc3b=0d0
real(real64), parameter:: qw3a=5d0/9d0
real(real64), parameter:: qw3b=8d0/9d0

integer(int32):: jq,counterq
real(real64):: rq,sq,distq
real(real64),dimension(nqperdim),parameter::qcoords=(/-qc3a,qc3b,qc3a/)
real(real64),dimension(nqperdim),parameter::qweights=(/qw3a,qw3b,qw3a/)

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{quadrature\_setup.f90}
!@@ This subroutine computes/fills all quadrature-related arrays for each element.
!@@ It further computes the real coordinates $(x_q,y_q,z_q)$ and reduced 
!@@ coordinates $(r_q,s_q,t_q)$ of these points, and assigns them their weights.
!@@ The required constants for the quadrature schemes are in 
!@@ {\filenamefont module\_quadrature.f90}.
!@@ Note that JxWq only receives Wq and will be multiplied by J in the sanity subroutine.
!==================================================================================================!

call tic() 

!----------------------------------------------------------

write(*,'(a,i5)') shift//'nqperdim=',nqperdim
write(*,'(a,i5)') shift//'nqel=',nqel

!----------------------------------------------------------
! mapping functions should be used instead of Q2 basis fcts

do iel=1,nel
   counterq=0
   do iq=1,nqperdim
   do jq=1,nqperdim
      counterq=counterq+1
      rq=qcoords(iq)
      sq=qcoords(jq)
      call NNN(rq,sq,NVq(:,counterq),mV,orderV)
      call NNN(rq,sq,NPq(:,counterq),mP,orderP)
      call dNNNdr(rq,sq,dNVqdr(:,counterq),mV,orderV)
      call dNNNds(rq,sq,dNVqds(:,counterq),mV,orderV)
      mesh(iel)%xq(counterq)=sum(mesh(iel)%xV*NVq(:,counterq))
      mesh(iel)%yq(counterq)=0
      mesh(iel)%zq(counterq)=sum(mesh(iel)%zV*NVq(:,counterq))
      mesh(iel)%JxWq(counterq)=qweights(iq)*qweights(jq)

      distq=sqrt(mesh(iel)%xq(counterq)**2+mesh(iel)%zq(counterq)**2)
      mesh(iel)%gxq(counterq)=-mesh(iel)%xq(counterq)/distq*g0
      mesh(iel)%gyq(counterq)=0
      mesh(iel)%gzq(counterq)=-mesh(iel)%zq(counterq)/distq*g0

   end do
   end do
end do

!----------------------------------------------------------
if (debug) then

open(unit=123,file='OUTPUT/qpts.dat',status='replace')
do iel=1,nel
   do iq=1,nqel
      write(123,*) mesh(iel)%xq(iq),mesh(iel)%zq(iq)
   end do
end do

write(2345,*) limit//'quadrature_setup'//limit
write(2345,*) 'nqperdim=',nqperdim
write(2345,*) 'nqel=',nqel
write(2345,*) minval(mesh(1)%xq),maxval(mesh(1)%xq)
write(2345,*) minval(mesh(1)%yq),maxval(mesh(1)%yq)
write(2345,*) minval(mesh(1)%zq),maxval(mesh(1)%zq)
end if

!==============================================================================!

call toc('____quadrature_setup')

end subroutine

!==================================================================================================!
!==================================================================================================!
