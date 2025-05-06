!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! SURGE                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine sanity_check

use iso_fortran_env
use module_core
use module_parameters,only: iel,iq,nel,R1,R2
use module_mesh 
use module_timing
use module_constants,only: pi

implicit none

real(real64):: dNdx(mV),dNdz(mV),volume,jcob

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{sanity\_check}
!@@ This subroutine loops over all elements and quadrature points and compute the 
!@@ total volume of the domain as $V=\sum_{iel=1}^{nel} \sum_{iq=1}^{nqel} |J|_{iq} w_{iq}$.
!==================================================================================================!

call tic() 

!==============================================================================!

volume=0
do iel=1,nel
   mesh(iel)%volume=0.
   do iq=1,nqel
      call compute_dNdx_dNdz_at_qpt(iel,iq,dNdx,dNdz,jcob)
      mesh(iel)%JxWq(iq)=jcob*mesh(iel)%JxWq(iq)
      volume=volume+mesh(iel)%JxWq(iq)
      mesh(iel)%volume=mesh(iel)%volume+mesh(iel)%JxWq(iq)
   end do
end do

write(*,'(a,es9.3)') shift//'computed volume ',volume
write(*,'(a,es9.3)') shift//'analytical volume ',pi*(R2**2-R1**2)/2

!==============================================================================!

call toc('________sanity_check')

end subroutine

!==================================================================================================!
!==================================================================================================!
