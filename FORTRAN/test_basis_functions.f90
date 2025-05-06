!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine test_basis_functions

use iso_fortran_env
use module_core
use module_parameters,only: nel,debug,iel,iq
use module_arrays,only: NVq 
use module_mesh 
use module_timing

implicit none

real(real64):: dNdx(mV),dNdz(mV)
real(real64):: dudxq,dudzq,jcob,uq

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{test\_basis\_functions}
!@@ This subroutine tests the consistency of the basis functions. 
!@@ An analytical velocity field is prescribed (constant, linear or quadratic) and the 
!@@ corresponding values are computed onto the quadrature points via the 
!@@ (derivatives of the) basis functions.
!@@ It generates three ascii files in the {\foldernamefont OUTPUT} folder.
!==================================================================================================!

call tic() 

!==============================================================================!

if (debug) then

open(unit=123,file='OUTPUT/TEST/test_basis_functions_constant.ascii',action='write')
open(unit=124,file='OUTPUT/TEST/test_basis_functions_linear.ascii',action='write')
open(unit=125,file='OUTPUT/TEST/test_basis_functions_quadratic.ascii',action='write')

do iel=1,nel
   !-------------------------
   mesh(iel)%u=1
   do iq=1,nqel
      uq=sum(NVq(1:mV,iq)*mesh(iel)%u(1:mV)) 
      call compute_dNdx_dNdz_at_qpt(iel,iq,dNdx,dNdz,jcob)
      dudxq=sum(dNdx(1:mV)*mesh(iel)%u(1:mV)) 
      dudzq=sum(dNdz(1:mV)*mesh(iel)%u(1:mV)) 
      write(123,*) mesh(iel)%xq(iq),mesh(iel)%zq(iq),uq,dudxq,dudzq,jcob
   end do
      write(123,*) '#--------------'
   !-------------------------
   mesh(iel)%u=mesh(iel)%xV
   do iq=1,nqel
      uq=sum(NVq(1:mV,iq)*mesh(iel)%u(1:mV)) 
      call compute_dNdx_dNdz_at_qpt(iel,iq,dNdx,dNdz,jcob)
      dudxq=sum(dNdx(1:mV)*mesh(iel)%u(1:mV)) 
      dudzq=sum(dNdz(1:mV)*mesh(iel)%u(1:mV)) 
      write(124,*) mesh(iel)%xq(iq),mesh(iel)%zq(iq),uq,dudxq,dudzq,jcob
   end do
   !-------------------------
   mesh(iel)%u=mesh(iel)%xV**2
   do iq=1,nqel
      uq=sum(NVq(1:mV,iq)*mesh(iel)%u(1:mV)) 
      call compute_dNdx_dNdz_at_qpt(iel,iq,dNdx,dNdz,jcob)
      dudxq=sum(dNdx(1:mV)*mesh(iel)%u(1:mV)) 
      dudzq=sum(dNdz(1:mV)*mesh(iel)%u(1:mV)) 
      write(125,*) mesh(iel)%xq(iq),mesh(iel)%zq(iq),uq,dudxq,dudzq,jcob
   end do
end do
close(123) ; write(*,'(a)') shift//'-> OUTPUT/TEST/test_basis_functions_constant.ascii'
close(124) ; write(*,'(a)') shift//'-> OUTPUT/TEST/test_basis_functions_linear.ascii'
close(125) ; write(*,'(a)') shift//'-> OUTPUT/TEST/test_basis_functions_quadratic.ascii'

! reset fields to zero
do iel=1,nel
   mesh(iel)%u(:)=0
   mesh(iel)%v(:)=0
   mesh(iel)%w(:)=0
end do

end if ! debug

!==============================================================================!

call toc('test_basis_functions')

end subroutine

!==================================================================================================!
!==================================================================================================!
