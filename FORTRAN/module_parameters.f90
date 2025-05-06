!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

module module_parameters

use iso_fortran_env

implicit none

integer(int32):: nelr,nelt        ! number of elements in each direction
integer(int32):: nel              ! total number of elements
integer(int32):: NfemP,NfemV,Nfem ! dofs
integer(int32):: NV,NP         ! number of velocity and pressure nodes
integer(int32):: iel,iq
integer(int32):: xi 

logical:: debug                  ! triggers lots of additional checks & prints
logical:: export_to_vtu          ! whether vtu files are created 
logical:: export_to_ascii        ! whether solution ascii files are created 

real(real64):: R1      ! inner radius of annulus or shell 
real(real64):: R2      ! outer radius of annulus or shell 
real(real64):: g0 

character(len=2):: mapping

contains

subroutine write_params
use iso_fortran_env
use module_core 
implicit none
write(*,'(a,2i10)')    ' NV,NP             =',NV,NP
write(*,'(a,i10)')     ' NfemV             =',NfemV
write(*,'(a,i10)')     ' NfemP             =',NfemP
write(*,'(a,i10)')     ' Nfem              =',Nfem
write(*,'(a,2i10)')    ' nqel              =',nqel
write(*,'(a,l10)')     ' debug             =',debug
end subroutine

end module

!==================================================================================================!
