!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

module module_core

use iso_fortran_env

implicit none

integer(int32),parameter:: orderV=2         ! basis functions polynomial order
integer(int32),parameter:: orderP=1         ! basis functions polynomial order
integer(int32),parameter:: nqperdim=3       ! number of quadrature points per dimension
integer(int32),parameter:: nqel=nqperdim**2 ! number of quadrature points per element
integer(int32),parameter:: mV=9             ! number of velocity nodes per element
integer(int32),parameter:: mP=4             ! number of pressure nodes per element
integer(int32),parameter:: ndofV=2          ! 
integer(int32),parameter:: ndofP=1          ! 
integer(int32),parameter:: mVel=mV*ndofV    ! number of velocity dofs per element

real(real64),dimension(mV),parameter:: rVnode=(/-1d0,0d0,1d0,-1d0,0d0,1d0,-1d0,0d0,1d0/)
real(real64),dimension(mV),parameter:: sVnode=(/-1d0,-1d0,-1d0, 0d0,0d0,0d0, 1d0,1d0,1d0/)

real(real64),dimension(mP),parameter:: rPnode=(/-1d0,1d0,1d0,-1d0/)
real(real64),dimension(mP),parameter:: sPnode=(/-1d0,-1d0,1d0,1d0/)

end module

!==================================================================================================!
