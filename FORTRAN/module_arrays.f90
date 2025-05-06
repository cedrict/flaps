!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

module module_arrays

use iso_fortran_env
use module_core

implicit none

real(real64):: NVq(mV,nqel)
real(real64):: NPq(mP,nqel)
real(real64):: dNVqdr(mV,nqel)
real(real64):: dNVqds(mV,nqel)
real(real64):: dNVqdt(mV,nqel)
real(real64),allocatable:: rhs_f(:)
real(real64),allocatable:: rhs_h(:)
real(real64),allocatable:: rnode(:),snode(:),tnode(:)

integer(4), dimension(:,:), allocatable :: Vnode_belongs_to
integer(4), dimension(:,:), allocatable :: Pnode_belongs_to

end module

!==================================================================================================!
