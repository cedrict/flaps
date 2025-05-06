!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

module module_mesh

use module_core
use iso_fortran_env

implicit none

type element

  ! V node based arrays/quantities
  integer(int32):: iconV(mV)          ! connectivity array for V nodes
  real(real64):: xV(mV),yV(mV),zV(mV) ! coordinates of V nodes 
  real(real64):: u(mV),v(mV),w(mV)    ! velocity degrees of freedom 
  real(real64):: q(mV)                ! projected pressure 
  real(real64):: rV(mV),thetaV(mV)    !  
  real(real64):: bc_valV(2*mV)    ! value of the V bc on nodes
  real(real64):: nx(mV),nz(mV) 
  logical(1):: bc_fixV(2*mV)    ! whether V bc are precribed on the node
  logical(1):: hullV(mV)              ! whether V node is on hull 
  logical(1):: surfaceV(mV)           ! whether V node is on surface
  logical(1):: cmbV(mV)               ! whether V node is on cmb

  ! P node based arrays/quantities
  integer(int32):: iconP(mP)          ! connectivity array for P nodes
  real(real64):: xP(mP),yP(mP),zP(mP) ! coordinates of P nodes 
  real(real64):: p(mP)                ! pressure degrees of freedom 

  !qpoints based arrays/quantities
  real(real64):: xq(nqel),yq(nqel),zq(nqel)    ! coordinates of q. points inside elt
  real(real64):: JxWq(nqel)                    ! jacobian*weight at q. point
  real(real64):: rhoq(nqel)                    ! density at q. points
  real(real64):: etaq(nqel)                    ! viscosity at q. points
  real(real64):: gxq(nqel),gyq(nqel),gzq(nqel) ! gravity vector at q. point

  !elemental quantities
  real(real32):: volume    ! volume of the element
  real(real32):: xc,yc,zc  ! coordinates of element center
  logical(1):: bnd1_elt
  logical(1):: bnd2_elt
  logical(1):: bnd3_elt
  logical(1):: bnd4_elt
  logical(1):: bnd5_elt
  logical(1):: bnd6_elt
  logical(1):: surface           ! whether elt is on surface
  logical(1):: cmb               ! whether elt is on cmb
end type element

type(element), dimension(:), allocatable :: mesh

end module

!==================================================================================================!
