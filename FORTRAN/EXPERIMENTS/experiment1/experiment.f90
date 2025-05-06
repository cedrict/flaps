!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! SURGE                                                                                            ! 
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!
! available modules to build with:
! use module_arrays
! use module_constants
! use module_mesh
! use module_parameters
! use module_quadrature
!==================================================================================================!

subroutine experiment_declare_main_parameters

use iso_fortran_env
use module_parameters

implicit none

!----------------------------------------------------------

R1=3400e3
R2=6400e3

nelr=16

mapping='Q2'

xi=6

g0=10

debug=.True.

export_to_vtu=.true.

!----------------------------------------------------------

end subroutine

!==================================================================================================!

subroutine assign_material_properties_to_qpoints

use iso_fortran_env
use module_parameters,only: iel,iq,nel
use module_mesh

implicit none

!----------------------------------------------------------

do iel=1,nel

   if (sqrt(mesh(iel)%xc**2 + (mesh(iel)%zc-5200e3)**2)<300e3) then
      mesh(iel)%rhoq(:)=4000
      mesh(iel)%etaq(:)=1e22
   else
      mesh(iel)%rhoq(:)=3900
      mesh(iel)%etaq(:)=1e23
   end if

end do

!----------------------------------------------------------

end subroutine

!==================================================================================================!

subroutine experiment_define_bc_V

use iso_fortran_env
use module_mesh
use module_parameters,only: iel,nel

implicit none

integer i

!----------------------------------------------------------
!no slip on the hull

do iel=1,nel
   mesh(iel)%bc_fixV(:)=.false.
   do i=1,mV
      if (mesh(iel)%hullV(i)) then
         ! x-component
         mesh(iel)%bc_fixV(ndofV*(i-1)+1)=.True.  
         mesh(iel)%bc_valV(ndofV*(i-1)+1)=0
         ! x-component
         mesh(iel)%bc_fixV(ndofV*(i-1)+2)=.True.  
         mesh(iel)%bc_valV(ndofV*(i-1)+2)=0 
      end if
   end do
end do

!----------------------------------------------------------

end subroutine

!==================================================================================================!

subroutine analytical_solution(x,y,z,p,u,v,w)
use iso_fortran_env
implicit none
real(real64),intent(in) :: x,y,z
real(real64),intent(out) :: p,u,v,w
end subroutine

!==================================================================================================!

subroutine prescribe_velocity
use iso_fortran_env
implicit none
end subroutine

!==================================================================================================!

subroutine postprocess
use iso_fortran_env
implicit none
end subroutine

