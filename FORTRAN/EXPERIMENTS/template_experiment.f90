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

Lx=
Ly=
Lz=

nelx=
nely=
nelz=

geometry=

debug=

FEspace='__Q1'

solve_T=

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
   do iq=1,nqel

   end do
end do

!----------------------------------------------------------

end subroutine

!==================================================================================================!

subroutine experiment_define_bc_P
use iso_fortran_env
implicit none
end subroutine

!==================================================================================================!

subroutine experiment_define_bc_T

use iso_fortran_env
use module_mesh
use module_parameters,only: iel,nel

implicit none

!----------------------------------------------------------

do iel=1,nel
   if (mesh(iel)%bnd5_elt) then
      mesh(iel)%fix_T(1:4)=
      mesh(iel)%bc_T(1:4)=
   end if
   if (mesh(iel)%bnd6_elt) then
      mesh(iel)%fix_T(5:8)=
      mesh(iel)%bc_T(5:8)=
   end if
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

