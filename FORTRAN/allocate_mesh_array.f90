!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine allocate_mesh_array

use iso_fortran_env
use module_parameters,only: nel,iel,nelr,nelt,xi
use module_mesh 
use module_timing

implicit none

integer(int32):: nb_int_per_cell,nb_real32_per_cell,nb_real64_per_cell,nb_logical_per_cell

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{allocate\_mesh\_arrays}
!@@ this subroutine allocates the mesh variable to the nel size 
!@@ and estimates how much memory the mesh is going to take.
!==================================================================================================!

call tic() 

!==============================================================================!

nelt=xi*nelr

nel=nelr*nelt

allocate(mesh(nel)) ! there are nel elements/cells in the mesh

!----------------------------------------------------------

nb_int_per_cell=mV+mP
nb_real32_per_cell=0
nb_real64_per_cell=0
nb_logical_per_cell=0

write(*,'(a,f7.1,a)') shift//'(max)memory needed=',&
nel*(nb_int_per_cell*4. &
    +nb_real32_per_cell*4. &
    +nb_real64_per_cell*8. &
    +nb_logical_per_cell)/1024./1024.,'Mb'

!==============================================================================!

call toc('allocate_mesh_arrays')

end subroutine

!==================================================================================================!
!==================================================================================================!
