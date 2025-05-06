!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine output_solution_to_ascii

use iso_fortran_env
use module_core
use module_parameters,only: iel,nel,export_to_ascii
use module_mesh 
use module_timing

implicit none

integer(int32):: k

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{output\_solution\_to\_ascii}
!@@ This subroutine export the nodal fields to the file {solution.ascii} in the OUTPUT folder.
!==================================================================================================!

call tic() 

!==============================================================================!

if (.not.export_to_ascii) return

write(*,'(a)') shift//'-> '//'OUTPUT/solution.ascii'

open(unit=123,file='OUTPUT/solutionV.ascii',action='write')
do iel=1,nel
   do k=1,mV
      write(123,'(4es7.3)') &
            mesh(iel)%xV(k),mesh(iel)%zV(k),&
            mesh(iel)%u(k),mesh(iel)%w(k)
   end do
end do
close(123)

open(unit=123,file='OUTPUT/solutionP.ascii',action='write')
do iel=1,nel
   do k=1,mP
      write(123,'(3es7.3)') &
            mesh(iel)%xP(k),mesh(iel)%zP(k),mesh(iel)%P(k)
   end do
end do
close(123)




!==============================================================================!

call toc('_output_sol_to_ascii')

end subroutine

!==================================================================================================!
!==================================================================================================!
