!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! ELEFANT                                                                        C. Thieulot       !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine compute_belongs

use iso_fortran_env
use module_parameters, only: nel,iel,debug,NV,NP
use module_mesh 
use module_arrays, only: Vnode_belongs_to,Pnode_belongs_to
use module_timing

implicit none

integer(int32) :: inode,i,iV,iP

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{compute\_belongs}
!@@ This subroutine allocates and fills the {\sl U,V,Wnode\_belongs\_to} and {\sl Pnode\_belongs\_to} 
!@@ arrays. For a given Unode {\sl ip},
!@@ {\sl Unode\_belongs\_to(1,ip)} is the number of elements that {\sl ip} belongs to.
!@@ Furthermore, {\sl Unode\_belongs\_to(2:9,ip)} is the actual list of elements.
!==================================================================================================!

call tic() 

!==============================================================================!



write(*,'(a,i5)') shift//'NV=',NV

allocate(Vnode_belongs_to(5,NV)) ; Vnode_belongs_to=0

do iel=1,nel
   do i=1,mV
      inode=mesh(iel)%iconV(i)
      Vnode_belongs_to(1,inode)=Vnode_belongs_to(1,inode)+1
      if (Vnode_belongs_to(1,inode)>5) then
         print *, 'compute_belongs: Vnode_belongs_to array too small'
         stop
      end if
      Vnode_belongs_to(1+Vnode_belongs_to(1,inode),inode)=iel
   end do
end do

!----------------------------------------------------------

allocate(Pnode_belongs_to(5,NV)) ; Pnode_belongs_to=0

do iel=1,nel
   !print *,'elt:',iel
   do i=1,mP
      inode=mesh(iel)%iconP(i)
      !print *,'->',inode
      Pnode_belongs_to(1,inode)=Pnode_belongs_to(1,inode)+1
      if (Pnode_belongs_to(1,inode)>5) then
         print *, 'compute_belongs: Pnode_belongs_to array too small'
         stop
      end if
      Pnode_belongs_to(1+Pnode_belongs_to(1,inode),inode)=iel
   end do
end do

!-----------------------------------------------------------------------------!

if (debug) then
write(2345,'(a)') limit//'compute_belongs'//limit
do iV=1,NV
write(2345,'(a,i6,a,i2,a,4i6)') 'V node',iV,' belongs to ',Vnode_belongs_to(1,iV),' elts:  ',Vnode_belongs_to(2:,iV)
end do
do iP=1,NP
write(2345,'(a,i6,a,i2,a,4i6)') 'P node',iP,' belongs to ',pnode_belongs_to(1,iP),' elts:  ',pnode_belongs_to(2:,iP)
end do
end if

!==============================================================================!

call toc('_____compute_belongs')

end subroutine

!==================================================================================================!
!==================================================================================================!
