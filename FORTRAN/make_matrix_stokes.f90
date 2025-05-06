!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! SURGE                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine make_matrix_stokes

use iso_fortran_env
use module_parameters
use module_arrays,only: rhs_f,rhs_h
use module_sparse,only: csrK,csrGT
use module_mesh 
use module_timing

implicit none

integer(int32):: inode,jnode,k,k1,k2,m1,m2,ikk,jkk,i1,i2
real(real64):: Kel(mVel,mVel),fel(mVel),Gel(mVel,mP),hel(mP)

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{make\_matrix\_stokes}
!@@ This subroutine builds the pressure linear system for the flow equation. 
!@@ It loops over each element, builds its elemental matrix ${\bm A}_{el}$
!@@ and right hand side $\vec{b}_{el}$, applies boundary conditions, 
!@@ and assembles these into the global matrix csrA and the corresponding right hand side rhs. 
!==================================================================================================!

call tic() 

!==============================================================================!

csrK%mat=0
csrGT%mat=0
rhs_f=0
rhs_h=0

do iel=1,nel

   !call compute_elemental_matrix_P(Kel,bel)

   !call impose_boundary_conditions_P(Kel,bel)

   do k1=1,mV
   do i1=1,ndofV
      inode=mesh(iel)%iconV(k1)
      m1=ndofV*inode+i1
      ikk=ndofV*k1+i1
      do k2=1,mV
      do i2=1,ndofV
         jnode=mesh(iel)%iconV(k2)
         m2=ndofV*inode+i2
         jkk=ndofV*k2+i2
         do k=csrK%ia(inode),csrK%ia(inode+1)-1
            if (csrK%ja(k)==jnode) then
               !csrK%mat(k)=csrK%mat(k)+Kel(k1,k2)
               exit
            end if
         end do
      end do
      end do
      !rhs(inode)=rhs(inode)+bel(k1)
   end do
   end do

end do

write(*,'(a,2es12.4)') shift//'K (m/M)',minval(csrK%mat),maxval(csrK%mat)
write(*,'(a,2es12.4)') shift//'f (m/M)',minval(rhs_f),maxval(rhs_f)

!==============================================================================!

call toc('_______make_matrix_P',nel)

end subroutine

!==================================================================================================!
!==================================================================================================!
