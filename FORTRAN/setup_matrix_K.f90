!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine setup_matrix_K

use iso_fortran_env
use module_core
use module_parameters,only: NfemV,NV,iel
use module_sparse,only: csrK
use module_mesh
use module_arrays,only: Vnode_belongs_to
use module_timing

implicit none

integer(int32):: ip,jp,i,j,k,i1,i2,j1,j2,k1,k2,nsees,nz
logical, dimension(:), allocatable :: alreadyseen

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{setup\_matrix\_K}
!@@ This subroutine allocates arrays ia, ja, and mat of csrK, and builds arrays ia and ja.
!==================================================================================================!

call tic() 

!==============================================================================!

csrK%N=NfemV

write(*,'(a,i11,a)') shift//'csrK%N       =',csrK%N,' '

!----------------------------------------------------------
! compute NZ
!----------------------------------------------------------

call cpu_time(t3)
allocate(alreadyseen(NV))
NZ=0

!xx,xy,xz
do ip=1,NV
   alreadyseen=.false.
   do k=1,Vnode_belongs_to(1,ip)
      iel=Vnode_belongs_to(1+k,ip)
      !-------- 
      do i=1,mV
         jp=mesh(iel)%iconV(i)
         if (.not.alreadyseen(jp)) then
            NZ=NZ+1
            alreadyseen(jp)=.true.
         end if
      end do
   end do
end do

print *,'1dof NZ=',NZ

NZ=NZ*ndofV**2

csrK%NZ=NZ

print *,'NZ=',NZ

deallocate(alreadyseen)
call cpu_time(t4)


!----------------------------------------------------------
! fill arrays ia,ja
!----------------------------------------------------------

write(*,'(a,i11,a,f7.3,a)') shift//'csrK%NZ      =',csrK%NZ,' | ',t4-t3,'s'

allocate(csrK%ia(csrK%N+1)) ; csrK%ia=0 
allocate(csrK%ja(csrK%NZ))  ; csrK%ja=0 
allocate(csrK%mat(csrK%NZ)) ; csrK%mat=0 

call cpu_time(t3)
allocate(alreadyseen(2*NV))
NZ=0
csrK%ia(1)=1

do ip=1,NV
   nsees=0
   alreadyseen=.false.
   do k=1,Vnode_belongs_to(1,ip)
      iel=Vnode_belongs_to(1+k,ip)
      !-------- 
      do i=1,mV
         jp=mesh(iel)%iconV(i)
         if (.not.alreadyseen(jp)) then
            NZ=NZ+1
            csrK%ja(NZ)=jp
            nsees=nsees+1
            alreadyseen(jp)=.true.
         end if
      end do
      !-------- 
      do i=1,mV
         jp=mesh(iel)%iconV(i)+NV
         if (.not.alreadyseen(jp)) then
            NZ=NZ+1
            csrK%ja(NZ)=jp
            nsees=nsees+1
            alreadyseen(jp)=.true.
         end if
      end do
   end do
   csrK%ia(ip+1)=csrK%ia(ip)+nsees
end do

print *,csrK%ia





!==============================================================================!

call toc('______setup_matrix_K')

end subroutine

!==================================================================================================!
!==================================================================================================!
