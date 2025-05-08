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

integer(int32):: ip,jp,i,j,k,i1,i2,j1,j2,k1,k2,nsees,nz,inode,jnode,idof,jdof,irow,jcol
logical, dimension(:), allocatable :: alreadyseen

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{setup\_matrix\_K}
!@@ This subroutine allocates arrays ia, ja, and mat of csrK, and builds arrays ia and ja.
!==================================================================================================!

call tic() 

!==============================================================================!

csrK%N=NfemV

write(*,'(a,i11,a)') shift//'csrK%N=',csrK%N,' '

!----------------------------------------------------------
! compute NZ
! I here only compute NZ assuming there is only 1 dof 
! per node. I then multiply this number by ndofV^2 to 
! get the correct number.
!----------------------------------------------------------

call cpu_time(t3)
allocate(alreadyseen(NV))
NZ=0

do inode=1,NV
   alreadyseen=.false.
   do k=1,Vnode_belongs_to(1,inode)
      iel=Vnode_belongs_to(1+k,inode)
      do i=1,mV
         jnode=mesh(iel)%iconV(i)
         if (.not.alreadyseen(jnode)) then
            NZ=NZ+1
            alreadyseen(jnode)=.true.
         end if
      end do
   end do
end do

write(*,'(a,i7)') shift//'NZ (ndof=1)=',NZ

NZ=NZ*ndofV**2

csrK%NZ=NZ

write(*,'(a,i7)') shift//'NZ (ndof=2)=',NZ

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

do inode=1,NV                            ! loop over V nodes
   do idof=1,ndofV                       ! loop over dofs on node
      irow=ndofV*(inode-1)+idof          ! row in K matrix
      nsees=0
      alreadyseen(:)=.false.
      !print *,'--------------------'
      do k=1,Vnode_belongs_to(1,inode)   ! loop over elts to which inode belongs
         iel=Vnode_belongs_to(1+k,inode)
         do i=1,mV                       ! loop over nodes of element iel
            jnode=mesh(iel)%iconV(i)
            !print *,jnode
            do jdof=1,ndofV
               jcol=ndofV*(jnode-1)+jdof ! column in K matrix

               if (.not.alreadyseen(jcol)) then
                  NZ=NZ+1
                  csrK%ja(NZ)=jcol
                  nsees=nsees+1
                  alreadyseen(jcol)=.true.
               end if
            end do

         end do
      end do
      !print *,'irow=',irow,'nsees:',nsees
      csrK%ia(irow+1)=csrK%ia(irow)+nsees
   end do
end do

!print *,NZ
!print *,csrK%ia
!print *,csrK%ja

deallocate(alreadyseen)    

call cpu_time(t4) ; write(*,'(a,f10.3,a)') shift//'ia,ja time:',t4-t3,'s'

write(*,'(a,i9)' ) shift//'NZ=',NZ
write(*,'(a,2i9)') shift//'csrK%ia',minval(csrK%ia), maxval(csrK%ia)
write(*,'(a,2i9)') shift//'csrK%ja',minval(csrK%ja), maxval(csrK%ja)

!==============================================================================!

call toc('______setup_matrix_K')

end subroutine

!==================================================================================================!
!==================================================================================================!
