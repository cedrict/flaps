!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine geometry_setup

use iso_fortran_env
use module_parameters
use module_mesh 
use module_constants,only: pi,eps
!use module_arrays
use module_timing

implicit none

integer(int32):: nnr,nnt,i,j,counter,k
real(real64):: r,theta

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{geometry\_setup}
!@@
!==================================================================================================!

call tic() 

!==============================================================================!


nnr=2*nelr+1
nnt=2*nelt+1

NV=nnr*nnt

NP=(nelr+1)*(nelt+1)

NfemV=ndofV*NV
NfemP=NP

Nfem=NfemV+NfemP

write(*,'(a,i6)') shift//'nelr =',nelr
write(*,'(a,i6)') shift//'nelt =',nelt
write(*,'(a,i6)') shift//'nel =',nel
write(*,'(a,i6)') shift//'NV =',NV
write(*,'(a,i6)') shift//'NP =',NP

!----------------------------------------------------------
! velocity nodes
!----------------------------------------------------------
! 473
! 896
! 152
!----------------------------------------------------------

counter=0
do j=1,nelr
   do i=1,nelt
      counter=counter+1
      mesh(counter)%iconV(1)=(i-1)*2+1+(j-1)*2*nnt
      mesh(counter)%iconV(2)=(i-1)*2+2+(j-1)*2*nnt
      mesh(counter)%iconV(3)=(i-1)*2+3+(j-1)*2*nnt
      mesh(counter)%iconV(4)=(i-1)*2+1+(j-1)*2*nnt+nnt
      mesh(counter)%iconV(5)=(i-1)*2+2+(j-1)*2*nnt+nnt
      mesh(counter)%iconV(6)=(i-1)*2+3+(j-1)*2*nnt+nnt
      mesh(counter)%iconV(7)=(i-1)*2+1+(j-1)*2*nnt+nnt*2
      mesh(counter)%iconV(8)=(i-1)*2+2+(j-1)*2*nnt+nnt*2
      mesh(counter)%iconV(9)=(i-1)*2+3+(j-1)*2*nnt+nnt*2

      ! node 1
      r=R1+(j-1)*(R2-R1)/nelr ; theta=pi/2-(i-1)*pi/nelt
      mesh(counter)%xV(1)=r*cos(theta)
      mesh(counter)%zV(1)=r*sin(theta)
      mesh(counter)%rV(1)=r
      mesh(counter)%thetaV(1)=theta

      ! node 5
      r=R1+(j-1)*(R2-R1)/nelr ; theta=pi/2-(i-1+0.5)*pi/nelt
      mesh(counter)%xV(5)=r*cos(theta)
      mesh(counter)%zV(5)=r*sin(theta)
      mesh(counter)%rV(5)=r
      mesh(counter)%thetaV(5)=theta

      ! node 2
      r=R1+(j-1)*(R2-R1)/nelr ; theta=pi/2-(i-1+1)*pi/nelt
      mesh(counter)%xV(2)=r*cos(theta)
      mesh(counter)%zV(2)=r*sin(theta)
      mesh(counter)%rV(2)=r
      mesh(counter)%thetaV(2)=theta

      ! node 8
      r=R1+(j-1+0.5)*(R2-R1)/nelr ; theta=pi/2-(i-1)*pi/nelt
      mesh(counter)%xV(8)=r*cos(theta)
      mesh(counter)%zV(8)=r*sin(theta)
      mesh(counter)%rV(8)=r
      mesh(counter)%thetaV(8)=theta

      ! node 9
      r=R1+(j-1+0.5)*(R2-R1)/nelr ; theta=pi/2-(i-1+0.5)*pi/nelt
      mesh(counter)%xV(9)=r*cos(theta)
      mesh(counter)%zV(9)=r*sin(theta)
      mesh(counter)%rV(9)=r
      mesh(counter)%thetaV(9)=theta

      ! node 6
      r=R1+(j-1+0.5)*(R2-R1)/nelr ; theta=pi/2-(i-1+1)*pi/nelt
      mesh(counter)%xV(6)=r*cos(theta)
      mesh(counter)%zV(6)=r*sin(theta)
      mesh(counter)%rV(6)=r
      mesh(counter)%thetaV(6)=theta

      ! node 4
      r=R1+(j-1+1)*(R2-R1)/nelr ; theta=pi/2-(i-1)*pi/nelt
      mesh(counter)%xV(4)=r*cos(theta)
      mesh(counter)%zV(4)=r*sin(theta)
      mesh(counter)%rV(4)=r
      mesh(counter)%thetaV(4)=theta

      ! node 7
      r=R1+(j-1+1)*(R2-R1)/nelr ; theta=pi/2-(i-1+0.5)*pi/nelt
      mesh(counter)%xV(7)=r*cos(theta)
      mesh(counter)%zV(7)=r*sin(theta)
      mesh(counter)%rV(7)=r
      mesh(counter)%thetaV(7)=theta

      ! node 3
      r=R1+(j-1+1)*(R2-R1)/nelr ; theta=pi/2-(i-1+1)*pi/nelt
      mesh(counter)%xV(3)=r*cos(theta)
      mesh(counter)%zV(3)=r*sin(theta)
      mesh(counter)%rV(3)=r
      mesh(counter)%thetaV(3)=theta

      !do k=1,mV
      !   print *,mesh(counter)%xV(k),mesh(counter)%zV(k)
      !end do

      mesh(counter)%xc=mesh(counter)%xV(9)
      mesh(counter)%zc=mesh(counter)%zV(9)

   end do
end do

do iel=1,nel
   do i=1,mV
      if (mesh(iel)%rV(i)/R2>1-eps) then
         mesh(iel)%hullV(i)=.True.
         mesh(iel)%surfaceV(i)=.True.
         mesh(iel)%surface=.True.
      end if
      if (mesh(iel)%rV(i)/R1<1+eps) then
         mesh(iel)%hullV(i)=.True.
         mesh(iel)%cmbV(i)=.True.
         mesh(iel)%cmb=.True.
      end if
      if (mesh(iel)%xV(i)/R1<eps) mesh(iel)%hullV(i)=.True.
   end do
end do


!----------------------------------------------------------
! pressure nodes
!----------------------------------------------------------
! 4-3
! | |
! 1-2
!----------------------------------------------------------

counter=0
do j=1,nelr
   do i=1,nelt
      counter=counter+1
      mesh(counter)%iconP(1)=i+(j-1)*(nelt+1)  
      mesh(counter)%iconP(2)=i+1+(j-1)*(nelt+1)
      mesh(counter)%iconP(3)=i+1+j*(nelt+1)
      mesh(counter)%iconP(4)=i+j*(nelt+1)

      ! node 1
      r=R1+(j-1)*(R2-R1)/nelr ; theta=pi/2-(i-1)*pi/nelt
      mesh(counter)%xP(1)=r*cos(theta)
      mesh(counter)%zP(1)=r*sin(theta)

      ! node 2
      r=R1+(j-1)*(R2-R1)/nelr ; theta=pi/2-(i-1+1)*pi/nelt
      mesh(counter)%xP(2)=r*cos(theta)
      mesh(counter)%zP(2)=r*sin(theta)

      ! node 3
      r=R1+(j-1+1)*(R2-R1)/nelr ; theta=pi/2-(i-1+1)*pi/nelt
      mesh(counter)%xP(3)=r*cos(theta)
      mesh(counter)%zP(3)=r*sin(theta)

      ! node 4
      r=R1+(j-1+1)*(R2-R1)/nelr ; theta=pi/2-(i-1)*pi/nelt
      mesh(counter)%xP(4)=r*cos(theta)
      mesh(counter)%zP(4)=r*sin(theta)

   end do
end do

!----------------------------------------------------------

if (debug) then

open(unit=123,file='OUTPUT/iconV.dat',status='replace')
do iel=1,nel
   write(123,'(a)') '----------------------------'
   write(123,'(a,i4,a)') '---element #',iel,' -----------'
   write(123,'(a)') '----------------------------'
   write(123,'(a,i8,a,2f20.10)') ' node 1 ', mesh(iel)%iconV(1),' at pos. ',mesh(iel)%xV(1),mesh(iel)%zV(1)
   write(123,'(a,i8,a,2f20.10)') ' node 2 ', mesh(iel)%iconV(2),' at pos. ',mesh(iel)%xV(2),mesh(iel)%zV(2)
   write(123,'(a,i8,a,2f20.10)') ' node 3 ', mesh(iel)%iconV(3),' at pos. ',mesh(iel)%xV(3),mesh(iel)%zV(3)
   write(123,'(a,i8,a,2f20.10)') ' node 4 ', mesh(iel)%iconV(4),' at pos. ',mesh(iel)%xV(4),mesh(iel)%zV(4)
   write(123,'(a,i8,a,2f20.10)') ' node 5 ', mesh(iel)%iconV(5),' at pos. ',mesh(iel)%xV(5),mesh(iel)%zV(5)
   write(123,'(a,i8,a,2f20.10)') ' node 6 ', mesh(iel)%iconV(6),' at pos. ',mesh(iel)%xV(6),mesh(iel)%zV(6)
   write(123,'(a,i8,a,2f20.10)') ' node 7 ', mesh(iel)%iconV(7),' at pos. ',mesh(iel)%xV(7),mesh(iel)%zV(7)
   write(123,'(a,i8,a,2f20.10)') ' node 8 ', mesh(iel)%iconV(8),' at pos. ',mesh(iel)%xV(8),mesh(iel)%zV(8)
   write(123,'(a,i8,a,2f20.10)') ' node 9 ', mesh(iel)%iconV(9),' at pos. ',mesh(iel)%xV(9),mesh(iel)%zV(9)
end do
close(123)

open(unit=123,file='OUTPUT/iconP.dat',status='replace')
do iel=1,nel
   write(123,'(a)') '----------------------------'
   write(123,'(a,i4,a)') '---element #',iel,' -----------'
   write(123,'(a)') '----------------------------'
   write(123,'(a,i8,a,2f20.10)') ' node 1 ', mesh(iel)%iconP(1),' at pos. ',mesh(iel)%xP(1),mesh(iel)%zP(1)
   write(123,'(a,i8,a,2f20.10)') ' node 2 ', mesh(iel)%iconP(2),' at pos. ',mesh(iel)%xP(2),mesh(iel)%zP(2)
   write(123,'(a,i8,a,2f20.10)') ' node 3 ', mesh(iel)%iconP(3),' at pos. ',mesh(iel)%xP(3),mesh(iel)%zP(3)
   write(123,'(a,i8,a,2f20.10)') ' node 4 ', mesh(iel)%iconP(4),' at pos. ',mesh(iel)%xP(4),mesh(iel)%zP(4)
end do
close(123)

end if

!==============================================================================!

call toc('______geometry_setup')

end subroutine

!==================================================================================================!
!==================================================================================================!
