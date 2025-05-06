!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

module module_timing

use iso_fortran_env

implicit none

character(len=42), parameter :: shift='                                        |'
character(len=42), parameter :: limit='*****************************************'

real(real64) :: t1,t2,t3,t4

contains

subroutine tic()
  implicit none
  call cpu_time(t1)
end subroutine tic

subroutine toc(thing,nel)
  implicit none
  character(20)::thing
  integer(int32),optional::nel
  call cpu_time(t2)
  if (present(nel)) then
     write(*,'(a,f8.3,a,i8,a)') trim(thing), real(t2-t1),'s',nel,'   |'
  else
     write(*,'(a,f8.2,a)') trim(thing), real(t2-t1),'s           |'
  end if
end subroutine toc

end module

!==================================================================================================!
