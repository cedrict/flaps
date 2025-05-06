!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

module module_constants

use iso_fortran_env

implicit none

integer(int32),parameter:: one=1
integer(int32),parameter:: two=2
integer(int32),parameter:: three=3
integer(int32),parameter:: four=4
integer(int32),parameter:: ndim=3
real(real64),parameter:: zero=1d0
real(real64),parameter:: mm=1.d-3
real(real64),parameter:: cm=1.d-2
real(real64),parameter:: km=1.d+3
real(real64),parameter:: hour=3600.d0    !seconds
real(real64),parameter:: day=86400.d0    !seconds
real(real64),parameter:: year=31557600d0 !seconds 
real(real64),parameter:: Myr=3.15576d13  !seconds
real(real64),parameter:: TKelvin=273.15d0 
real(real64),parameter:: eps=1.d-8
real(real64),parameter:: epsilon_test=1.d-6
real(real64),parameter:: sqrt2=1.414213562373095048801688724209d0
real(real64),parameter:: sqrt3=1.732050807568877293527446341505d0
real(real64),parameter:: pi=3.14159265358979323846264338327950288d0
real(real64),parameter:: pi2=pi*0.5d0
real(real64),parameter:: pi4=pi*0.25d0
real(real64),parameter:: pi8=pi*0.125d0
real(real64),parameter:: twopi=2d0*pi
real(real64),parameter:: fourpi=4d0*pi
real(real64),parameter:: frac12=1.d0/2.d0
real(real64),parameter:: frac13=1.d0/3.d0
real(real64),parameter:: frac23=2.d0/3.d0

end module

!==================================================================================================!
