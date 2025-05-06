!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

module module_sparse 

use iso_fortran_env

type compressedrowstorage_sqr    
   integer(int32) :: N                 ! size of square matrix 
   integer(int32) :: NZ                ! number of nonzeros
   integer(int32),allocatable :: ia(:)  
   integer(int32),allocatable :: ja(:)
   integer(int32),allocatable :: ua(:) ! index of diagonal term
   real(real64),dimension(:),allocatable :: mat
end type compressedrowstorage_sqr

type(compressedrowstorage_sqr) csrK 
type(compressedrowstorage_sqr) csrGT 

end module

!==================================================================================================!
