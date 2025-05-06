!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! SURGE                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine output_solution_to_vtu

use iso_fortran_env
use module_core
use module_constants,only: zero
use module_parameters,only: nel,iel,debug,export_to_vtu
use module_mesh
use module_timing

implicit none

integer(int32),external:: conv_l1_to_int
integer(int32):: k,cell_type=28
real(real64):: u,v,w

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{output\_solution\_to\_vtu}
!@@ This subroutine exports the solution field and derivatives in vtu format in a file found in 
!@@ the OUTPUT folder. 
!==================================================================================================!

call tic() 

!==============================================================================!

if (.not.export_to_vtu) return

!----------------------------------------------------------

write(*,'(a)') shift//'-> '//'OUTPUT/solution.vtu'

open(unit=123,file='OUTPUT/solution.vtu',status='replace',form='formatted')
write(123,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
write(123,'(a)') '<UnstructuredGrid>'
write(123,'(a,i8,a,i7,a)') '<Piece NumberOfPoints="',mV*nel,'" NumberOfCells="',nel,'">'
!------------------
write(123,'(a)') '<Points>'
write(123,'(a)') '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,'(3es12.4)') mesh(iel)%xV(k),mesh(iel)%yV(k),mesh(iel)%zV(k)
   end do
end do
write(123,'(a)') '</DataArray>'
write(123,'(a)') '</Points>'
!------------------
write(123,'(a)') '<CellData Scalars="scalars">'
!-----
if (debug) then
write(123,*) '<DataArray type="Float32" Name="surface" Format="ascii">'
do iel=1,nel
   write(123,*) conv_l1_to_int(mesh(iel)%surface)
end do
write(123,*) '</DataArray>'
end if
!-----
if (debug) then
write(123,*) '<DataArray type="Float32" Name="cmb" Format="ascii">'
do iel=1,nel
   write(123,*) conv_l1_to_int(mesh(iel)%cmb)
end do
write(123,*) '</DataArray>'
end if
!-----
write(123,'(a)') '<DataArray type="Float32" Name="rho" Format="ascii">'
do iel=1,nel 
   write(123,'(es12.4)') sum(mesh(iel)%rhoq(1:nqel))/nqel
end do
write(123,'(a)') '</DataArray>'
!-----
write(123,'(a)') '<DataArray type="Float32" Name="vol" Format="ascii">'
do iel=1,nel 
   write(123,'(es12.4)') mesh(iel)%volume 
end do
write(123,'(a)') '</DataArray>'
!-----
write(123,*) '</CellData>'
!-------------------------
write(123,*) '<PointData Scalars="scalars">'
!-----
write(123,*) '<DataArray type="Float32" Name="pressure (q)" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,'(es12.4)') mesh(iel)%q(k)
   end do
end do
write(123,*) '</DataArray>'
!-----
if (debug) then
write(123,*) '<DataArray type="Float32" Name="r" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,*) mesh(iel)%rV(k)
   end do
end do
write(123,*) '</DataArray>'
end if
!-----
if (debug) then
write(123,*) '<DataArray type="Float32" Name="theta (polar)" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,*) mesh(iel)%thetaV(k)
   end do
end do
write(123,*) '</DataArray>'
end if
!-----
if (debug) then
write(123,*) '<DataArray type="Int32" Name="hull" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,*) conv_l1_to_int(mesh(iel)%hullV(k))
   end do
end do
write(123,*) '</DataArray>'
end if
!-----
if (debug) then
write(123,*) '<DataArray type="Int32" Name="surface" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,*) conv_l1_to_int(mesh(iel)%surfaceV(k))
   end do
end do
write(123,*) '</DataArray>'
end if
!-----
if (debug) then
write(123,*) '<DataArray type="Int32" Name="cmb" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,*) conv_l1_to_int(mesh(iel)%cmbV(k))
   end do
end do
write(123,*) '</DataArray>'
end if
!-----
if (debug) then
write(123,*) '<DataArray type="Int32" Name="fix_u" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,*) conv_l1_to_int(mesh(iel)%bc_fixV(2*(k-1)+1))
   end do
end do
write(123,*) '</DataArray>'
!-----
write(123,*) '<DataArray type="Int32" Name="fix_w" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,'(i2)') conv_l1_to_int(mesh(iel)%bc_fixV(2*(k-1)+2))
   end do
end do
write(123,*) '</DataArray>'
end if
!-----
if (debug) then
write(123,*) '<DataArray type="Int32" Name="bnd1" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,*) conv_l1_to_int(mesh(iel)%bnd1_elt)
   end do
end do
write(123,*) '</DataArray>'
end if
!-----
if (debug) then
write(123,*) '<DataArray type="Int32" Name="bnd2" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,*) conv_l1_to_int(mesh(iel)%bnd2_elt)
   end do
end do
write(123,*) '</DataArray>'
end if
!-----
if (debug) then
write(123,*) '<DataArray type="Int32" Name="bnd3" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,*) conv_l1_to_int(mesh(iel)%bnd3_elt)
   end do
end do
write(123,*) '</DataArray>'
end if
!-----
if (debug) then
write(123,*) '<DataArray type="Int32" Name="bnd4" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,*) conv_l1_to_int(mesh(iel)%bnd4_elt)
   end do
end do
write(123,*) '</DataArray>'
end if
!-----
!if (debug) then
!write(123,*) '<DataArray type="Float32" Name="analytical pressure" Format="ascii">'
!do iel=1,nel
!   do k=1,m
!      call analytical_solution(mesh(iel)%x(k),mesh(iel)%y(k),mesh(iel)%z(k),Panal,u,v,w)
!      write(123,'(es11.3)') Panal
!   end do
!end do
!write(123,*) '</DataArray>'
!end if
!-----
!if (debug) then
!write(123,*) '<DataArray type="Float32" NumberOfComponents="3" Name="analytical velocity" Format="ascii">'
!do iel=1,nel
!   do k=1,m
!      call analytical_solution(mesh(iel)%x(k),mesh(iel)%y(k),mesh(iel)%z(k),Panal,u,v,w)
!      write(123,'(3es11.3)') u,v,w
!   end do
!end do
!write(123,*) '</DataArray>'
!end if
!-----
write(123,*) '<DataArray type="Float32" NumberOfComponents="3" Name="velocity" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,'(3es12.4)') mesh(iel)%u(k),mesh(iel)%v(k),mesh(iel)%w(k)
   end do
end do
write(123,*) '</DataArray>'

!-----
write(123,*) '<DataArray type="Float32" NumberOfComponents="3" Name="normal vector" Format="ascii">'
do iel=1,nel
   do k=1,mV
      write(123,'(3es12.4)') mesh(iel)%nx(k),0.,mesh(iel)%nz(k)
   end do
end do
write(123,*) '</DataArray>'







!-----
write(123,*) '</PointData>'
!------------------
write(123,*) '<Cells>'
!-----
write(123,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
do iel=1,nel
   write(123,*) ( (iel-1)*mV+k-1,k=1,mV) 
end do
write(123,*) '</DataArray>'
!-----
write(123,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
write(123,*) (iel*mV,iel=1,nel)
write(123,*) '</DataArray>'
!-----
write(123,*) '<DataArray type="Int32" Name="types" Format="ascii">'
write(123,*) (cell_type,iel=1,nel)
write(123,*) '</DataArray>'
!-----
write(123,*) '</Cells>'
!------------------
write(123,*) '</Piece>'
write(123,*) '</UnstructuredGrid>'
write(123,*) '</VTKFile>'
close(123)

!==============================================================================!

call toc('___output_sol_to_vtu')

end subroutine

!==================================================================================================!
!==================================================================================================!

pure function conv_l1_to_int(b)
implicit none
logical(1), intent(in) :: b
integer conv_l1_to_int
if (b) then
conv_l1_to_int=1
else
conv_l1_to_int=0
end if
end function
