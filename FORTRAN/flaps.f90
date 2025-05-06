program flaps

use iso_fortran_env
use module_mesh
use module_parameters
use module_timing

implicit none

!==========================================================

print *,compiler_version()
print *,compiler_options()

!==========================================================
call header
call experiment_declare_main_parameters
!call read_command_line_options

call allocate_mesh_array
call geometry_setup
call quadrature_setup
call test_basis_functions
call compute_belongs
call setup_matrix_K
call sanity_check
call compute_normals
call experiment_define_bc_V
call spacer
call write_params
call spacer
call assign_material_properties_to_qpoints
call output_qpoints
call make_matrix_stokes
call output_solution_to_vtu
call output_solution_to_ascii
call postprocess

!==========================================================

call footer

end program
