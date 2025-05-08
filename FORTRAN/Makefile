.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

default: code

include Makefile.machine
include Makefile.experiment

OBJECTS = \
module_core.o\
module_constants.o\
module_arrays.o\
module_parameters.o\
module_timing.o\
module_mesh.o\
module_sparse.o\
allocate_mesh_array.o\
basis_functions.o\
compute_belongs.o\
compute_dNdx_dNdz_at_qpt.o\
compute_normals.o\
geometry_setup.o\
header.o\
make_matrix_stokes.o\
output_qpoints.o\
output_solution_to_ascii.o\
output_solution_to_vtu.o\
quadrature_setup.o\
setup_matrix_K.o\
sanity_check.o\
test_basis_functions.o\
flaps.o

.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

.f.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f

code:	$(OBJECTS) $(object_experiment)
	$(F90) $(OPTIONS) $(OBJECTS) experiment.o $(LIBS) -o flaps

clean:	
	rm -f *.o
	rm -f *.mod
	rm -f *.dat 
	rm -f *.out 
	rm -f *.ascii 
	rm -f opla
	rm -f fort.*
	rm -f flaps
	rm -f OUTPUT/*.vtu
	rm -f OUTPUT/*.ascii


cleandata:
	rm -f *.ascii 
	rm -f OUTPUT/*.vtu
	rm -f OUTPUT/*.ascii

	
