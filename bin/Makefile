EXEC = hydropad
DEFINES = -DUSEMPI -DTRACE
DEFACC = 
RM = /bin/rm

# Use MPI? [yes|no]
USEMPI = yes
# Use host OpenMP? [yes|no]
USEOMP = no
# Use OpenACC? [yes|no]
USEACC = no

#File .lst
LST = no

# Manual Profiling? [yes|no]
PROFILE = no

F90 = ftn
FFLAGS = -O3 -c
CC = CC
DEBUG = -g
LD = ftn
LFLAGS = 

#If debugging mode, else comment it out
#FFLAGS += -g

FFLAGS += $(DEFINES) $(DEFACC)

#############################################################################
MOD = mod
#############################################################################
# Librairies
LIBS =

# --- CUDA libraries ---
LIBCUDA = 

# Sources directories are searched in this exact order
#####ifeq ($(ACC),no)
VPATH = ../hydropad:../io:../hydro:../nbody:../gravity:../physics
IPAPATH =
#####else
####VPATH = 
####IPAPATH = 
######endif

#############################################################################
# Objects
# Modules
MODS = mpi_inc.o dimension.o scalar.o vector.o matrix.o io_mod.o
# Driver
MAIN = main.o dealloc_arrays.o evolve_sys.o init_sys.o indata_simple.o initialize_hydro.o alloc_arrays.o dealloc_vectors.o alloc_vectors.o savetn.o savetnp1.o expand.o calculate_thermo.o speed.o functions.o means.o exchange_mesh.o exchange_var.o shsearch.o
# I/O
IO = write_var.o write_var_mpi.o
# Hydro
HYDRO = boundary.o ppm_functs.o riemann.o internal.o ppm.o ppmsolver.o interp.o dda.o flatten.o ppm_lxlylz.o total.o
# Gravity
GRAVITY =
# N-body
NBODY = 

# All objects
HYDROPAD = $(MAIN) $(IO) $(HYDRO) $(GRAVITY) $(NBODY)

#Compile 
.SUFFIXES : .f .o .F .F90
#
#.src.o : ${MODS}
#	$(CPP) $(CPPFLAGS) $*.src $*.F90
#	$(F90) $(DEBUG) $(FFLAGS) $*.F90

.F90.o: ${MODS}
	$(F90) $(DEBUG) $(FFLAGS) $<

.f.o: $(MODS)
	$(F90) $(DEBUG) $(FFLAGS) $<

$(EXEC) : $(MODS) $(HYDROPAD)
	$(LD) $(DEBUG) $(LFLAGS) $(MODS) $(HYDROPAD) -o $(EXEC).MPI$(USEMPI).OMP$(USEOMP).ACC$(USEACC) $(LIBS)

clean:
	$(RM) *.o
	${RM} *.mod
	${RM} *.lst