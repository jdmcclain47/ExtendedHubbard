include ../makefile.inc

CC=mpic++

CFLAGS=-w

INC = -I./../lib/libcommon -I./../ -I${EIGENDIR} #-I/usr/local/include/eigen3

LIB = \
-L./../lib/libcommon \
-L/usr/lib -lcommon \
-lm -lgfortran 
#-L/lib/libblas \
#-llapack -lblas 

SRCS = ao_ewald_io.cpp scf_init.cpp scf_dens_gamma.cpp scf_dens_restart.cpp scf_dens_kpoint.cpp scf_diis.cpp \
       scf.cpp scf_fill_matrix.cpp scf_fourier.cpp scf_make_matrices.cpp scf_print.cpp scf_write_mo.cpp \
       ppp.cpp ao_ints.cpp ao_atomic_dist.cpp ao_make_full_gamma.cpp moint_kpoint_twobody.cpp \
       ao_kernel_types.cpp ao_hopping.cpp ao_ewald.cpp ao_integral_cut_ws.cpp ao_non_coulomb_kern.cpp create_header.cpp \
       davidson.cpp davidson_gamma.cpp mointegrals.cpp moint_write_file.cpp moint_read_file.cpp binary_search.cpp \
       mp2.cpp sorter.cpp main.cpp davidson_new.cpp

FSRCS = write_mo_file.f

OBJS = $(SRCS:.cpp=.o)
FOBJS = $(FSRCS:.f=.o)

MAIN = a.out

all: $(MAIN)

$(MAIN): $(OBJS) $(FOBJS)
	$(CC) $(CFLAGS) $(INC) -o $(MAIN) $(OBJS) $(FOBJS) $(LIB);
	cp $(MAIN) ../tests;

.f.o:
	gfortran -c $< -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INC) $(LIB) -c $< -o $@

clean:
	rm *.o
