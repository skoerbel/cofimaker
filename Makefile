####################################
# begin system-specific settings
####################################
#
######################
# begin use on HP-Notebook or meta
######################
#
FC = gfortran
FCFLAGS = -O3 -fbounds-check -fbacktrace -fdump-core -ffpe-trap='invalid','zero','overflow','underflow','precision','denormal' -std=legacy
LIBS = -llapack
#
######################
# end use on michel 
######################
#
#
######################
# begin use on omega 
######################
#
#FC = gfortran
#FCFLAGS = -O3 -fbounds-check -fbacktrace -fdump-core -ffpe-trap='invalid','zero','overflow','underflow','precision','denormal' -std=gnu
#LIBS = -L ~/.local/lib -llapack
#
######################
# end use on omega 
######################
#
######################
# begin use in Halle & Fionn & lonsdale
######################
#
#FC=ifort
#FCFLAGS=-O3 -align none -mkl 
#
######################
# end use in Halle 
######################
#
#######################
# begin use on supermuc (LRZ)
# ######################
##
#FC=ifort
#FCFLAGS=-O3 -align none -mkl 
#LIBS=-L/lrz/sys/intel/compiler/composer_xe_2015.5.223/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
##
#######################
## end use on supermuc (LRZ) 
#######################
#
#
#LIBS = -L ~/.local/lib/liblapack.a ~/.local/lib/libblas.a ~/.local/lib/liblapack.so.3
#
##### f95 #############################
#FC = f95
#FCFLAGS = -O3 
#
####################################
# end system-specific settings
####################################
#
EXECUTABLE=cofima
OBJECTS= defs.o combi.o integr.o deriv.o misc.o linalg.o ipol.o transform.o writecoords.o readcoords.o postproc.o convert.o crysanpero.o for_lammps.o for_mbpp.o for_gulp.o for_dlpoly.o replace.o cutmod.o distri.o aver.o md.o mult.o cluster.o gaussians.o madelung.o symmetry.o help_cofima.o cofima.o

all:	$(OBJECTS)
	$(FC) $(FCFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIBS)
	
cofima: $(OBJECTS)	
	$(FC) $(FCFLAGS) -o $@ $(OBJECTS) $(LIBS)

clean:
	    rm -rf *o *.mod $(EXECUTABLE)

# dependencies:
misc.o: defs.f combi.f
combi.o: defs.f
deriv.o: defs.f
integr.o: combi.f defs.f
ipol.o: defs.f linalg.f
for_lammps.o: defs.f misc.f
for_mbpp.o: defs.f misc.f
for_gulp.o: defs.f misc.f
for_dlpoly.o: defs.f
crysanpero.o: defs.f
transform.o: defs.f  
writecoords.o: defs.f transform.f
readcoords.o: defs.f writecoords.f misc.f
mult.o: defs.f readcoords.f writecoords.f
postproc.o: defs.f misc.f readcoords.f linalg.o transform.o
convert.o: defs.f readcoords.f writecoords.f
cutmod.o: defs.f readcoords.f writecoords.f misc.f replace.f
distri.o: defs.f readcoords.f writecoords.f
aver.o: readcoords.f
linalg.o: defs.f combi.f misc.f
md.o: defs.f misc.f
replace.o:defs.f readcoords.f misc.f writecoords.f
cluster.o: misc.f replace.f
gaussians.o: defs.f integr.f
madelung.o: defs.f
symmetry.o: defs.f misc.f transform.f writecoords.f
cofima.o: defs.f convert.f deriv.f integr.f ipol.f linalg.f misc.f postproc.f crysanpero.f for_lammps.f for_mbpp.f for_gulp.f for_dlpoly.f transform.f mult.f cutmod.f distri.f aver.f md.f replace.f cluster.o writecoords.f readcoords.f gaussians.f madelung.f help_cofima.f symmetry.f


#%.o:    %.f
#	$(FC) -c $(FCFLAGS) $<

.o:    
	$(FC) -c $(FCFLAGS) $<

