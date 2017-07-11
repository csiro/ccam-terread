#Makefile for GLOBPE on NEC SX4   f90 (from mrd)

# Specify any desired preprocessor or compiler flags here; -R2 for .L files.

ifneq ($(CUSTOM),yes)
FC = ifort
XFLAGS = -O -assume byterecl -fp-model precise -traceback
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf -lnetcdff
INC = -I $(NETCDF_ROOT)/include
PPFLAG90 = -fpp
PPFLAG77 = -fpp
endif

ifeq ($(GFORTRAN),yes)
FC = gfortran
XFLAGS = -O2 -mtune=native -march=native -I $(NETCDF_ROOT)/include
PPFLAG90 = -x f95-cpp-input
PPFLAG77 = -x f77-cpp-input
endif

OBJT = terread.o read_ht.o read1km.o read250.o read10km.o amap.o rwork.o ncdf.o \
       ccinterp.o latltoij_m.o setxyz_m.o xyzinfo_m.o newmpar_m.o \
       indices_m.o parm_m.o precis_m.o ind_m.o jimco_m.o jimcc_m.o \
       jim_utils.o nfft_m.o netcdf_m.o stacklimit.o

terread :$(OBJT)
	$(FC) $(XFLAGS) $(OBJT) $(LIBS) -o terread

clean:
	rm *.o *.mod core terread

# This section gives the rules for building object modules.

.SUFFIXES:.f90

stacklimit.o: stacklimit.c
	cc -c stacklimit.c
version.h: FORCE
	rm -f brokenver tmpver
	echo "      character(len=*), parameter :: version ='TERREAD r'" > brokenver
	echo "      character(len=*), parameter :: version ='TERREAD r`svnversion .`'" > tmpver
	grep exported tmpver || grep Unversioned tmpver || cmp tmpver brokenver || cmp tmpver version.h || mv tmpver version.h
FORCE:

.f90.o:
	$(FC) -c $(XFLAGS) $(INC) $(PPFLAG90) $<
.f.o:
	$(FC) -c $(XFLAGS) $(INC) $(PPFLAG77) $<

# Remove mod rule from Modula 2 so GNU make doesn't get confused
%.o : %.mod

terread.o : ccinterp.o rwork.o netcdf_m.o version.h
read250.o : ccinterp.o rwork.o netcdf_m.o
read1km.o : ccinterp.o rwork.o netcdf_m.o
read10km.o : ccinterp.o rwork.o
read_ht.o : netcdf_m.o
ccinterp.o : ccinterp.f90 setxyz_m.o xyzinfo_m.o latltoij_m.o newmpar_m.o indices_m.o precis_m.o
latltoij_m.o : latltoij_m.f90 xyzinfo_m.o newmpar_m.o precis_m.o
setxyz_m.o : setxyz_m.f90 newmpar_m.o indices_m.o parm_m.o precis_m.o ind_m.o xyzinfo_m.o jimco_m.o jimcc_m.o 
xyzinfo_m.o : xyzinfo_m.f90 precis_m.o
newmpar_m.o : newmpar_m.f90 
precis_m.o : precis_m.f90
indices_m.o : indices_m.f90
parm_m.o : parm_m.f90 precis_m.o 
ind_m.o : ind_m.f90 newmpar_m.o 
jimcc_m.o : jimcc_m.f90 parm_m.o precis_m.o 
jimco_m.o : jimco_m.f90 precis_m.o jim_utils.o nfft_m.o 
jim_utils.o : jim_utils.f90 precis_m.o 
nfft_m.o : nfft_m.f90 precis_m.o
ncdf.o : netcdf_m.o


