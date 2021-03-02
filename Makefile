
#### Makefile for Cur+ ####

# CC= cc
CFLAG= 

# FC= gfortran
#FFLAG= -w -fbounds-check -O2
FFLAG= -w -O2

###------------------ Amber netCDF input support ------------------###
## Amber netCDF support requires an external installation of the
## netCDF Fortran Libraries. These can be installed using your
## favourite package manager, where available (e.g. on Debian:
## "apt-get install netcdf-fortran", on MacOS X: "port install
## netcdf-fortran").  Otherwise, the netCDF Fortran Libraries can
## be installed from the source tarball (at
## http://www.unidata.ucar.edu/downloads/netcdf).  Finally, the
## netCDF Fortran Libraries can also be found within
## installations of the AmberTools package (http://ambermd.org).
## 
## If you have AmberTools installed (version >=13), please set
## NETCDF to the location of your AmberTools installation; e.g.:
# NETCDF=$(AMBERHOME)/AmberTools
NETCDF=$(AMBERHOME)
##
## If you installed netcdf-fortran using macports on MACOSX,
## set NETCDF to the location of your macports tree; e.g.:
## NETCDF=/opt/local
##
## If you installed netcdf-fortran using the official Unidata 
## tarball, set NETCDF to the location of your installation,
## e.g.:
## NETCDF=/usr/local
##

###------------------ Gromacs XTC input support -------------------###
## To enable XTC support, uncomment the following line:
XTC= yes

## Please do not edit below this line.
###----------------------------------------------------------------###
ifdef NETCDF
FFLAG_NC= -I$(NETCDF)/include/ -DNETCDF
ifneq ("$(wildcard $(NETCDF)/lib/libnetcdff.a)","")
libnetcdf= netcdff
else
libnetcdf= netcdf
endif
LDFLAG_NC= -L$(NETCDF)/lib -l$(libnetcdf)
endif

ifdef XTC
XTC_OBJ= f77_molfile.o
FFLAG_XTC= -DXTC
LDFLAG_XTC= -L. -lmolfile_plugin -lstdc++
endif

OBJ= aacur.o axis.o axref.o backbo.o bisection.o dotdelta.o eigen.o \
findaxis.o getdate.o ionbld.o ionpar.o input.o intaxe.o intop.o \
locate.o lsfit.o manta.o ncerror.o nml.o \
params.o pdbout.o curpar.o screw.o setup.o smooth.o torp.o xtcerror.o $(XTC_OBJ)

all: Cur+ Canal Canion

Cur+ : $(OBJ) 
	$(FC) $(OBJ) $(LDFLAG_NC) $(LDFLAG_XTC) -o Cur+
Canal : canal.o
	$(FC) $(FFLAG) canal.o -o Canal
Canion : canion.o
	$(FC) $(FFLAG) canion.o -o Canion

.SUFFIXES : .o .c .f
.c.o :
	$(CC) -cpp $(CFLAG) -c $<

.f.o :
	$(FC) -cpp $(FFLAG) $(FFLAG_NC) $(FFLAG_XTC) -c $<

clean: 
	rm -f Cur+ Canal Canion *.o *.a

f77_molfile.o: f77_molfile.c plugins
	$(CC) $(CFLAG) -Iinclude -c $<

plugins: gromacsplugin.C
	$(CC) -D"STATIC_PLUGIN" -Iinclude -D"VMDPLUGIN=molfile_gromacsplugin" -c gromacsplugin.C -o gromacsplugin-s.o
	ar cr libmolfile_plugin.a gromacsplugin-s.o

