# Makefile for various platforms
# Execute using Build csh-script only!
# Used together with Perl scripts in SRC/SCRIPT 
# (C) 2005 Marat Khairoutdinov
#------------------------------------------------------------------
# uncomment to disable timers:
#
#NOTIMERS=-DDISABLE_TIMERS
#-----------------------------------------------------------------

# Name of executable
SAM = SAM_$(ADV_DIR)_$(SGS_DIR)_$(RAD_DIR)_$(MICRO_DIR)

## Determine platform 
##PLATFORM := $(shell uname -s)


#ifeq ($(PLATFORM),Linux)


export NETCDF_INCDIR=/mnt/nfs/clustersw/Debian/bullseye/openmpi/4.0.2/usr/netcdf/4.8.1/include
export NETCDF_LIBDIR=/mnt/nfs/clustersw/Debian/bullseye/openmpi/4.0.2/usr/netcdf/4.8.1/lib

# export NETCDFFORTRAN_LIBDIR=/mnt/nfs/clustersw/Debian/bullseye/openmpi/4.0.2/usr/netcdf/4.8.1/lib #CMuller: try

#----------------------------------------------------------------------
#---------- LOC:=IST CLUSTER
#----------------------------------------------------------------------
  FF77 = mpif77
  FF90 = mpif90
  CC   = mpicc
  LD   = mpif90

#LDFLAGS=/mnt/nfs/clustersw/Debian/bullseye/openmpi/4.0.2/usr/netcdf/4.8.1/lib -I/mnt/nfs/clustersw/Debian/bullseye/openmpi/4.0.2/usr/netcdf/4.8.1/include 

#LDFLAGS=-L/mnt/nfs/clustersw/Debian/bullseye/openmpi/4.0.2/usr/netcdf/4.8.1/lib -lmpi -lnetcdff -lnetcdf
  #LDFLAGS = -L${NETCDFC_LIBDIR} -L${NETCDFFORTRAN_LIBDIR} -lmpi -lnetcdff -lnetcdf

  LDFLAGS = -L${NETCDF}/lib -lmpi -lnetcdff -lnetcdf 
  FFLAGS = -fallow-argument-mismatch -c -O2  -I${CPATH}  -I${NETCDF_INCDIR} 
  CFLAGS = -c -O2  -I${CPATH}  -I${NETCDF_INCDIR}  -DFORTRANUNDERSCORE

#  LDFLAGS = -L${LD_LIBRARY_PATH} -L${LD_RUN_PATH} -lmpi -lnetcdff -lnetcdf
#  FFLAGS = -c -O2 -axCOMMON-AVX512,CORE-AVX512 -I${NETCDFC_INCDIR} -I${NETCDFFORTRAN_INCDIR} 
#  CFLAGS = -c -O2 -axCOMMON-AVX512,CORE-AVX512 -DFORTRANUNDERSCORE -I${NETCDFC_INCDIR} -I${NETCDFFORTRAN_INCDIR} 





#----------------------------------------------------------------------
#---------- LOC:=IRENE
#----------------------------------------------------------------------
#
#  FF77 = mpif77
#  FF90 = mpif90
#  CC   = mpicc
#  LD   = mpif90
#  LDFLAGS = -L${NETCDFC_LIBDIR} -L${NETCDFFORTRAN_LIBDIR} -lmpi -lnetcdff -lnetcdf
#  FFLAGS = -c -O2 -axCOMMON-AVX512,CORE-AVX512 -I${NETCDFC_INCDIR} -I${NETCDFFORTRAN_INCDIR} 
#  CFLAGS = -c -O2 -axCOMMON-AVX512,CORE-AVX512 -DFORTRANUNDERSCORE -I${NETCDFC_INCDIR} -I${NETCDFFORTRAN_INCDIR} 




#----------------------------------------------------------------------
#---------- LOC:=CURIE
#----------------------------------------------------------------------
#FF77 = mpif77
#FF90 = mpif90
#CC   = mpicc
#LD   = mpif90
#LDFLAGS = -L${NETCDF_LIBDIR} -lnetcdff -lnetcdf
#FFLAGS = -c -O2 -axAVX,SSE4.1 -I${NETCDF_INCDIR}
#CFLAGS = -c -O2 -axAVX,SSE4.1 -DFORTRANUNDERSCORE -I${NETCDF_INCDIR}

#----------------------------------------------------------------------
#-------  Linux, Intel Compiler (Old example that came with the code)
#----------------------------------------------------------------------
# LIB_MPI = /usr/local/pkg/iopenmpi/lib
# INC_MPI = /usr/local/pkg/iopenmpi/include
# INC_NETCDF = /nfs/user08/marat/local/include
# LIB_NETCDF = /nfs/user08/marat/local/lib
# FF77 = /usr/local/pkg/iopenmpi/bin/mpif90 -c -fixed -extend_source
# FF90 = /usr/local/pkg/iopenmpi/bin/mpif90 -c
# CC = mpicc -c -DLINUX
# FFLAGS = -O3 
# #FFLAGS = -g -ftrapuv -check all
# FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
# LD = /usr/local/pkg/iopenmpi/bin/mpif90
# LDFLAGS = -L${LIB_NETCDF} -lnetcdf


#----------------------------------------------------------------------
#------- LOC := GFDL    OLD VERSION that ran at GFDL 2012-2014
#----------------------------------------------------------------------
#  FF77 = ftn -c -O2 -fpp -fno-alias -stack_temps -safe_cray_ptr -ftz -i_dynamic -assume byterecl -g -Wp,-w
#  FF90 = ftn -c -O2 -fpp -fno-alias -stack_temps -safe_cray_ptr -ftz -i_dynamic -assume byterecl -g -Wp,-w
#  CC   = cc -c -DLINUX
#  LD      = ftn
#  FFLAGS  = -I$(CRAY_NETCDF_DIR)/include
#  LDFLAGS = -L$(CRAY_NETCDF_DIR)/lib


#----------------------------------------------------------------------
#---------- LOC:=JADE version that ran on JADE
#----------------------------------------------------------------------
#  FF77 = ifort 
#  FF90 = ifort
#  CC   = icc
#  LD   = ifort 
#  LDFLAGS = -lnetcdff -lnetcdf -lmpi
#  FFLAGS = -c -O2 -xSSE4.1
#  CFLAGS = -c -O2 -xSSE4.1 -DFORTRANUNDERSCORE 




## endif #<-endif platform=Linux


#----------------------------------------------
# you dont need to edit below this line


#compute the search path
dirs := . $(shell cat Filepath)
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir))) 

.SUFFIXES:
.SUFFIXES: .f .f90 .c .o



all: $(SAM_DIR)/$(SAM)


SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
	$(SAM_SRC)/SCRIPT/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	$(SAM_SRC)/SCRIPT/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES))) 

$(SAM_DIR)/$(SAM): $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)


.f90.o:
	${FF90}  ${FFLAGS} $<
.f.o:
	${FF77}  ${FFLAGS} $<
.c.o:
	${CC}  ${CFLAGS} -I$(SAM_SRC)/TIMING $(NOTIMERS) $<



include Depends



clean: 
	rm ./OBJ/*



