##################################################
# makefile for minimal Anisotropy Toolkit        #
# last edit: 2021-11-25                          #
# edited. by: J. Biteau                          #
##################################################


SYSTEM=$(shell uname)
CC = gcc 
CXX = g++
OPTFLAGS = #-O
DBGFLAGS = -g3
WFLAGS = -D__USE_FIXED_PROTOTYPES__ -Wall
OBJ = ./

.PHONY: clean
#------------------- defs ------------------------------
CFITSIO_INCDIR = $(CFITSIO_DIR)/include
CFITSIO_LIBDIR = $(CFITSIO_DIR)/lib

CHEALPIX_INCDIR = $(HEALPIX_DIR)/include
CHEALPIX_LIBDIR = $(HEALPIX_DIR)/lib

INCDIR = -I$(CFITSIO_INCDIR) \
	-I$(CHEALPIX_INCDIR) \
	$(shell root-config --cflags)

LIBDIR = -L${CHEALPIX_LIBDIR} -lchealpix -lm \
	-L${CFITSIO_LIBDIR} -lcfitsio \
	$(shell root-config --glibs)
#-------------------------------------------------------

#------- alias -----------------------------------------
libobjs = \
	STClibrary.o \
	LatLong-UTMConversion.o \
	healpixmap.o \
	harmotools.o \
	maptools.o \
	projmap.o \
	common.o 

HEADERS = $(patsubst %.o,%.h,$(libobjs))

thelib = libtkit.a
#-------------------------------------------------------

#-------- rules ----------------------------------------
# rules for the library sources
$(OBJ)%.o:%.cc %.h
	$(COMPILE.cc) $(DBGFLAGS) $(OPTFLAGS) $(WFLAGS) $(INCDIR) -o $@ $<
#-------------------------------------------------------

#------- targets ---------------------------------------
all : lib 
lib: $(thelib) 
clean :
	@echo "Deleting library objects, executables and associated objects."
	@/bin/rm -f $(thelib) $(execs) $(libobjs) $(exeobjs)
	@/bin/rm -f *~
	@echo "Done"
#-------------------------------------------------------

#-------- specific rules -------------------------------
$(thelib): $(libobjs)
	@echo "Making "$(thelib)
	@ar r $@ $^
	@ranlib $@
	@echo "Done."
#-------------------------------------------------------

