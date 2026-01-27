##################################################
# makefile for exposure calculation              #
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
SRC_DIR = ../Utilities

CFITSIO_INCDIR = $(CFITSIO_DIR)/include
CFITSIO_LIBDIR = $(CFITSIO_DIR)/lib

CHEALPIX_INCDIR = $(HEALPIX_DIR)/include
CHEALPIX_LIBDIR = $(HEALPIX_DIR)/lib

INCDIR = -I$(CFITSIO_INCDIR) \
	-I$(CHEALPIX_INCDIR) \
	-I$(SRC_DIR) \
	$(shell root-config --cflags)

LIBDIR = -L${CHEALPIX_LIBDIR} -lchealpix -lm \
	-L${CFITSIO_LIBDIR} -lcfitsio -lMinuit \
	$(shell root-config --glibs)
#-------------------------------------------------------

#------- alias -----------------------------------------
execs = \
	exposure.exe
	 
libobjs = \
	exposure.o
	 	
thelib = $(SRC_DIR)/libtkit.a

#-------------------------------------------------------

#-------- rules ----------------------------------------
# rules for the executable sources
$(OBJ)%.o:%.cc
	$(COMPILE.cc) $(DBGFLAGS) $(OPTFLAGS) $(WFLAGS) $(INCDIR) -o $@ $<
#-------------------------------------------------------

#------- targets ---------------------------------------
all :	lib $(execs)
lib :
	cd $(SRC_DIR)/; make
clean :
	cd $(SRC_DIR)/; make clean;
	@echo "Deleting library objects, executables and associated objects."
	@/bin/rm -f $(theutils) $(execs)
	@/bin/rm -f *~
	@echo "Done"
deepclean :
	make clean
	@echo "Deleting Toolkit library"
	@/bin/rm -f $(thelib)

#-------------------------------------------------------

#-------- specific rules -------------------------------
exposure.exe : exposure.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
#-------------------------------------------------------
