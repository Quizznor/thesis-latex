##################################################
# makefile for catalog-based analysis            #
# last edit: 2022-04-14                          #
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
CAT_DIR = ../Catalog_Based

CFITSIO_INCDIR = $(CFITSIO_DIR)/include
CFITSIO_LIBDIR = $(CFITSIO_DIR)/lib

CHEALPIX_INCDIR = $(CHEALPIX_DIR)/include
CHEALPIX_LIBDIR = $(CHEALPIX_DIR)/lib

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
	evolution.exe

exeobjs = $(patsubst %.exe,%.o,$(execs))

HEADERS = $(patsubst %.o,%.h)

thelib = $(SRC_DIR)/libtkit.a
theutils = $(CAT_DIR)/utils.a
#-------------------------------------------------------

#-------- rules ----------------------------------------
# rules for the library sources
$(OBJ)%.o:%.cc %.h
	$(COMPILE.cc) $(DBGFLAGS) $(OPTFLAGS) $(WFLAGS) $(INCDIR) -o $@ $<
# rules for the executable sources
$(OBJ)%.o:%.cc
	$(COMPILE.cc) $(DBGFLAGS) $(OPTFLAGS) $(WFLAGS) $(INCDIR) -o $@ $<
#-------------------------------------------------------

#------- targets ---------------------------------------
all :	lib $(theutils) $(execs)
lib :
	cd $(SRC_DIR)/; make
	cd $(CAT_DIR)/; make
clean :
	cd $(SRC_DIR)/; make clean;
	cd $(CAT_DIR)/; make clean;
	@echo "Deleting library objects, executables and associated objects."
	@/bin/rm -f $(execs) $(exeobjs)
	@/bin/rm -f *~
	@echo "Done"
deepclean :
	make clean
	@echo "Deleting Toolkit library"
	@/bin/rm -f $(thelib) $(theutils)

#-------------------------------------------------------

#-------- specific rules -------------------------------
evolution.exe : evolution.o $(theutils) $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
#-------------------------------------------------------
