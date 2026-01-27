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
SRC_DIR = ../../Utilities

INCDIR = -I$(SRC_DIR) \
	$(shell root-config --cflags)

LIBDIR = -lm \
	-lMinuit \
	$(shell root-config --glibs)
#-------------------------------------------------------

#------- alias -----------------------------------------
execs = \
	propagate.exe 
	
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
propagate.exe : propagate.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
#-------------------------------------------------------
