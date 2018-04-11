#
#   HPCCS Make File
#   Software to Calculate Mobilities - Version 1.0
#   Developed by Leandro Zanotto - Center for Computational Engineering & Sciences#   
#	Unicamp - University of Campinas
#	Use this file to compile using GCC
#
#  For Debug CCFLAGS= -g -O2 -std=c++0x -fopenmp -m64 -mtune=native -ftree-vectorize -mavx -ffast-math -fopt-info-vec-optimized
SRC_DIR=src
CCFLAGS= -O3 -std=c++0x -fopenmp -m64 -mtune=native -mavx -mfma -ffast-math
INC_DIR = ${HOME}/HPCCS
CC=g++
LDFLAGS=
OBJ= molecule.o atomMLJ.o fhandler.o globals.o hpccs.o mobil2.o potentialHe.o potentialN2.o rotate.o diffeq_deriv.o gsang.o mt.o

all: hpccs

hpccs: ${OBJ}
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -o $@ ${OBJ}

mt.o: $(SRC_DIR)/mt.cc
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

fhandler.o: $(SRC_DIR)/fhandler.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

molecule.o: $(SRC_DIR)/molecule.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

hpccs.o: $(SRC_DIR)/hpccs.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

atomMLJ.o: $(SRC_DIR)/atomMLJ.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

globals.o: $(SRC_DIR)/globals.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

mobil2.o: $(SRC_DIR)/mobil2.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

potentialHe.o: $(SRC_DIR)/potentialHe.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

potentialN2.o: $(SRC_DIR)/potentialN2.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

rotate.o: $(SRC_DIR)/rotate.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

gsang.o: $(SRC_DIR)/gsang.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

diffeq_deriv.o: $(SRC_DIR)/diffeq_deriv.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

clean:
	rm *.o hpccs output.out
