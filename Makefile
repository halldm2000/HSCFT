#*******************************************************************************
# Makefile for HSCFT Simulation
#*******************************************************************************
SOURCES = \
	CODE/FFTWManager.c\
	CODE/FileIO.c\
	CODE/FileManager.c\
	CODE/HydroModel.c\
	CODE/MathManager.c\
	CODE/MemoryManager.c\
	CODE/SimulationManager.c\
	CODE/StressModel.c\
	CODE/ThermoModel.c\
	CODE/WallManager.c

#PARAMETER DEFINITIONS
EXE =		hscft
MPI_LIBS =      -lrfftw_mpi -lfftw_mpi
LIBS =		-lm -lrfftw -lfftw 
MPICC =		mpicc
CC =		gcc
DEFINES =       -DUSE_MPI

#ARCHITECTURE SPECIFIC COMPILE FLAGS
CFLAGS_G4 = -funroll-loops -O3 -mtune=G4
CFLAGS_G5 = -funroll-loops -O3 -mtune=G5 -L$(HOME)/lib -I$(HOME)/include
CFLAGS_ATHLON = -march=athlon  -funroll-loops -O3 -ffast-math -L/usr/local/lib
CFLAGS_ATHLON64 = -march=athlon64  -funroll-loops -O3 -ffast-math -L/usr/local/lib -L/usr/local/fftw/2.1.5/gcc/lib/ -I/usr/local/fftw/2.1.5/gcc/include/
CFLAGS_PENTIUM4 = -march=pentium4 -mcpu=pentium4 -O3 -funroll-loops -I/users/halldm/fftw-2.1.5/include -L/users/halldm/fftw-2.1.5/lib -I$(MPI_ROOT)/include -L$(MPI_ROOT)/lib
CFLAGS_XRAID2 = -funroll-loops -O3 -mtune=G5 -I/Volumes/XRaid2/daily/halldm/INCLUDES -L/Volumes/XRaid2/daily/halldm/LIBS

#BUILD TARGETS FOR VARIOUS ARCHITECTURES WITH AND WITHOUT MPI
athlon:
	$(CC) $(CFLAGS_ATHLON) $(SOURCES) $(LIBS) -o hscft_athlon
athlon64:
	$(CC) $(CFLAGS_ATHLON64) $(SOURCES) $(LIBS) -o hscft_athlon64
athlonmpi:
	$(MPICC) $(CFLAGS_ATHLON) $(SOURCES) $(LIBS) $(MPI_LIBS) $(DEFINES) -o hscft_athlonmpi
g4:
	$(CC) $(CFLAGS_G4) $(SOURCES) $(LIBS) -o hscft_g4
g4mpi:
	$(MPICC) $(CFLAGS_G4) $(SOURCES) $(LIBS) $(MPI_LIBS) $(DEFINES) -o hscft_g4mpi
g5:
	$(CC) $(CFLAGS_G5) $(SOURCES) $(LIBS) -o hscft_g5
pentium:
	$(CC) $(CFLAGS_PENTIUM4) $(SOURCES) $(LIBS) -o hscft_pentium
pentiummpi:
	$(MPICC) $(CFLAGS_PENTIUM4) $(SOURCES) $(MPI_LIBS) $(LIBS) $(DEFINES) -o hscft_pentiummpi
xraid2:
	$(CC) $(CFLAGS_XRAID2) $(SOURCES) $(LIBS) -o hscft_g5
