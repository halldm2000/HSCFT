//******************************************************************************
//  FFTWManager.h
//  HSCFT_Multiblock
//
//  Manages construction of fast-Fourier transform routines using
//  FFTW2.1.5 (Fastest Fourier Tranform in the West)
//******************************************************************************
#ifndef FFTWManager_H
#define FFTWManager_H

#include "Constants.h"

void FFTWManager_CreatePlans();
void FFTWManager_Finalize();
void FFTWManager_Initialize();
void FFTWManager_MakeWaveTable();

RFFTWND_PLAN plan;          // Plan for transforming to Fourier space
RFFTWND_PLAN iplan;         // Plan for transforming to real space

int totalSize;              // number of doubles in scalar fftwmatrix allocation
int nxLocal;                // number of x points assigned to this processor
int nyTrans;                // number of y frequencies assigned to this processor
int startX,startY,startZ;   // leftmost point assigned to this processor
int startFX,startFY,startFZ;// leftmost point assigned to this processor
int maxX,maxY,maxZ;         // maximum gripoint assigned to this processor
int maxFX,maxFY,maxFZ;      // max gridpoint in Fourier space
int strideX,strideY;        // width of data blocks in real space
int strideFX,strideFY;      // width of data blocks in Fourier space

double k2Max;               // maximum wavenumber, square magnitude
double kMax;                // maximum wavenumber magnitude

VEC waveTableX;             // table of wavenumbers
VEC waveTableY;             // table of wavenumbers
VEC waveTableZ;             // table of wavenumbers

SF work1;                   // scalar temporary data storage
SF work2;                   // scalar temporary data storage
SF work3;                   // scalar temporary data storage
VF workVector;

#endif
