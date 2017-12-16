//******************************************************************************
//  FFTWManager.c
//******************************************************************************
#include "Constants.h"
#include "FileIO.h"
#include "FileManager.h"
#include "FFTWManager.h"
#include "MathManager.h"
#include "MemoryManager.h"
#include "SimulationManager.h"

//******************************************************************************
// FFTWManager_Initialize
//******************************************************************************
void FFTWManager_Initialize()
{
    SimulationManager_Trace("FFTWManager_Initialize",2);

    FFTWManager_CreatePlans();
    waveTableX= calloc(maxFX,sizeof(double));
    waveTableY= calloc(maxFY,sizeof(double));
    waveTableZ= calloc(maxFZ,sizeof(double));
    SimulationManager_DisplayMessage("Creating Frequency Lookup Tables.");
}

//******************************************************************************
// FFTW_CreatePlans
//******************************************************************************
void FFTWManager_CreatePlans()
{
    SimulationManager_DisplayMessage("Creating FFTW Plans.");
    
    startX = startY = startZ =0;
    startFX= startFY= startFZ=0;
    plan = CREATE_FORWARD_TRANSFORM(dim);
    iplan= CREATE_INVERSE_TRANSFORM(dim);

    if (dim==3)
    {
        GET_DATA_SIZES_3D(plan,nxLocal,startX,nyTrans,startFY,totalSize);
        GET_REAL_LIMITS_3D(maxX,maxY,maxZ);
        GET_FOURIER_LIMITS_3D(maxFX,maxFY,maxFZ);
        GET_REAL_STRIDE_3D(strideX,strideY);
        GET_FOURIER_STRIDE_3D(strideFX,strideFY);
    }
    if (dim==2)
    {      
        GET_DATA_SIZES_2D(plan,nxLocal,startX,nyTrans,startFY,totalSize);
        GET_REAL_LIMITS_2D(maxX,maxY,maxZ);
        GET_FOURIER_LIMITS_2D(maxFX,maxFY,maxFZ);
        GET_REAL_STRIDE_2D(strideX,strideY);
        GET_FOURIER_STRIDE_2D(strideFX,strideFY);
    }
    if (nxLocal <4) 
    {
        SimulationManager_DisplayMessage("Too many nodes. nxLocal cannot be less than 4.");
        exit(1);
    }
    
    // Display the which blocks of data are allocated to which processor
    printf("%d: nxLocal=%d  startX=%d  nyTrans=%d  startFY=%d  totalSize=%d\n",
           processRank, nxLocal, startX, nyTrans,startFY,totalSize);
    printf("maxX=%d mayY=%d maxZ=%d\n"      ,maxX,maxY,maxZ);
    printf("maxFX=%d mayFY=%d maxFZ=%d\n"   ,maxFX,maxFY,maxFZ);
    printf("strideX=%d strideY=%d\n"        ,strideX,strideY);
    printf("strideFX=%d strideFY=%d\n"      ,strideFX,strideFY);

    work1 = MemoryManager_AllocateScalarField();
    work2 = MemoryManager_AllocateScalarField();
    work3 = MemoryManager_AllocateScalarField();
    workVector = MemoryManager_AllocateVectorField(D);

    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FFTW_DestroyPlans
//******************************************************************************
void FFTWManager_Finalize()
{
    DESTROY_TRANSFORM(plan);
    DESTROY_TRANSFORM(iplan);
    
    MemoryManager_FreeScalarField(work1);
    MemoryManager_FreeScalarField(work2);
    MemoryManager_FreeScalarField(work3);
    MemoryManager_FreeVectorField(workVector,D);
    
    free(waveTableX);
    free(waveTableY);
    free(waveTableZ);
}

//******************************************************************************
// FFTWManager_MakeWaveTable
// Compute the wavenumbers for each frequency
//******************************************************************************
void FFTWManager_MakeWaveTable()
{
    int fx,fy,fz;
    FORFX(fx)
    {
        int localPos    = fx + startFX;
        int freq        = (localPos <= nx/2) ? localPos : localPos-nx;
        double waveNum  = TWO_PI*freq/L[0];
        waveTableX[fx]= waveNum;
        //printf("%d wavetableX[%d]=%f\n",processRank,fx,waveTableX[fx]);
    }
    
    FORFY(fy)
    {
        int localPos    = fy + startFY;
        int freq        = (localPos <= ny/2) ? localPos : localPos-ny;
        double waveNum  = TWO_PI*freq/L[1];
        waveTableY[fy]= waveNum;
    }
    
    FORFZ(fz)
    {
        int localPos    = fz + startFZ;
        int freq        = (localPos <= nz/2) ? localPos : localPos-nz;
        double waveNum  = TWO_PI*freq/L[2];
        waveTableZ[fz]= waveNum;
    }
    
    k2Max = SQR(TWO_PI*nx/L[0]/2.0) + SQR(TWO_PI*ny/L[1]/2) + SQR(TWO_PI*nz/L[2]/2);
    kMax = sqrt(k2Max);    
    WAIT_FOR_ALL_PROCESSES;
}
