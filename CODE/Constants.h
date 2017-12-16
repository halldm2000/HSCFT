#ifndef Constants_H
#define Constants_H

#include <math.h>
#include <rfftw.h>
#include <fftw.h>
#include <stdio.h>
#include <sys/stat.h>
#include <time.h>
#include <string.h>

// General use constants
#define PI                  3.141592653589793238
#define TWO_PI              (2.0*PI)
#define MAX_STRING          256
#define GOLDEN_SECTION      0.61803399
#define GOLDEN_RATIO        (1.0+GOLDEN_SECTION)

// Physical constants and conversion factors
#define NA                  6.022e+23       // molecules/mol
#define Kb                  1.3807e-16      // erg/K
#define M0                  100.0           // g/mol
#define MONOMER_SIZE        1.8e-7          // cm
#define DENSITY             1.0             // g/cm3
#define NDENSITY            (NA*DENSITY/M0) // molecules/cm3
#define R                   1.987           // cal/mol K in calories,not ergs
#define T                   473

// Mathematical Macros
#define MAX(a,b)            ( (a)>(b)?(a):(b)   )
#define MIN(a,b)            ( (a)<(b) ? (a):(b) )
#define MIN4(a,b,c,d)       ( MIN(MIN(a,b), MIN(c,d))   )
#define MAX4(a,b,c,d)       ( MAX(MAX(a,b), MAX(c,d))   )
#define SQR(a)              ( (a)*(a)           )
#define SQUARE_MAG(v)       ( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] )
#define CUBE(a)             ( (a)*(a)*(a)       )
#define DELTA_FCN(i,j)      ( (i)==(j)          )
#define THETA(a,b)          ( (a)>(b) ? 1:0     )
#define NONZERO(denom)      ( denom + (denom==0))
#define ROUND(value)        ( (int)(value+0.5)  )
#define IMOD(a,na)          ( ((a%na)+na)%na    )

// Model-specific macros and constants (for static ram allocation)s
#define D               3   // Number of spatial dimensions
#define MAX_B           10  // Maximum number of blocks per copolymer
#define MAX_C           10  // Maximum number of distinct copolymers in the system
#define MAX_M           10  // Maximum number of distinct monomers per copolymer
#define MAX_S           500 // Max number of steps in the s direction
//#define MAX_P           20  // Max number of particles of a given type

#define DEFAULT_TRACE_LEVEL 5

// Switch out blocks of code so we can run with or without MPI installed
#ifdef USE_MPI  
    // The need for MPI forces the use of FFTW2.1.5 vs FFTW3
    #include <mpi.h>
    #include <fftw_mpi.h>
    #include <rfftw_mpi.h>
    #define BROADCAST_STRING(str)               MPI_Bcast(str, MAX_STRING, MPI_CHAR,0, MPI_COMM_WORLD)
    #define BROADCAST_DOUBLES(values,num)       MPI_Bcast(values, num, MPI_DOUBLE,0, MPI_COMM_WORLD)
    #define BROADCAST_INTEGERS(values,num)      MPI_Bcast(result, num, MPI_INT , 0, MPI_COMM_WORLD)
    #define CREATE_FORWARD_TRANSFORM(d)         rfftwnd_mpi_create_plan(MPI_COMM_WORLD, d, gridSize, FFTW_REAL_TO_COMPLEX, FFTW_MEASURE)
    #define CREATE_INVERSE_TRANSFORM(d)         rfftwnd_mpi_create_plan(MPI_COMM_WORLD, d, gridSize, FFTW_COMPLEX_TO_REAL, FFTW_MEASURE)
    #define DESTROY_TRANSFORM(tplan)            rfftwnd_mpi_destroy_plan(tplan)
    #define FOURIER_TRANSFORM(field,work3d)    rfftwnd_mpi(plan,1, field[0][0],work3d[0][0],FFTW_TRANSPOSED_ORDER)
    #define INVERSE_TRANSFORM(field,work3d)    rfftwnd_mpi(iplan,1,field[0][0],work3d[0][0],FFTW_TRANSPOSED_ORDER)
    #define FIND_GLOBAL_MAX(localMax,globalMax) MPI_Allreduce(&localMax,&globalMax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD)
    #define FIND_GLOBAL_MIN(localMin,globalMin) MPI_Allreduce(&localMin,&globalMin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD)
    #define FIND_GLOBAL_SUM(localSum,globalSum) MPI_Allreduce(&localSum,&globalSum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD)
    #define FINALIZE_ALL_PROCESSES              MPI_Finalize()
    #define FINDEX(f) ((f[1]*strideFX + f[0])*strideFY + f[2]) // fx,fy are transposed
   
    #define GET_DATA_SIZES_2D(pl,nxL,stX,nyT,stFY,tot) rfftwnd_mpi_local_sizes(pl, &nxL, &stX, &nyT, &stFY, &tot);
    #define GET_DATA_SIZES_3D(pl,nxL,stX,nyT,stFY,tot) rfftwnd_mpi_local_sizes(pl, &nxL, &stX, &nyT, &stFY, &tot);
    #define GET_REAL_LIMITS_3D(maxX,maxY,maxZ) maxX=nxLocal; maxY=ny; maxZ=nz;
    #define GET_REAL_LIMITS_2D(maxX,maxY,maxZ) maxX=nxLocal; maxY=ny; maxZ=1;
    #define GET_FOURIER_LIMITS_3D(maxFX,maxFY,maxFZ) maxFX=nx; maxFY=nyTrans; maxFZ=nz/2+1;
    #define GET_FOURIER_LIMITS_2D(maxFX,maxFY,maxFZ) maxFX=nx; maxFY=nyTrans; maxFZ=1;
    #define GET_REAL_STRIDE_3D(strideX,strideY) strideX=ny; strideY=2*(nz/2+1);
    #define GET_REAL_STRIDE_2D(strideX,strideY) strideX=2*(ny/2+1); strideY=1;
    #define GET_FOURIER_STRIDE_3D(strideFX,strideFY) strideFX=nx; strideFY=nz/2+1;
    #define GET_FOURIER_STRIDE_2D(strideFX,strideFY) strideFX=nx; strideFY=1;

    #define GET_PROCESS_RANK(rank)              MPI_Comm_rank(MPI_COMM_WORLD, &processRank)
    #define GET_NUM_PROCESSES(num)              MPI_Comm_size(MPI_COMM_WORLD, &num)
    #define GET_TIME(time)                      time=MPI_Wtime();
    #define INITIALIZE_ALL_PROCESSES(argc,argv) MPI_Init(&argc,&argv)
    #define RFFTWND_PLAN                        rfftwnd_mpi_plan
    #define WAIT_FOR_ALL_PROCESSES              MPI_Barrier(MPI_COMM_WORLD); 
#else        
    // For simplicity we employ FFTW2.1.5 for the non-MPI case as well
    #define BROADCAST_STRING(str)               // do nothing
    #define BROADCAST_DOUBLES(values,num)       // do nothing
    #define BROADCAST_INTEGERS(values,num)      // do nothing
    #define CREATE_FORWARD_TRANSFORM(d)         rfftwnd_create_plan(d,gridSize,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE | FFTW_IN_PLACE)
    #define CREATE_INVERSE_TRANSFORM(d)         rfftwnd_create_plan(d,gridSize,FFTW_COMPLEX_TO_REAL,FFTW_MEASURE | FFTW_IN_PLACE)
    #define DESTROY_TRANSFORM(tplan)            rfftwnd_destroy_plan(tplan)
    #define FOURIER_TRANSFORM(field,work3d)     rfftwnd_one_real_to_complex(plan , field[0][0], (fftw_complex*)field[0][0]);
    #define INVERSE_TRANSFORM(field,work3d)     rfftwnd_one_complex_to_real(iplan, (fftw_complex*)field[0][0], field[0][0]);
    #define FIND_GLOBAL_MAX(localMax,globalMax) globalMax=localMax;
    #define FIND_GLOBAL_MIN(localMin,globalMin) globalMin=localMin;
    #define FIND_GLOBAL_SUM(localSum,globalSum) globalSum=localSum;
    #define FINALIZE_ALL_PROCESSES              // do nothing

    #define FINDEX(f)   ( (f[0]*strideFX + f[1])*strideFY + f[2] )
    #define GET_DATA_SIZES_3D(pl,nxL,stX,nyT,stFY,tot) nxL=nx; stX=0; nyTrans=ny; stFY=0; totalSize=2*nx*ny*(nz/2+1);
    #define GET_DATA_SIZES_2D(pl,nxL,stX,nyT,stFY,tot) nxL=nx; stX=0; nyTrans=(ny/2+1); stFY=0; totalSize=2*nx*(ny/2+1);
    #define GET_REAL_LIMITS_3D(maxX,maxY,maxZ)          maxX=nx; maxY=ny; maxZ=nz;
    #define GET_REAL_LIMITS_2D(maxX,maxY,maxZ)          maxX=nx; maxY=ny; maxZ=1;
    #define GET_FOURIER_LIMITS_3D(maxFX,maxFY,maxFZ)    maxFX=nx; maxFY=ny; maxFZ=(nz/2+1);
    #define GET_FOURIER_LIMITS_2D(maxFX,maxFY,maxFZ)    maxFX=nx; maxFY=(ny/2+1); maxFZ=1;
    #define GET_REAL_STRIDE_3D(strideX,strideY)         strideX=ny; strideY=2*(nz/2+1);
    #define GET_REAL_STRIDE_2D(strideX,strideY)         strideX=2*(ny/2+1); strideY=1;
    #define GET_FOURIER_STRIDE_3D(strideFX,strideFY)    strideFX=ny; strideFY=(nz/2+1);
    #define GET_FOURIER_STRIDE_2D(strideFX,strideFY)    strideFX=(ny/2+1); strideFY=1;

    #define GET_PROCESS_RANK(rank)              rank=0;
    #define GET_NUM_PROCESSES(num)              num=1;
    #define GET_TIME(time) time=clock()/CLOCKS_PER_SEC
    #define INITIALIZE_ALL_PROCESSES(argc,argv) // do nothing
    #define RFFTWND_PLAN                        rfftwnd_plan
    #define WAIT_FOR_ALL_PROCESSES              // do nothing
#endif //USE_MPI

#define GET_RUN_DIMENSION(dim) if (nz>1) dim=3; else dim=2;

enum runStates  {RUN,STOP};
enum MATERIAL_TYPES {POLYMER,SOLID,SIMPLE_LIQUID};

// Define some dynamic storage types to emphasize their meaning
typedef int           GP[3];        //  Vector of three integers
typedef double        V3[3];        //  Vector of three doubles
typedef double*       VEC;          //  General Vector (1st order tensor)
typedef double**      T2;           //  2nd order tensor
typedef double***     T3;           //  3nd order tensor
typedef double***     SF;           //  Scalar Field
typedef double****    VF;           //  Vector Field
typedef double*****   TF;           //  Tensor Field
typedef double******  TF3;          //  3rd order Tensor Field
typedef double******* TF4;          //  4th order Tensor Field

// STRUCTURES
struct ibox {GP ul,ur,ll,lr; };
struct box  {V3 ul,ur,ll,lr; };

typedef struct ibox IBOX;
typedef struct box BOX;

// Macros for loops in Fourier space
#define FORFX(fx) for(fx=0; fx<maxFX; fx++)
#define FORFY(fy) for(fy=0; fy<maxFY; fy++)
#define FORFZ(fz) for(fz=0; fz<maxFZ; fz++)
#define FOREACHFREQ(f) FORFX(f[0]) FORFY(f[1]) FORFZ(f[2])

#define FORX(x) for (x=0; x<maxX; x++)
#define FORY(y) for (y=0; y<maxY; y++)
#define FORZ(z) for (z=0; z<maxZ; z++)
#define FOREACH(x,y,z) FORX(x) FORY(y) FORZ(z)
#define FOREACH_VINDEX(i,x,y,z) FOREACH_DIM(i) FOREACH(x,y,z)
#define FOREACH_TINDEX(i,j,x,y,z) FOREACH_DIM(i) FOREACH_VINDEX(j,x,y,z)
#define FOREACH_GRIDPOINT(r) FORX(r[0]) FORY(r[1]) FORZ(r[2])

// Macros for other common loops
#define FOREACH_S(s)            for(s=0; s<ns; s++)
#define FOREACH_BLOCK(b,c)      for(b=0; b<numB[c]; b++)
#define FOREACH_MATERIAL(c)     for(c=0; c<numC; c++)
#define FOREACH_PARTICLE(i,c)   for(i=0; i<numParticles[c]; i++)
#define FOREACH_LIQUID(c)       for(c=0; c<numLiquids; c++)
#define FOREACH_POLYMER(c)      FOREACH_LIQUID(c) if(materialType[c]==POLYMER)
#define FOREACH_SIMPLE_LIQUID(c) FOREACH_LIQUID(c) if(materialType[c]==SIMPLE_LIQUID)
#define FOREACH_SOLID(c)        for(c=numLiquids; c<numC; c++)
#define FOREACH_DIM(d)          for(d=0; d<D; d++)
#define FOREACH_MONOMER(m,c)    for(m=0; m<numM[c]; m++)
#define FOREACH_LOWER(r,c)      for(r=0; r<D; r++) for (c=r; c<D; c++)
#define FOREACH_UPPER(r,c)      for(r=0; r<D; r++) for (c=0; c<=r; c++)
#define FORALL_LIQUID_COMPONENTS(m,c) FOREACH_LIQUID(c) FOREACH_MONOMER(m,c)
#define FORALL_POLYMER_COMPONENTS(m,c) FOREACH_POLYMER(c) FOREACH_MONOMER(m,c)
#define FORALL_SIMPLE_LIQUID_COMPONENTS(m,c) FOREACH_SIMPLE_LIQUID(c) FOREACH_MONOMER(m,c)
#define FORALL_SOLID_COMPONENTS(m,c)  FOREACH_SOLID(c)  FOREACH_MONOMER(m,c)
#define FORALL_MONOMERS(m,c) FOREACH_MATERIAL(c) FOREACH_MONOMER(m,c)
#define FORALL_BLOCKS(b,c)   FOREACH_MATERIAL(c) FOREACH_BLOCK(b,c)

// Macros for Fourier space computations
#define GETWAVENUMS(f,k)  k[0]=waveTableX[f[0]]; k[1]=waveTableY[f[1]]; k[2]=waveTableZ[f[2]]; 
#define INDEX(x,y,z) ((x*strideX + y)*strideY + z)

#endif
