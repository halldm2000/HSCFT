//******************************************************************************
//  WallManager.h
//******************************************************************************
#ifndef FileManager_WallManager_H
#define FileManager_WallManager_H
#include "Constants.h"

void WallManager_AllocateMemory();
void WallManager_AllocateMemoryForParticle(int c,int numP);
void WallManager_FreeMemoryForParticle(int c);
void WallManager_Finalize();

void WallManager_ComputePivotPoint(int c);
void WallManager_ComputeExternalForce();
void WallManager_ComputeLiquidFraction();
void WallManager_ComputeWallVelocityFields();
void WallManager_ComputeParticleBounds(int c);
void WallManager_ConstructParticleField(int c);
void WallManager_ConstructParticleMoment(int c);
void WallManager_ConstructParticlePESurface(int c);
void WallManager_ConstructSingleParticlePotential();
void WallManager_ConstructSystemFromFile();
void WallManager_GetParticleField(int c, int i, SF field);
void WallManager_GetParticlePEField(int c, int i, SF field);
void WallManager_ImportBitmap();
void WallManager_Initialize();
void WallManager_FillGeometry(char* geometryFile,char* compositionFile);
void WallManager_FillBox(SF field,double val, double* center,double* width);
void WallManager_FillCylinder(SF field,double val, double* center, double radius);
void WallManager_FillSphere(SF field,double val, double* center, double radius);
void WallManager_FillTriangle(SF field,double val, double* center,double* width);
void WallManager_FitDataInMask(TF field,int numComponents, int total);
void WallManager_MakeGaussianFilter();
void WallManager_MakeGeometryFromFile(SF field,char* fileName);
void WallManager_MeasureForcesAndTorques(int c, int j);
void WallManager_SetObjectLocations(int c);
void WallManager_SmoothEdges(SF wallField);

int numLiquids;             // Number of liquid fields in simulation
int numSolids;              // Number of solid fields in the simulation
int numConstrained;         // Number of externally driven materials
int numParticles[MAX_C];    // Number of solid particles of type c

double efMagnitude;         // Magnitude of externally applied force
double filterWidth;         // Radius of Gaussian smoothing filter
double fLiquid;             // Total volume fraction of polymeric materials
double fSolid;              // Total volume fraction of solid materials
double radius;              // Radius of circular shape

double fM[MAX_C][MAX_M];    // Volume Fraction of each monomer on a copolymer
double tConstrained[MAX_C]; // Translation constraint flag
double aConstrained[MAX_C]; // Angular constraint flag
double particleMass[MAX_C]; // Total mass of a particle of type c
double particleMoment[MAX_C];// Angular momement of a particle of type c
double angSpeedMax[MAX_C];  // Maximum angular speed of each driven object
double angFreq[MAX_C];      // Oscillation frequency for rotating objects
double angPhase[MAX_C];     // Initial phase angle for oscillating objecys
double transVMax[MAX_C][D]; // Maximum velocity for translating objects
double transFreq[MAX_C][D]; // Oscillation frequency for translating objects
double transPhase[MAX_C][D];// Initial phase angle for translating objects
double pivot[MAX_C][D];     // Center of mass of each particle template

T2 netTorque[MAX_C];        // Net center of mass torque on each particle
T2 netForce[MAX_C];         // Net force on each particle 
T2 rCM[MAX_C];              // Center of mass position for each particle
T2 vCM[MAX_C];              // Center of mass velocity for each particle
T2 avCM[MAX_C];             // Angular velocity of each particle
T2 aCM[MAX_C];              // Angular displacement of each particle

V3 center;                  // Center position of shape (x,y,z)
V3 width;                   // Dimensions of current shape (wx,wy,wz)
V3 rprime;                  // Position relative to particle center of mass

IBOX rbounds;               // bounding box for rotated particle
IBOX bounds[MAX_C];         // bounding box for each particle
IBOX peBounds[MAX_C];       // bounding box for each particle

SF filter;                  // Gaussian filter used to smooth geometries
SF liquid;                  // Local sum of all liquid volume fractions
SF mask;                    // Mask for construction of geometric shapes
SF netParticlePE;           // Total potential energy due to particle-particle interactions
SF potential;               // Particle-paticle potential energy function
SF solid;                   // Local sum of solid volume fractions 
SF total;                   // Local sum of all volume fractions 
SF temp1;                   // temporary storage for intermediate calculations
SF temp2;                   // temporary storage for intermediate calculations
SF temp3;                   // temporary storage for intermediate calculations
VF particle;                // vector of particle shapes
VF particlePE;              // vector of particle potential energy surfaces
VF vWall;                   // Velocity field for all solids
TF phi;                     // monomer volume fraction fields

char bitmapPath[MAX_STRING];// Path to text file describing a bitmapped object

#endif
