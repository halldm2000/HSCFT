//******************************************************************************
//  SimulationManager.h
//  Created:    halldm  Mar 06 2005.
//******************************************************************************
#ifndef SimulationManager_H
#define SimulationManager_H
#include "Constants.h"

void SimulationManager_AdvanceFields();
void SimulationManager_AntialiasFields();
void SimulationManager_ComputeDerivedQuantities();
void SimulationManager_ComputeDimensionlessGroups();
void SimulationManager_ComputeAverageN();
void SimulationManager_ComputeMonomerVolumeFractions(int c);
void SimulationManager_DisplayStatus();
void SimulationManager_DisplayMessage(char* message);
void SimulationManager_DisplayStatus();
void SimulationManager_FinalizeModel();
void SimulationManager_InitializeData();
void SimulationManager_InitializeModel(int argc, char** argv);
void SimulationManager_InitializeSystemSize();
void SimulationManager_MeasureMonomerVolumeFractions();
void SimulationManager_MeasureCopolymerVolumeFractions();
void SimulationManager_MeasureSystemProperties();
void SimulationManager_ReadCommandLineParams(int argc,char** argv);
void SimulationManager_SeedRandomNumberGenerator();
void SimulationManager_Trace(char* message,int level);

int appendInterval;             // Number of steps to skip between data appends
int numBlocks;                  // Total Number of blocks on all copolymers
int numC;                       // Number of distinct copolymer species
int dim;                        // Dimension of this run: 2d or 3d
int displayInterval;            // Number of steps to skip between data displays
int gridSize[D];                // Number of gridpoints in each dimension
int oldGridSize[D];             // Number of gridpoints in initial condition data
int newExperiment;              // New experiement flag
int nx,ny,nz;                   // Convenient names for gridsize dimensions
int nmin;                       // smallest gridsize dimension
int numProcesses;               // Total number of processes spawned by simulation
int programStatus;              // Status flag for use in adaptive timestepping
int processRank;                // Unique identifier for this mpi process
int restart;                    // Indicates continuation of a stopped run
int timeStep;                   // Number of time steps that have elapsed thus far
int restartTimestep;            // Time step index at which the restarted run began
int numB[MAX_C];                // Number of blocks on each copolymer
int numM[MAX_C];                // Number of  monomers on each copolymers

// Dimensional and dimensionless conversion parameters
double Ca;                      // capillary number
double De;                      // Deborah number
double RE;                      // Reynolds number
double lC;                      // characteristic length scale;
double Gamma;                   // Friction factor
double Gc;                      // characteristic stress
double kT;                      // thermal energy
double muC;                     // characteristic chemical potential
double tauD;                    // characteristic diffusion time (smallest tauE)
double tC;                      // characteristic time scale;
double vC;                      // charcateristic convective velocity
double wC;                      // characteristic diffusion velocity
double zeta0C;                  // characteristic monomer friction
double zetaC;                   // characteristic monomer friction

double  deltaT;                 // time step for outer loop
double  deltaT_v;               // time step for v fixed point iteration
double  L[3];                   // Domain size
double  lambda[MAX_C];          // copolymer fraction in whole system
double  lambdaC[MAX_C];         // copolymer fraction in its composition
double  meanV;                  // Average magnitude of the mean velocity field
double  norm;                   // commonly used scaling parameter 1/(nx*ny*nz)
double  runEndTime;             // Time at which the simulation should cease
double  runTime;                // Total elapsed dimensionless time
double  volume;                 // Volume of the box
double  timeEnd;                // Clock time at which the simulation completed
double  timeStart;              // Clock time at which the simulation began
double  writeInterval;          // run time to skip between write events
double  restartInterval;        // run time to skip between restart saves
double  restartTime;            // run time at which run should be restarted

double  Rg[MAX_C];              // Unperturbed radius of gyration of each copolymer
double  Z[MAX_C];               // Number of "entanglements" per copolymer

double  density[MAX_C][MAX_M];  // density of each component
double  f[MAX_C][MAX_M];        // Average volume fraction of each monomer
double  mViscosity[MAX_C][MAX_M];// viscosity of each material
double  meanPhi[MAX_C][MAX_M];  // Average fraction of each monomer type
double  meanPhiAux[MAX_C][MAX_M];// Spatial average of auxillary monomer fractions

#endif
