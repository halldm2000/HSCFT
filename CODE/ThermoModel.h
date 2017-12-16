//******************************************************************************
//  ThermoModel_SCFTDiblock.h
//  HSCFT_Diblock
//  Created by halldm on Fri Jul 23 2004.
//******************************************************************************
#ifndef ThermoModel_SCFTDiblock_H
#define ThermoModel_SCFTDiblock_H
#include "Constants.h"

void ThermoModel_ApplyExponentialLaplcian(SF field,SF result,double step,int c);
void ThermoModel_ComputeChemicalPotentials();
void ThermoModel_ComputeAuxilliaryField(int c);
void ThermoModel_ComputePartitionFcn(int c);
void ThermoModel_ComputePropagator(int c);
void ThermoModel_ComputeResidualErrors();
void ThermoModel_DisplayStatus();
void ThermoModel_Finalize();
void ThermoModel_Initialize();
void ThermoModel_MeasureSystemProperties();
void ThermoModel_ReduceResidualErrors();
void ThermoModel_SeekLocalEquilibrium();

int materialType[MAX_C];                // Material identifier: solid, liquid, polymer
int maxSCFTSteps;                       // Max SCFT steps before bailing out
int stepsPerBlock;                      // Number of steps per block in q calc
int ns[MAX_C];                          // Number gridpoints along polymer index
int scftStep;                           // Counts number of scftSteps needed to converge
int scftInterval;                       // Number of step between status displays

double blockLength[MAX_C][MAX_B];       // Degree of polymerization of each diblock
double chi[MAX_C][MAX_M][MAX_C][MAX_M]; // Flory Chi parameter table
double CahnNumber;                      // Dimensionless fluid tension
double deltaS[MAX_C][MAX_B];            // Size of steps along s in q calculation
double epsilon;                         // steepest descent initial step size
double F;                               // Total free energy for the system
double lineDensity[MAX_C][MAX_M][MAX_B];// Linear density of material on each copolymer
double maxSCFTResidual;                 // Maximum acceptable SCFT residual error
double maxThermoStep;                   // Upper limit on stepsize in error reduction
double meanN;                           // Characteristic polymerization index
double meanOmega[MAX_C][MAX_M];         // Average value of conjugate potential fields
double minN;                            // Shortest polmyerization index (for polymers)
double N[MAX_C];                        // Total degree of polymerization of all blocks
double Ne[MAX_C];                       // Entanglement length of a given copolymer
double rmsOmegaErr[MAX_C][MAX_M];       // Residual error in each conjugate field

double Q[MAX_C];                // Partition function for a single copolymer
double scftErr;                 // Max residual error in pseudo-equilbrium conditions
double scftErr1;                // SCFT Error sampled after first iteration step
double rhoTemp[MAX_S];          // Temporary variable used in quadrature calculation
double rmsStepLength;           // step length taken in line minimization
double rmsMu[MAX_C][MAX_M];     // RMS value of chemical potential fields
double rmsOmega[MAX_C][MAX_M];  // RMS value of conjugate potential fields
double rmsPhi[MAX_C][MAX_M];    // RMS value of target monomer concentrations
double rmsPhiAux[MAX_C][MAX_M]; // RMS value of axilliary monomer concentrations

TF mu;                          // chemical potential fields
TF phiAux;                      // auxiliary monomer concentrations
TF q;                           // propagator, starting at s=0 end
TF qdagger;                     // copropagator, starting at s=1 end
TF omega;                       // chemical field conjugate to monomer concentrations
TF omegaErr;                    // residual errors in conjuagte fields
TF omegaErrOld;                 // residual errors used in parabolic minimization

#endif
