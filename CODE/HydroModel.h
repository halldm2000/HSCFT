//******************************************************************************
//  HydroModel_MultiFluid.h
//  HSCFT_Diblock
//  Created by halldm on Fri Jul 23 2004.
//******************************************************************************
#ifndef HydroModel_MultiFluid_H
#define HydroModel_MultiFluid_H
#include "Constants.h"

void HydroModel_AdvancePhi();
void HydroModel_AdvanceV();
void HydroModel_AdvanceW();
void HydroModel_AdvanceParticleVelocities();
void HydroModel_ApplyProjectionMethod();
void HydroModel_ComputeExcessVolumeFractions();
void HydroModel_ComputeStressDivisionParams();
void HydroModel_ComputeTotalOsmoticPressure();
void HydroModel_ComputeViscousForce();
void HydroModel_ComputeWallForces();
void HydroModel_Finalize();
void HydroModel_Initialize();
void HydroModel_MeasureSystemProperties();
void HydroModel_RescalePolymerFields();
void HydroModel_SeekSteadyState();

int vDisplayInterval;           // number of steps between display in v iteration
int vCount;                     // number of steps in current v iteration
int vCountTotal;                // total number of v steps takenin run
int suppressVDrift;             // Flag indicating that mean v should be zeroed

double avf[MAX_C][MAX_M];       // apparant volume fraction for each component
double evf[MAX_C][MAX_M];       // excess volume fraction for each component
double maxGradPhi[MAX_C][MAX_M];// maximum volume fraction gradients
double maxPhi[MAX_C][MAX_M];    // maximum volume fractions
double maxElasticForce;         // measure of maximum elastic force
double maxOsmoticForce;         // measure of maximum osmotic pressure
double maxWallForce;            // maximum force of solid on fluid
double maxViscousForce;         // maximum viscous force
double maxV;                    // maximum mean velocity
double minPhi[MAX_C][MAX_M];    // minimum value of the concentraiton field
double meanPhi[MAX_C][MAX_M];   // average value of the concentraiton field
double meanZeta0[MAX_C];        // average friction coefficient for each copolymer
double meanP;                   // average pressure
double meanV;                   // average mean velocity magnitude
double monomerFriction[MAX_C][MAX_M];// dimensional friction coefficients
double noise;                   // amplitude of gaussian random noise
double rmsdVdL;                 // rms rate of change in v relaxation 
double rmsGradPhi[MAX_C][MAX_M];// rns value of phi gradient
double rmsV;                    // rms value of mean velocity field
double rmsVDiv;                 // rms residual error in pressure
double rmsVErr;                 // rms residual error in v field
double rmsVM[MAX_C][MAX_M];     // rms values of monomer velocity fields
double rmsW[MAX_C][MAX_M];      // rms values of relative velocity fields
double vThreshold;              // target residual error in force balance eqn
double zeta0[MAX_C][MAX_M];     // friction coefficient of each monomer

SF divVStar;                    // divergence of vstar (compressible v field)
SF P;                           // hydrodynamic pressure field
SF vDiv;                        // divergence of mean velocity field
SF vMag;                        // magnitude of mean velocity field
SF viscosity;                   // Effective viscosity at each point in the fluid

VF brownianForce;               // force due to random thermal fluctuations
VF externalForce;               // externally applied force field to induce flow
VF dVdL;                        // mean velocity field
VF flux;                        // local flux of a single monomer species
VF gradP;                       // hydrodynamic pressure gradient
VF gradPi;                      // total osmotic pressure gradient
VF osmoticForce;                // osmotic force due to a single species
VF v;                           // mean velocity field
VF vOld;                        // mean velocity field
VF vErr;                        // residual error remaining in force balance eqn
VF vForces;                     // forces acting on v field
VF viscousForce;                // osmotic force due to a single species
VF vLiquid;                     // mean velocity field
VF vStar;                       // velocity field used in projection method
VF wForces;                     // forces acting on v field
VF wallForces;                  // total force applied by walls
VF convectiveTerm;              // convective transport term for RE>0

TF3 fPsi;                        // solid object force fields
TF3 gradMu;                      // gradient of chemical potential fields
TF3 gradPhi;                     // gradient of conecntration fields
TF3 vPhi;                        // monomer velocity fields
TF3 w;                           // relative velocity fields

#endif

