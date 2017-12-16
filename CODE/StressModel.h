//******************************************************************************
//  StressModel_Maxwell.h
//  HSCFT_Diblock
//  Created by halldm on Fri Jul 23 2004.
//******************************************************************************
#ifndef StressModel_Maxwell_H
#define StressModel_Maxwell_H
#include "Constants.h"

void StressModel_AdvanceSigma();
void StressModel_ComputeElasticForce();
void StressModel_ComputeStressDivisionParams();
void StressModel_ComputeTubeVelocities();
void StressModel_ComputeUpperConvectedTerms();
void StressModel_Initialize();
void StressModel_Finalize();
void StressModel_MeasureSystemProperties();

double K0[MAX_C][MAX_M];            // Bulk modulus of each component
double G0[MAX_C][MAX_M];            // Shear modulus of each component
double meanVTube;                   // Average tube velocity
double rmsAlpha[MAX_C][MAX_M];      // stress division parameter diagnostic
double rmsVTube;                    // Average tube velocity
double tauRep[MAX_C];               // Reptation time for each copolymer
double tauNe[MAX_C];                // Rouse time per entanglement length

SF divVTube;                        // divergence of the tube velocity

VF vTemp1, vTemp2, vTemp3;          // Temporary vector fields
VF vTube;                           // network tube velocity
VF elasticForceSum;                 // elastic force due to all network deformations

TF alpha;                           // monomer stress division parameters
TF G;                               // concentration dep. shear modulus
TF K;                               // concentration dep. bulk modulus
TF zeta;                            // Friction fields for liquid monomers
TF gradVTube;                       // rate of deformation tensor for vtube
TF kappa;                           // rate of strain tensor in shear stress eqn

TF3 elasticForce;                   // elastic forces due to volumetric changes

TF4 sigma;                          // elastic stress tensor for each component
TF4 sigmaOld;                       // elastic stress tensor for each component
TF4 ucTerms;                        // upper convected terms

#endif
