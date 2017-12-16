//******************************************************************************
//  FileManager.h
//  HSCFT_Diblock
//  Created by halldm on Fri Jul 23 2004.
//******************************************************************************
#ifndef FileManager_DataWriter_H
#define FileManager_DataWriter_H
#include "Constants.h"

void FileManager_AppendScalarsToFile();
void FileManager_Finalize();
void FileManager_Initialize();
void FileManager_LookUpFloryParameters();
void FileManager_ReadOutputPrefs();
void FileManager_ReadExperimentFile();
void FileManager_ReadCompositionFile(char* fileName);
void FileManager_ReadCopolymerFile(int c);
void FileManager_ReadRestartFiles();
void FileManager_ReadRunParameters();
void FileManager_ReadMaterialFiles();
void FileManager_ReadMotionDescriptor(int c);
void FileManager_WriteDataToFile();
void FileManager_WriteDimensionalParameters();
void FileManager_WriteDomainFile();
void FileManager_WriteFieldsToFile();
void FileManager_WriteRestartFiles();

// Output preference flags
int antialias;          // Use antialiasing to manage Gibb's phenomena
int appendEVF;          // Append volume fraction data
int appendEnergy;       // Append statistics on free energy
int appendForces;       // Append vectors describing various forces
int appendMu;           // Append scalar measure of chemical potential fields
int appendOmegas;       // Append scalar measure of conjugate fields
int appendPhi;          // Append scalar measure of volume fraction fields
int appendPhiGrads;     // Append scalar measure of concentration gradients
int appendResiduals;    // Append measure of residual error
int appendVels;         // Append measures of velocity fields
int writeMu;            // Write chemical potential fields
int writeP;             // Write pressure field
int writePhi;           // Write volume fraction field
int writeOmegas;        // Write conjugate potentials
int writeProps;         // Write propagators
int writeV;             // Write mean velocity field
int writeVDiv;          // Write scalar velocity divergence field
int writeVM;            // Write velocity field for each component
int writeVSubsteps;     // Write velocity steps in fixed point iteration
int writeVTube;         // Write rheometric mean velocity field
int writeW;             // Write relative velocity fields for each component
int writeForces;        // Write various forces fields

// Run-time execution flags
int updateStresses;     // Employ viscoelastic stress-strain model in this run
int updateThermo;       // Compute chemical potential fields in this run
int useUpperConvected;  // Use upper convected time derivs in stress-strain model

int copolymersInComp;   // Number of distinct copolymers in the system
int outputStep;         // Counter for current output step
int lastOutputStep;     // Last step written to file
int lastRestartStep;    // Last step written for use in run restart
int resampleInput;      // input needs to be resampled from different size i.c.s
int traceLevel;         // Highest trace level to display

char description[MAX_STRING];       // String describing current experiment
char experimentFile[MAX_STRING];    // Path to experiment descriptor
char restartDir[MAX_STRING];        // Path to use if re-starting an old run
char compositionFile[MAX_STRING];   // Path to file describing fluid composition
char floryParamsFile[MAX_STRING];   // Path to file describing flory interaction params
char geometryFile[MAX_STRING];      // path to geometry specification file
char outputPrefsFile[MAX_STRING];   // path to preferences file
char runParamsFile[MAX_STRING];     // path to run parameters
char systemFile[MAX_STRING];        // path to system descriptor
char motionDescriptorFile[MAX_STRING];// path to motion descriptor file

char copolymerPath[MAX_C][MAX_STRING];// Path to molecule descriptor

char materialName[MAX_C][MAX_M][MAX_STRING]; // name of each material
char materialPath[MAX_C][MAX_M][MAX_STRING]; // path of each material

#endif
