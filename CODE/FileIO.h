//******************************************************************************
//  FileIO.h
//  HSCFT_Diblock
//  Created by halldm on Fri Jul 23 2004.
//
//  Lsow level file input-output routines
//******************************************************************************
#ifndef FileIO_H
#define FileIO_H
#include "Constants.h"

struct stat FileIO_statBuffer;

void FileIO_AppendScalar(char* fieldName,double value);
void FileIO_AppendVector(char* fieldName,double* values,int dim);
void FileIO_AppendParticleData(char* name,T2 val[MAX_C]);
void FileIO_AppendVariableSizeT3(char* field,T3 t3, int td3, int* td2, int td1);
void FileIO_AppendTensor(char* fieldName,double values[MAX_C][MAX_M]);
void FileIO_CopyToRunDir(char* fileDescriptor);
void FileIO_CreateDirectories();
void FileIO_GenerateUniqueRunPath();
void FileIO_GetStartupFile(char* fileName,char* pathToFile);
double FileIO_ReadLastScalar(char* fieldName);
double FileIO_ReadNthScalar(char* fileName,int index);
void FileIO_ReadScalarField(char* fieldName,SF sf);
void FileIO_ReadVectorField(char* field,VF vf,int numDims);
void FileIO_ReadInts(char*,int* result,FILE* file,int num);
void FileIO_ReadDoubles(char*,double* result,FILE* file,int num);
void FileIO_ReadParameters();
void FileIO_ReadParticleData(char* name, T2 tensor[MAX_C]);
void FileIO_ReadString(FILE* file, char* result);
void FileIO_ReadTensorField(char* field,TF tf, int td2, int td1);
void FileIO_ReadT2(char* fieldName,T2 tensor, int dim1, int dim2);
void FileIO_ReadVariableSizeTF(char* field,TF tf, int td2, int* td1);
void FileIO_ReadVariableSizeTF4(char* field,TF4 tf, int td4, int* td3, int td2, int td1);
int  FileIO_SeekScalarInFile(char* fileName,double targetValue);
int  FileIO_VerifyFile(char* fileName);
void FileIO_VerifyDirectory(char* path);
void FileIO_WriteAsDouble(FILE* file, double value, char* label );
void FileIO_WriteParameters();
void FileIO_WriteParticleData(char* name, T2 tensor[MAX_C],int step);
void FileIO_WriteScalarField(char* fieldName,SF sf,int step);
void FileIO_WriteString(FILE* file, char* value, char* label);
void FileIO_WriteT2(char* field,T2 tensor,int dim1, int dim2, int step);
void FileIO_WriteTensorField(char* field,TF tf,int td2, int td1,int step);
void FileIO_WriteVariableSizeTF(char* field,TF tf, int td2, int* td1,int step);
void FileIO_WriteVariableSizeTF3(char* field,TF3 t3, int td3, int* td2, int td1,int step);
void FileIO_WriteVariableSizeTF4(char* field,TF4 TF4, int td4, int* td3, int td2,int td1, int step);
void FileIO_WriteVectorField(char* fieldName,VF vf,int numDims,int step);
void FileIO_WriteVector(char* fieldName,VEC vec,int dim,int step);

char appendDir[MAX_STRING];         // directory for vectors
char buffer[MAX_STRING];            // string buffer for reading in labels
char dataDir[MAX_STRING];           // path to top data directory
char startupDest[MAX_STRING];       // data directory for startup files
char runDir[MAX_STRING];            // path to data dir for this run
char runName[MAX_STRING];           // auto-generated path for this run
char runStatsPath[MAX_STRING];      // path to the run statistics file
char experimentFile[MAX_STRING];    // path to experiement descriptor file

SF resampleBuffer;                  // buffer for holding initial data of diff size

#endif

