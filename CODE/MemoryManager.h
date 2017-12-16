//******************************************************************************
//  MemoryManager.h
//******************************************************************************
#ifndef MemoryManager_H
#define MemoryManager_H
#include "Constants.h"

// Memory allocation routines
VEC     MemoryManager_AllocateVector(int size);
T2      MemoryManager_AllocateT2(int td2,int td1);
SF      MemoryManager_AllocateArray(int* sizes);
SF      MemoryManager_AllocateScalarField();
VF      MemoryManager_AllocateVectorField(int td1);
TF      MemoryManager_AllocateTensorField(int td2,int td1);
TF      MemoryManager_AllocateVariableSizeTF(int td2,int* td1);
T3      MemoryManager_AllocateVariableSizeT3(int td3,int* td2,int td1);
TF3     MemoryManager_AllocateVariableSizeTF3(int td3,int* td2,int td1);
TF3     MemoryManager_Allocate3TensorField(int td3, int td2,int td1);
TF4     MemoryManager_AllocateVariableSizeTF4(int td4, int* td3,int td2,int td1);
TF      MemoryManager_ExpandVariableSizeTF(TF oldTF,int tDim, int vDim);

// Memory de-allocation routines
void    MemoryManager_FreeVector        (VEC vector);
void    MemoryManager_FreeT2            (T2 tensor, int d);
void    MemoryManager_FreeArray         (SF sf, int* sizes);
void    MemoryManager_FreeScalarField   (SF sf);
void    MemoryManager_FreeVectorField   (VF vf, int d);
void    MemoryManager_FreeTensorField   (TF tf, int td2, int td1);
void    MemoryManager_FreeVariableSizeTF(TF tf,int td2,int* td1);
void    MemoryManager_Free3TensorField  (TF3 t3, int td3, int td2,int td1);
void    MemoryManager_FreeVariableSizeTF3(TF3 t3, int td3, int* td2,int td1);
void    MemoryManager_FreeVariableSizeTF4(TF4 TF4, int td4, int* td3, int td2,int td1);

#endif
