//******************************************************************************
//  MemoryManager.c
//  Created by halldm on Fri Jul 23 2004.
//******************************************************************************
#include "MemoryManager.h"
#include "FFTWManager.h"
#include "SimulationManager.h"

VEC MemoryManager_AllocateVector(int size) {return calloc(size, sizeof(double));}
void MemoryManager_FreeVector(VEC vector) { free(vector);}

//******************************************************************************
// MemoryManager_AllocateT2
//******************************************************************************
T2  MemoryManager_AllocateT2(int td2,int td1)
{
    T2 tensor = calloc(td2,sizeof(VEC));
    int i; for(i=0; i<td2; i++) tensor[i]=MemoryManager_AllocateVector(td1);
    return tensor;
}

//******************************************************************************
// MemoryManager_FreeT2
//******************************************************************************
void MemoryManager_FreeT2(T2 tensor,int d)
{
    int i; for (i=0; i<d; i++) MemoryManager_FreeVector(tensor[d]);
    free(tensor);
}

//******************************************************************************
// MemoryManager_AllocateArray
// allocate a scalar field whose size differs from the usual size
// TODO: make this work in MPI environment
//******************************************************************************
SF MemoryManager_AllocateArray(int* size)
{
    int totalSize = (size[0]*size[1]*size[2]);
    printf("Allocating buffer array. size[0]=%d, size[1]=%d, size[2]=%d, Totalsize = %d.\n",
           size[0],size[1],size[2],totalSize);

    double*  dataBlock = calloc(totalSize, sizeof(double));
    
    // index the pointers into the data for easy access
    SF sf = calloc(size[0], sizeof(double**));
    int x,y; for (x=0; x<size[0]; x++)
    {
        sf[x] = calloc(size[1], sizeof(double*));
        
        for (y=0; y<size[1]; y++) 
        {
            int index = size[2]*(size[1]*x + y);
            sf[x][y] = &(dataBlock[index]);
        }
    }
    return sf;
}

//******************************************************************************
// MemoryManager_AllocateScalarField
//******************************************************************************
SF MemoryManager_AllocateScalarField()
{
    // Allocate a block of pointers and a block of data
    double*  dataBlock = calloc(totalSize, sizeof(double));
    if (dataBlock == NULL) 
    {
        SimulationManager_DisplayMessage("Failed to allocate memory block.");
        exit(1);
    }
    
    // index the pointers into the data for easy access
    SF sf  = calloc(nxLocal, sizeof(double**));
    int x,y; FORX(x) 
    {
        sf[x] = calloc(ny, sizeof(double*));
        FORY(y) sf[x][y] = &(dataBlock[INDEX(x,y,0)]);
    }

    return sf;
}

//******************************************************************************
// MemoryManager_AllocateVectorField
//******************************************************************************
VF MemoryManager_AllocateVectorField(int d)
{
    VF vf = calloc(d,sizeof(SF));
    int i; for(i=0; i<d; i++) vf[i]=MemoryManager_AllocateScalarField();
    return vf;
}


//******************************************************************************
// MemoryManager_AllocateTensorField
//******************************************************************************
TF MemoryManager_AllocateTensorField(int td2,int td1)
{
    TF tf = calloc(td2,sizeof(VF));
    int i; for(i=0; i<td2; i++) tf[i]=MemoryManager_AllocateVectorField(td1);
    return tf;
}

//******************************************************************************
// MemoryManager_AllocateVariableSizeTF
//******************************************************************************
TF  MemoryManager_AllocateVariableSizeTF(int td2,int* td1)
{
    TF tf = calloc(td2,sizeof(VF));
    int i; for(i=0; i<td2; i++) tf[i]=MemoryManager_AllocateVectorField(td1[i]);
    return tf;
}

//******************************************************************************
// MemoryManager_Allocate3TensorField
//******************************************************************************
TF3 MemoryManager_Allocate3TensorField(int td3, int td2,int td1)
{
    TF3 t3 = calloc(td3,sizeof(TF));
    int i; for(i=0; i<td3; i++) t3[i]=MemoryManager_AllocateTensorField(td2,td1);
    return t3;
}

//******************************************************************************
// MemoryManager_AllocateVariableSizeT3
//******************************************************************************
T3 MemoryManager_AllocateVariableSizeT3(int td3, int* td2,int td1)
{
    T3 t3 = calloc(td3,sizeof(T2));
    int i,j; for(i=0; i<td3; i++)
    {
        t3[i]=calloc(td2[i],sizeof(VEC));
        for(j=0; j<td2[i]; j++) t3[i][j]=calloc(td1,sizeof(double));
    }
    return t3;
}

//******************************************************************************
// MemoryManager_AllocateVariableSizeTF3
//******************************************************************************
TF3 MemoryManager_AllocateVariableSizeTF3(int td3, int* td2,int td1)
{
    TF3 t3 = calloc(td3,sizeof(TF));
    int i; for(i=0; i<td3; i++) t3[i]=MemoryManager_AllocateTensorField(td2[i],td1);
    return t3;
}

//******************************************************************************
// MemoryManager_Allocate4TensorField
//******************************************************************************
TF4  MemoryManager_AllocateVariableSizeTF4(int td4,int* td3, int td2,int td1)
{
    TF4 TF4 = calloc(td4,sizeof(TF3));
    int i; for(i=0; i<td4; i++) TF4[i]=MemoryManager_Allocate3TensorField(td3[i],td2,td1);
    return TF4;
}

//******************************************************************************
// MathManager_ExpandVariableSizeTF
// Append a new vector field onto the end of an existing tensor field
//******************************************************************************
TF MemoryManager_ExpandVariableSizeTF(TF oldTF,int tDim, int vDim)
{
    TF newTF =0;
    
    if (tDim==1)
    {
        printf("Constructing new Variable Size TF. tdim=%d vdim=%d\n",tDim,vDim);
        newTF = calloc(1,sizeof(VF));
    }
    
    if (tDim>1)
    {
        printf("Expanding Variable Size to have tdim=%d vdim=%d\n",tDim,vDim);
        newTF = calloc(tDim,sizeof(VF));
        
        // Copy data pointers from the old tensor to a new larger one
        int i; for (i=0; i<tDim-1; i++) newTF[i]=oldTF[i];
        free(oldTF);
    }
    
    // Allocate the new vector of data at the end of the tensor
    newTF[tDim-1]=MemoryManager_AllocateVectorField(vDim);
    return newTF;
}

//******************************************************************************
// MemoryManager_FreeScalarField
//******************************************************************************
void MemoryManager_FreeArray(SF sf,int* size)
{
    // free data block
    free(sf[0][0]);
    
    // free pointers into data
    int x; for (x=0; x<size[0]; x++) free(sf[x]);
    free(sf);
}


//******************************************************************************
// MemoryManager_FreeScalarField
//******************************************************************************
void MemoryManager_FreeScalarField(SF sf)
{
    // free data block
    free(sf[0][0]);
    
    // free pointers into data
    int x; FORX(x) free(sf[x]);
    free(sf);
}

//******************************************************************************
// MemoryManager_FreeVectorField
//******************************************************************************
void MemoryManager_FreeVectorField(VF vf, int d)
{
    int i; for(i=0; i<d; i++) MemoryManager_FreeScalarField(vf[i]);
    free(vf);
}

//******************************************************************************
// MemoryManager_FreeTensorField
//******************************************************************************
void MemoryManager_FreeTensorField(TF tf,int td2,int td1)
{
    int i; for(i=0; i<td2; i++) MemoryManager_FreeVectorField(tf[i],td1);
    free(tf);
}

//******************************************************************************
// MemoryManager_FreeTensorField
//******************************************************************************
void MemoryManager_FreeVariableSizeTF(TF tf,int td2,int* td1)
{
    int i; for(i=0; i<td2; i++) MemoryManager_FreeVectorField(tf[i],td1[i]);
    free(tf);
}

//******************************************************************************
// MemoryManager_Free3TensorField
//******************************************************************************
void MemoryManager_Free3TensorField(TF3 t3,int td3,int td2,int td1)
{
    int i; for(i=0; i<td3; i++) MemoryManager_FreeTensorField(t3[i],td2,td1);
    free(t3);
}

//******************************************************************************
// MemoryManager_FreeVariablSizeTF3
//******************************************************************************
void MemoryManager_FreeVariableSizeTF3(TF3 t3,int td3,int* td2,int td1)
{
    int i; for(i=0; i<td3; i++) MemoryManager_FreeTensorField(t3[i],td2[i],td1);
    free(t3);
}

//******************************************************************************
// MemoryManager_Free4TensorField
//******************************************************************************
void MemoryManager_FreeVariableSizeTF4(TF4 TF4,int td4,int* td3,int td2,int td1)
{
    int i; for(i=0; i<td4; i++) MemoryManager_Free3TensorField(TF4[i],td3[i],td2,td1);
    free(TF4);
}
