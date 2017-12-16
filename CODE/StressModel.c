//******************************************************************************
//  StressModel.c
//******************************************************************************
#include "Constants.h"
#include "FileManager.h"
#include "FFTWManager.h"
#include "HydroModel.h"
#include "MathManager.h"
#include "MemoryManager.h"
#include "SimulationManager.h"
#include "StressModel.h"
#include "ThermoModel.h"
#include "WallManager.h"

//******************************************************************************
// StressModel_MeasureSystemProperties
//******************************************************************************
void StressModel_MeasureSystemProperties()
{
    rmsVTube = MathManager_RMSVectorField(vTube);
}

//******************************************************************************
// StressModel_ComputeElasticForce
//******************************************************************************
void StressModel_ComputeElasticForce(int index)
{   
    MathManager_ClearVectorField(elasticForceSum,D);
    int c,m; FORALL_LIQUID_COMPONENTS(m,c)  
    {
        MathManager_DivTensorField(sigma[c][m],elasticForce[c][m]);
        int i,x,y,z; FOREACH_VINDEX(i,x,y,z)
        {
            elasticForceSum[i][x][y][z]+=elasticForce[c][m][i][x][y][z];
        }
    }
}

//******************************************************************************
// StressModel_ComputeUpperConvectedTerms
//******************************************************************************
void StressModel_ComputeUpperConvectedTerms()
{
    int c,m,i,j,k,x,y,z;
    FORALL_LIQUID_COMPONENTS(m,c) MathManager_ClearTensorField(ucTerms[c][m],D,D);
    
    if (useUpperConvected) 
    {
        FORALL_LIQUID_COMPONENTS(m,c) 
        {
            MathManager_VectorDotGradTensor(vTube,sigmaOld[c][m],ucTerms[c][m]);

            FOREACH_TINDEX(i,j,x,y,z) FOREACH_DIM(k)
            {
                ucTerms[c][m][i][j][x][y][z]+=
                -sigmaOld[c][m][i][k][x][y][z]*gradVTube[k][j][x][y][z]
                -gradVTube[i][k][x][y][z]*sigmaOld[c][m][k][j][x][y][z];
            }
        }
    }
}

//******************************************************************************
// StressModel_AdvanceSigma
//******************************************************************************
void StressModel_AdvanceSigma()
{
    int c,i,j,m,x,y,z;

    FOREACH_TINDEX(i,j,x,y,z)
    {
        kappa[i][j][x][y][z]=
        +gradVTube[i][j][x][y][z]
        +gradVTube[j][i][x][y][z]
        -(2.0/dim)*DELTA_FCN(i,j)*divVTube[x][y][z];
    }

    StressModel_ComputeUpperConvectedTerms();
    
    FORALL_LIQUID_COMPONENTS(m,c) FOREACH(x,y,z)
    {
        G[c][m][x][y][z]= G0[c][m]*SQR(phi[c][m][x][y][z]);
        K[c][m][x][y][z]= K0[c][m]*THETA(phi[c][m][x][y][z],f[c][m]);

        FOREACH_DIM(i) FOREACH_DIM(j)
        {
            // For the first 2 steps, make sigma part of the steady state loop
            if (timeStep>2) sigma[c][m][i][j][x][y][z]= sigmaOld[c][m][i][j][x][y][z];
            
            sigma[c][m][i][j][x][y][z]+= 
            + deltaT*G[c][m][x][y][z]*kappa[i][j][x][y][z]
            + deltaT*K[c][m][x][y][z]*divVTube[x][y][z]*DELTA_FCN(i,j)
            - deltaT*ucTerms[c][m][i][j][x][y][z];
                        
            sigma[c][m][i][j][x][y][z]/=(1.0 + deltaT/tauRep[c]);
        }
    }
}

//******************************************************************************
// StressModel_ComputeStressDivisionParams
//******************************************************************************
void StressModel_ComputeStressDivisionParams()
{
    int c,m,x,y,z; 
    FOREACH(x,y,z)
    {
        double sum = 0;
        FORALL_MONOMERS(m,c) 
        {
            zeta[c][m][x][y][z]=zeta0[c][m]*phi[c][m][x][y][z];//*Z[c];
            sum+=zeta[c][m][x][y][z];
        }
        
        FORALL_MONOMERS(m,c) alpha[c][m][x][y][z]=zeta[c][m][x][y][z]/NONZERO(sum);
    }
    FORALL_MONOMERS(m,c) rmsAlpha[c][m] = MathManager_RMSScalarField(alpha[c][m]);
}

//******************************************************************************
// StressModel_ComputeTubeVelocities
//******************************************************************************
void StressModel_ComputeTubeVelocities()
{
    int c,i,m,x,y,z;
    
    MathManager_ClearVectorField(vTemp1,D);
    MathManager_ClearVectorField(vTemp2,D);

    FORALL_LIQUID_COMPONENTS(m,c) FOREACH_VINDEX(i,x,y,z)
    {
        // sum of relative velocity w terms
        vTemp1[i][x][y][z]+=(alpha[c][m][x][y][z]-phi[c][m][x][y][z])*w[c][m][i][x][y][z];
       
        // sum of alpha terms
        vTemp2[i][x][y][z]+=(alpha[c][m][x][y][z]-phi[c][m][x][y][z]);
    }
    FORALL_SOLID_COMPONENTS(m,c) FOREACH_VINDEX(i,x,y,z)
    {
        // sum of wall terms
        vTemp1[i][x][y][z]+=(alpha[c][m][x][y][z]-phi[c][m][x][y][z])*vPhi[c][m][i][x][y][z];
    }
    
    FOREACH_VINDEX(i,x,y,z) 
    {
        vTube[i][x][y][z]=(v[i][x][y][z]+vTemp1[i][x][y][z])/NONZERO(1.0-vTemp2[i][x][y][z]);
    }
        
    MathManager_Divergence(vTube,divVTube);
    MathManager_GradVector(vTube,gradVTube);
}

//******************************************************************************
// StressModel_Finalize
//******************************************************************************
void StressModel_Finalize()
{
    MemoryManager_FreeScalarField(divVTube);
    
    MemoryManager_FreeVectorField(vTube,D);
    MemoryManager_FreeVectorField(elasticForceSum,D);
    MemoryManager_FreeVectorField(vTemp1,D);
    MemoryManager_FreeVectorField(vTemp2,D);
    MemoryManager_FreeVectorField(vTemp3,D);

    MemoryManager_FreeTensorField(gradVTube,D,D);
    MemoryManager_FreeTensorField(kappa,D,D);

    MemoryManager_FreeVariableSizeTF(alpha,numC,numM);
    MemoryManager_FreeVariableSizeTF(G,numC,numM);
    MemoryManager_FreeVariableSizeTF(K,numC,numM);
    MemoryManager_FreeVariableSizeTF(zeta,numC,numM);

    MemoryManager_FreeVariableSizeTF3(elasticForce,numC,numM,D);

    MemoryManager_FreeVariableSizeTF4(ucTerms,numC,numM,D,D);
    MemoryManager_FreeVariableSizeTF4(sigma,numC,numM,D,D);
    MemoryManager_FreeVariableSizeTF4(sigmaOld,numC,numM,D,D);
}

//******************************************************************************
// StressModel_Initialize
//******************************************************************************
void StressModel_Initialize()
{
    SimulationManager_Trace("StressModel_Initialize",2);
    
    divVTube        = MemoryManager_AllocateScalarField();
    
    elasticForceSum = MemoryManager_AllocateVectorField(D);
    vTube           = MemoryManager_AllocateVectorField(D);
    vTemp1           = MemoryManager_AllocateVectorField(D);
    vTemp2           = MemoryManager_AllocateVectorField(D);
    vTemp3           = MemoryManager_AllocateVectorField(D);

    gradVTube       = MemoryManager_AllocateTensorField(D,D);
    kappa           = MemoryManager_AllocateTensorField(D,D);
    
    alpha           = MemoryManager_AllocateVariableSizeTF(numC,numM);
    G               = MemoryManager_AllocateVariableSizeTF(numC,numM);
    K               = MemoryManager_AllocateVariableSizeTF(numC,numM);
    zeta            = MemoryManager_AllocateVariableSizeTF(numC,numM);
    
    elasticForce    = MemoryManager_AllocateVariableSizeTF3(numC,numM,D);

    sigma           = MemoryManager_AllocateVariableSizeTF4(numC,numM,D,D);
    sigmaOld        = MemoryManager_AllocateVariableSizeTF4(numC,numM,D,D);
    ucTerms         = MemoryManager_AllocateVariableSizeTF4(numC,numM,D,D);

}
