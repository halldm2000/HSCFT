//******************************************************************************
//  HydroModel_MultiFluid.c
//  Created by halldm on Fri Jul 23 2004.
//******************************************************************************
#include "Constants.h"
#include "FileIO.h"
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
// HydroModel_ApplyProjectionMethod
// Decouple P and v fields using a projection method in Fourier space
//******************************************************************************
void HydroModel_ApplyProjectionMethod()
{
    MathManager_Divergence(vStar,divVStar);
    int x,y,z; FOREACH(x,y,z) divVStar[x][y][z]*=1.0/deltaT_v;
    MathManager_SolvePoissonEquation(divVStar,P);
}

//******************************************************************************
// HydroModel_ComputeExcessVolumeFractions
//******************************************************************************
void HydroModel_ComputeVolumeFractions()
{
    int c,m,x,y,z; FORALL_MONOMERS(m,c)
    {
        double localEVF=0;
        double localAVF=0;
        FOREACH(x,y,z) localEVF+= (phi[c][m][x][y][z]>f[c][m]);
        FOREACH(x,y,z) localAVF+= (phi[c][m][x][y][z]>0.5);
        evf[c][m]=0; FIND_GLOBAL_SUM(localEVF,evf[c][m]);
        avf[c][m]=0; FIND_GLOBAL_SUM(localAVF,avf[c][m]);
        evf[c][m]*=norm;
        avf[c][m]*=norm;
    }
}

//******************************************************************************
// HydroModel_MeasureSystemProperties
//******************************************************************************
void HydroModel_MeasureSystemProperties()
{
    int c,m; 
    MathManager_VectorFieldMagnitude(v,vMag);
    
    if (appendPhi) FORALL_MONOMERS(m,c)
    {
        maxPhi[c][m]    = MathManager_Maximum(phi[c][m]);
        meanPhi[c][m]   = MathManager_Average(phi[c][m]);
        minPhi[c][m]    = MathManager_Minimum(phi[c][m]);
    }

    if (appendPhiGrads) FORALL_MONOMERS(m,c)
    {
        MathManager_Gradient(phi[c][m],gradPhi[c][m]);
        rmsGradPhi[c][m]=MathManager_RMSVectorField(gradPhi[c][m]);
    }
    
    if (appendVels) 
    {
        rmsV= MathManager_RMSVectorField(v);
        meanV= MathManager_AverageMagnitude(v);
        FORALL_MONOMERS(m,c) rmsVM[c][m] = MathManager_RMSVectorField(vPhi[c][m]);
        FORALL_MONOMERS(m,c) rmsW[c][m] = MathManager_RMSVectorField(w[c][m]);
    }
    
    if (appendResiduals) rmsVDiv = MathManager_RMSScalarField(vDiv);        
    if (appendEVF) HydroModel_ComputeVolumeFractions();
    if (appendForces)
    {
        maxElasticForce = MathManager_MaximumMagnitude(elasticForceSum);
        maxViscousForce = MathManager_MaximumMagnitude(viscousForce);
        maxOsmoticForce = MathManager_MaximumMagnitude(gradPi);
        maxWallForce    = MathManager_MaximumMagnitude(wallForces);
    }
}

//******************************************************************************
// HydroModel_Initialize
//******************************************************************************
void HydroModel_Initialize()
{
    vCountTotal=0;
    SimulationManager_Trace("HydroModel_Initialize",2);

    divVStar        = MemoryManager_AllocateScalarField();
    P               = MemoryManager_AllocateScalarField();
    vDiv            = MemoryManager_AllocateScalarField();
    vMag            = MemoryManager_AllocateScalarField();
    viscosity       = MemoryManager_AllocateScalarField();

    brownianForce   = MemoryManager_AllocateVectorField(D);
    convectiveTerm  = MemoryManager_AllocateVectorField(D);
    dVdL            = MemoryManager_AllocateVectorField(D);
    externalForce   = MemoryManager_AllocateVectorField(D);
    flux            = MemoryManager_AllocateVectorField(D);
    gradP           = MemoryManager_AllocateVectorField(D);
    gradPi          = MemoryManager_AllocateVectorField(D);
    osmoticForce    = MemoryManager_AllocateVectorField(D);
    v               = MemoryManager_AllocateVectorField(D);
    viscousForce    = MemoryManager_AllocateVectorField(D);
    vLiquid         = MemoryManager_AllocateVectorField(D);
    vErr            = MemoryManager_AllocateVectorField(D);
    vForces         = MemoryManager_AllocateVectorField(D);
    vOld            = MemoryManager_AllocateVectorField(D);
    vStar           = MemoryManager_AllocateVectorField(D);
    vWall           = MemoryManager_AllocateVectorField(D);
    wallForces      = MemoryManager_AllocateVectorField(D);
    wForces         = MemoryManager_AllocateVectorField(D);
    
    fPsi            = MemoryManager_AllocateVariableSizeTF3(numC,numM,D);
    gradMu          = MemoryManager_AllocateVariableSizeTF3(numC,numM,D);
    gradPhi         = MemoryManager_AllocateVariableSizeTF3(numC,numM,D);
    vPhi            = MemoryManager_AllocateVariableSizeTF3(numC,numM,D);
    w               = MemoryManager_AllocateVariableSizeTF3(numC,numM,D);
}

//******************************************************************************
// HydroModel_Finalize
//******************************************************************************
void HydroModel_Finalize()
{
    MemoryManager_FreeScalarField(divVStar);
    MemoryManager_FreeScalarField(P);
    MemoryManager_FreeScalarField(vDiv);
    MemoryManager_FreeScalarField(vMag);
    
    MemoryManager_FreeVectorField(brownianForce,D);
    MemoryManager_FreeVectorField(convectiveTerm,D);
    MemoryManager_FreeVectorField(dVdL,D);
    MemoryManager_FreeVectorField(externalForce,D);
    MemoryManager_FreeVectorField(flux,D);
    MemoryManager_FreeVectorField(gradP,D);
    MemoryManager_FreeVectorField(gradPi,D);
    MemoryManager_FreeVectorField(osmoticForce,D);
    MemoryManager_FreeVectorField(v,D);
    MemoryManager_FreeVectorField(viscousForce,D);
    MemoryManager_FreeVectorField(vLiquid,D);
    MemoryManager_FreeVectorField(vErr,D);
    MemoryManager_FreeVectorField(vForces,D);
    MemoryManager_FreeVectorField(vOld,D);
    MemoryManager_FreeVectorField(vStar,D);
    MemoryManager_FreeVectorField(vWall,D);
    MemoryManager_FreeVectorField(wallForces,D);
    MemoryManager_FreeVectorField(wForces,D);
    
    // TODO: free net force
    MemoryManager_FreeVariableSizeTF(phi,numC,numM);
    MemoryManager_FreeVariableSizeTF3(gradMu,numC,numM,D);
    MemoryManager_FreeVariableSizeTF3(gradPhi,numC,numM,D);
    MemoryManager_FreeVariableSizeTF3(vPhi,numC,numM,D);
    MemoryManager_FreeVariableSizeTF3(w,numC,numM,D);
}

//******************************************************************************
// HydroModel_RescalePolymerFields
//******************************************************************************
void HydroModel_RescalePolymerFields()
{
    int c,m,x,y,z;
    
    // Make sure none of them are below zero
    FORALL_LIQUID_COMPONENTS(m,c) FOREACH(x,y,z) phi[c][m][x][y][z]=MAX(0.0,phi[c][m][x][y][z]);
    
    // Make sure solid field does not exceed one
    FOREACH(x,y,z) solid[x][y][z]=MIN(1.0,solid[x][y][z]);

    // reset the mean to the expected volume fraction, if it drifts
/*    FORALL_POLYMER_COMPONENTS(m,c)
    {
        double mean=0; 
        mean = MathManager_Average(phi[c][m]);
        FOREACH(x,y,z) phi[c][m][x][y][z]+=f[c][m]-mean;
    }
*/

    // To keep chemical potentials finite make sure polymers don't hit zero or 1.
    double PHI_MIN=0;
    double PHI_MAX=1.0-PHI_MIN;
    FORALL_POLYMER_COMPONENTS(m,c) FOREACH(x,y,z) 
    {
        phi[c][m][x][y][z]=MIN(phi[c][m][x][y][z],PHI_MAX);
        phi[c][m][x][y][z]=MAX(phi[c][m][x][y][z],PHI_MIN);
    }
    
    //Compute the sum of all liquid volume fraction fields
    MathManager_ClearScalarField(liquid);
    FORALL_LIQUID_COMPONENTS(m,c) FOREACH(x,y,z) liquid[x][y][z]+=phi[c][m][x][y][z];

    //If they don't add up to the correct amount, rescale them so they do
    FORALL_LIQUID_COMPONENTS(m,c) FOREACH(x,y,z) 
    {
        phi[c][m][x][y][z]*=(1.0-solid[x][y][z])/NONZERO(liquid[x][y][z]);
    }
    
    
}

//******************************************************************************
// HydroModel_AdvancePhi
// Compute the motion of all concentration fields
//******************************************************************************
void HydroModel_AdvancePhi()
{       
    SimulationManager_Trace("HydroModel_AdvancePhi",2);

    int c,i,j,m,x,y,z;
    
    FORALL_LIQUID_COMPONENTS(m,c)
    {   
        // Transport Conserved Conventration Fields
        MathManager_SumVectorFields(w[c][m],1.0,vTube,vPhi[c][m]);    
        FOREACH_VINDEX(i,x,y,z) flux[i][x][y][z]=-phi[c][m][x][y][z]*vPhi[c][m][i][x][y][z];
        MathManager_Divergence(flux,temp1);
        FOREACH(x,y,z) phi[c][m][x][y][z]+= deltaT*temp1[x][y][z];
    }
    
    FOREACH_SOLID(c)     
    {
        FOREACH_PARTICLE(j,c) FOREACH_DIM(i) 
        {
            rCM[c][j][i]+=vCM[c][j][i]*deltaT;
            aCM[c][j][i]+=avCM[c][j][i]*deltaT;
        }
        WallManager_ConstructParticleField(c);
    }

    WallManager_ComputeLiquidFraction();
    HydroModel_RescalePolymerFields();
}

//******************************************************************************
// HydroModel_ComputeTotalOsmoticPressure
//******************************************************************************
void HydroModel_ComputeTotalOsmoticPressure()
{
    int c,i,m,x,y,z;
    MathManager_ClearVectorField(gradPi,D);
    
    FORALL_MONOMERS(m,c) FOREACH_VINDEX(i,x,y,z)
        gradPi[i][x][y][z]+= phi[c][m][x][y][z]*gradMu[c][m][i][x][y][z];
}

//******************************************************************************
// HydroModel_ComputeViscousForce
// Note: This is a temporary function for testing hydro without viscoelasticity
//******************************************************************************
void HydroModel_ComputeViscousForce()
{
    // Compute effective viscosity at this location
    MathManager_ClearScalarField(viscosity);
    int c,i,m,x,y,z; 
    FOREACH(x,y,z) FORALL_LIQUID_COMPONENTS(m,c)
        viscosity[x][y][z]+=mViscosity[c][m]*phi[c][m][x][y][z]/NONZERO(liquid[x][y][z]);
        
    FOREACH_DIM(i) 
    {
        MathManager_Laplacian(v[i],viscousForce[i]);
        int x,y,z; FOREACH(x,y,z) viscousForce[i][x][y][z]*=viscosity[x][y][z];
    }
}

//******************************************************************************
// HydroModel_ComputeResidualErrors
// measure error in velocity relative to the maximum in the problem
//******************************************************************************
void HydroModel_ComputeResidualErrors()
{
    SimulationManager_Trace("HydroModel_ComputeResidualErrors",2);

    int i,x,y,z;
    maxV = MathManager_MaximumMagnitude(v);
    double vscale = MAX(maxV,1.0e-2);
    
    FOREACH_DIM(i) FOREACH(x,y,z)
    {
        vErr[i][x][y][z]= (vForces[i][x][y][z] - gradP[i][x][y][z])/vscale;
        dVdL[i][x][y][z]= v[i][x][y][z] - vOld[i][x][y][z];
    }
    
    MathManager_Divergence(v,vDiv);
    rmsVErr = MathManager_RMSVectorField(vErr);
    rmsdVdL = MathManager_RMSVectorField(dVdL);
    rmsVDiv = MathManager_RMSScalarField(vDiv);
}

//******************************************************************************
// HydroModel_AdvanceV
// Update the coupled mean velocity and pressure fields
//******************************************************************************
void HydroModel_AdvanceV()
{    
    SimulationManager_Trace("HydroModel_AdvanceV",2);
    int c,i,m,x,y,z; 
    double invCap=1.0/Ca;
    
    if (RE>0) MathManager_VectorDotGradVector(v,v,convectiveTerm);
    
    FOREACH_DIM(i) FOREACH(x,y,z)
    {   
        double localDensity=0.0;
        FORALL_LIQUID_COMPONENTS(m,c) localDensity+=phi[c][m][x][y][z]*density[c][m];
        
        vForces[i][x][y][z]= 
            - invCap*gradPi[i][x][y][z]
            + elasticForceSum[i][x][y][z]
            + viscousForce[i][x][y][z]
            - wallForces[i][x][y][z]
            + localDensity*externalForce[i][x][y][z]
            - RE*convectiveTerm[i][x][y][z];
        
        vStar[i][x][y][z]= v[i][x][y][z] + deltaT_v*vForces[i][x][y][z];
    }
    
    HydroModel_ApplyProjectionMethod();
    MathManager_Gradient(P,gradP);

    FOREACH_DIM(i) FOREACH(x,y,z) v[i][x][y][z] = vStar[i][x][y][z] - deltaT_v*gradP[i][x][y][z];
    if (suppressVDrift) FOREACH_DIM(i)
    {
        //printf("suppressing drift.\n");
        double mean = MathManager_Average(v[i]);
        FOREACH(x,y,z) v[i][x][y][z]-=mean;
    }
}

//******************************************************************************
// HydroModel_AdvanceW
// Compute relative velocities induced by force imbalance 
//******************************************************************************
void HydroModel_AdvanceW()
{
    SimulationManager_Trace("HydroModel_AdvanceW",2);
    int c,i,m,x,y,z; 
    double InvCa=1.0/Ca;
    
    FORALL_LIQUID_COMPONENTS(m,c) FOREACH_DIM(i) FOREACH(x,y,z)
    {   

        wForces[i][x][y][z]=
        -InvCa*phi[c][m][x][y][z]*gradMu[c][m][i][x][y][z]
        +alpha[c][m][x][y][z]*elasticForceSum[i][x][y][z]
        +phi[c][m][x][y][z]*
        (-gradP[i][x][y][z]+viscousForce[i][x][y][z]+density[c][m]*externalForce[i][x][y][z]);
       
        w[c][m][i][x][y][z]= wForces[i][x][y][z]*liquid[x][y][z]/NONZERO(Gamma*zeta[c][m][x][y][z]);
    }

    FORALL_SOLID_COMPONENTS(m,c) FOREACH_DIM(i) FOREACH(x,y,z)
        w[c][m][i][x][y][z]=vPhi[c][m][i][x][y][z]-vTube[i][x][y][z];
}

//******************************************************************************
// HydroModel_ComputeWallForces
//******************************************************************************
void HydroModel_ComputeWallForces()
{
    int c,i,m,x,y,z;  
    MathManager_ClearVectorField(wallForces,D);
    double InvCa=1.0/Ca;
    
    FORALL_SOLID_COMPONENTS(m,c)     
    {        
        MathManager_ClearVectorField(fPsi[c][m],D);

        FOREACH_DIM(i) FOREACH(x,y,z) if (phi[c][m][x][y][z]!=0)
        {
            fPsi[c][m][i][x][y][z]=
                phi[c][m][x][y][z]*(-InvCa*gradMu[c][m][i][x][y][z]-gradP[i][x][y][z])
                +alpha[c][m][x][y][z]*elasticForceSum[i][x][y][z]
                +phi[c][m][x][y][z]*viscousForce[i][x][y][z]
                -Gamma*(zeta[c][m][x][y][z]/NONZERO(liquid[x][y][z]))*
                (vPhi[c][m][i][x][y][z]-vTube[i][x][y][z]);
            
           // Add solid-solid frictional force
           // if (RE>0) fPsi[c][m][i][x][y][z]+=
           //     -10.0*phi[c][m][x][y][z]*
           //     (vPhi[c][m][i][x][y][z]-vWall[i][x][y][z]);
                      
            wallForces[i][x][y][z]+=fPsi[c][m][i][x][y][z];
            fPsi[c][m][i][x][y][z]+=phi[c][m][x][y][z]*density[c][m]*externalForce[i][x][y][z];
        }
    }
}

//******************************************************************************
// HydroModel_AdvanceParticleVelocities
// TODO: include moment of intertia tensor, I
//******************************************************************************
void HydroModel_AdvanceParticleVelocities()
{
    int i,j,c; 
    double mass = 1.0e+2;
    double moment = 1.0e+3;
        
    FOREACH_SOLID(c) FOREACH_PARTICLE(j,c)
    {
        WallManager_MeasureForcesAndTorques(c,j);
        
        if (RE>0)
        {
            mass = particleMass[c];
            moment = particleMoment[c];
        }
        
        if (tConstrained[c]) 
        {
            FOREACH_DIM(i) 
            vCM[c][j][i]=transVMax[c][i]*cos((transFreq[c][i]*runTime + transPhase[c][i])*TWO_PI);
        }
        else
        {
            FOREACH_DIM(i) 
            vCM[c][j][i]+=1.0/NONZERO(mass)*deltaT_v*netForce[c][j][i];
        }
        
        if (aConstrained[c]) 
        {
             
            avCM[c][j][2]=angSpeedMax[c]*cos((angFreq[c]*runTime + angPhase[c])*TWO_PI);
        }
        else
        {
            FOREACH_DIM(i) 
            avCM[c][j][i]+=1.0/NONZERO(moment)*deltaT_v*netTorque[c][j][i];
        }
    }
}

//******************************************************************************
// HydroModel_DisplayStatus
//******************************************************************************
void HydroModel_DisplayStatus()
{
    if (vCount%vDisplayInterval==0 && timeStep%displayInterval==0) 
    {
        if (processRank==0)
            printf("vcount=%d vErr=%.3e dvdL=%.3e divV=%.3e maxV=%.3e\n",
                   vCount,rmsVErr,rmsdVdL,rmsVDiv,maxV);
        
        int c; FOREACH_SOLID(c) if (!aConstrained[c])
        {
            printf("\t torque[z]=%.2e avCM[z]=%.2e aCM[z]=%.2e\n",
                   netTorque[c][0][2],avCM[c][0][2],aCM[c][0][2]);
        }
    }
}

//******************************************************************************
// HydroModel_WriteSubsteps
//******************************************************************************
void HydroModel_WriteSubsteps()
{
    if((writeVSubsteps) && (vCountTotal%appendInterval==0))
    {
        outputStep = vCountTotal;
        FileIO_WriteVectorField("vSubstep",v,D,vCountTotal);
        FileIO_WriteVectorField("vTubeSubstep",v,D,vCountTotal);
        FileIO_WriteVectorField("eForceSubstep",elasticForceSum,D,vCountTotal);
        FileIO_WriteVectorField("vForceSubstep",viscousForce,D,vCountTotal);
        FileIO_WriteVectorField("wForceSubstep",wallForces,D,vCountTotal);
        FileIO_WriteVectorField("xForceSubstep",externalForce,D,vCountTotal);
        FileIO_WriteVectorField("gradPiSubstep",gradPi,D,vCountTotal);
        FileIO_WriteScalarField("pSubstep",P,vCountTotal);
    }    
}

//******************************************************************************
// HydroModel_ComputeSteadyState
//******************************************************************************
void HydroModel_SeekSteadyState()
{
    SimulationManager_Trace("HydroModel_SeekSteadyState",2);

    int c,m; 
    FORALL_LIQUID_COMPONENTS(m,c) MathManager_CopyTensorField(sigma[c][m],sigmaOld[c][m],D,D);
    MathManager_CopyVectorField(v,vOld,D);
        
    vCount=0;
    int maxVSteps = 1.0e4;
    if (RE>0) deltaT_v = deltaT;
    
    double residualErr=0;
    do
    {
        WallManager_ComputeWallVelocityFields();
        if (updateStresses && timeStep>0) StressModel_AdvanceSigma();
        if (updateStresses) StressModel_ComputeElasticForce();
        StressModel_ComputeStressDivisionParams();
        HydroModel_ComputeTotalOsmoticPressure();
        StressModel_ComputeTubeVelocities();
        HydroModel_ComputeWallForces();
        HydroModel_ComputeViscousForce();
        HydroModel_AdvanceW();
        HydroModel_AdvanceV();
        HydroModel_AdvanceParticleVelocities();
        HydroModel_ComputeResidualErrors();
        HydroModel_DisplayStatus();
        HydroModel_WriteSubsteps();
        vCount++; vCountTotal++;
        if (!finite(rmsVErr)) exit(1);
        if (RE==0) residualErr = rmsVErr;
    }
    while(residualErr>vThreshold && vCount<maxVSteps);
}
