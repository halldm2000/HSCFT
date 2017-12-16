//******************************************************************************
//  ThermoModel.c
//  Created by halldm on Fri Jul 23 2004.
//******************************************************************************
#include "Constants.h"
#include "FileManager.h"
#include "FFTWManager.h"
#include "HydroModel.h"
#include "MathManager.h"
#include "MemoryManager.h"
#include "SimulationManager.h"
#include "ThermoModel.h"
#include "WallManager.h"

//******************************************************************************
// ThermoModel_MeasureSystemProperties
//******************************************************************************
void ThermoModel_MeasureSystemProperties()
{
    int c,m;
    if (appendMu) FORALL_MONOMERS(m,c) rmsMu[c][m]=MathManager_RMSScalarField(mu[c][m]);
    if (appendOmegas) FORALL_MONOMERS(m,c) meanOmega[c][m]=MathManager_Average(omega[c][m]);    
}

//******************************************************************************
// ThermoModel_ApplyExponentialLaplacian
//******************************************************************************
void ThermoModel_ApplyExponentialLaplacian(double*** field, double*** result,double step, int c)
{
    int freqs[D],index;
    double k[D],k2,temp;
    MathManager_CopyScalarField(field,result);
    fftw_complex* resultTilda = (fftw_complex*)result[0][0];

    FOURIER_TRANSFORM(result,work1);
    FOREACHFREQ(freqs)
    {
        index = FINDEX(freqs);
        GETWAVENUMS(freqs,k); 	 
        k2 = SQUARE_MAG(k);
        temp = exp(-k2*step)*norm;
        resultTilda[index].re*=temp;
        resultTilda[index].im*=temp;
    }
    INVERSE_TRANSFORM(result,work1);
}

//******************************************************************************
// performQuadrature
//******************************************************************************
void performQuadrature(SF field,int m, int c)
{
    int b,s,x,y,z;
    double a1=5.0/12.0;
    double a2=13.0/12.0;
    double lambdaFrac = lambda[c]/Q[c];
    
    FOREACH(x,y,z)
    {
        field[x][y][z]=0;

        FOREACH_BLOCK(b,c)
        {
            // Compute the densities for each point on this block
            int low = b*stepsPerBlock;
            int high= (b+1)*stepsPerBlock;
            for (s=low; s<=high; s++) 
            {
                rhoTemp[s]=lineDensity[c][m][b]*q[c][s][x][y][z]*qdagger[c][s][x][y][z];
            }
            
            //Integrate the densities for this block
            double blockSum=0.0;
            for (s=low+2; s<=high-2; s++)
            {
                blockSum+=rhoTemp[s];
            }
            blockSum+=a1*(rhoTemp[low]  + rhoTemp[high]);
            blockSum+=a2*(rhoTemp[low+1]+ rhoTemp[high-1]);
            
            //Add the block sum to the total, scaled by delta S step size
            field[x][y][z]+=blockSum*deltaS[c][b];
        }
        field[x][y][z]*=lambdaFrac;
    }
}

//******************************************************************************
// ThermoModel_ComputeAuxilliaryFields
//******************************************************************************
void ThermoModel_ComputeAuxilliaryField(int c)
{
    int m; FOREACH_MONOMER(m,c) performQuadrature(phiAux[c][m],m,c);
}

//******************************************************************************
// ThermoModel_ComputeChemicalPotentials
//******************************************************************************
void ThermoModel_ComputeChemicalPotentials()
{
    SimulationManager_Trace("ThermoModel_ComputeChemicalPotentials",2);
    int c1,c2,i,x,y,z,m1,m2;
        
    FORALL_MONOMERS(m1,c1) MathManager_ClearScalarField(mu[c1][m1]);
    
    // Add the enthalpic mixing terms
    FORALL_MONOMERS(m1,c1) FORALL_MONOMERS(m2,c2) FOREACH(x,y,z)
        mu[c1][m1][x][y][z]+=chi[c1][m1][c2][m2]*phi[c2][m2][x][y][z];

    // Add simple fluid entropy terms
    double C2=SQR(CahnNumber);//Squared Cahn number
    FORALL_SIMPLE_LIQUID_COMPONENTS(m1,c1)
    {
        MathManager_Laplacian(phi[c1][m1],temp1);
        FOREACH(x,y,z) mu[c1][m1][x][y][z]-=C2*temp1[x][y][z];
    }
    
    // Add polymer entropic terms
    FORALL_POLYMER_COMPONENTS(m1,c1) FOREACH(x,y,z) mu[c1][m1][x][y][z]+=-omega[c1][m1][x][y][z]/minN;
    
    // compute the gradient
    FORALL_MONOMERS(m1,c1) MathManager_Gradient(mu[c1][m1],gradMu[c1][m1]);
    
    //Add a noise term
    if (noise>0) 
    {
        double mean=0;
        double amplitude=1.0; //thermal noise is of order kT
        FORALL_SOLID_COMPONENTS(m1,c1)
        {
            MathManager_GaussianRandomField(temp1,mean,amplitude);
            MathManager_Gradient(temp1,brownianForce);
            FOREACH_DIM(i) FOREACH(x,y,z) gradMu[c1][m1][i][x][y][z]+=brownianForce[i][x][y][z];
        }
    }
}

//******************************************************************************
// ThermoModel_OffsetDrift
//******************************************************************************
void ThermoModel_OffsetDrift()
{
    int m,c,x,y,z;
    
    FORALL_POLYMER_COMPONENTS(m,c) 
    {
        double meanOmega = MathManager_Average(omega[c][m]);
        FOREACH(x,y,z) omega[c][m][x][y][z]-=meanOmega;
    }
}

//******************************************************************************
// ThermoModel_SeekLocalEquilbrium
//******************************************************************************
void ThermoModel_SeekLocalEquilibrium()
{
    SimulationManager_Trace("ThermoModel_SeekLocalEquilibrium",2);
    scftStep= 0;
    scftErr = 1.0;
    int maxSteps = maxSCFTSteps;
    if (timeStep<2) maxSteps=1000;
    do
    {
        WAIT_FOR_ALL_PROCESSES;
        ThermoModel_OffsetDrift();
        ThermoModel_ComputeResidualErrors();
        ThermoModel_ReduceResidualErrors();
        ThermoModel_DisplayStatus();
        scftStep++;
        if (scftStep==1) scftErr1 = scftErr;
    }
    while((scftStep < maxSteps) && (scftErr > maxSCFTResidual));
}

//******************************************************************************
// ThermoModel_ComputePartitionFcn
//******************************************************************************
void ThermoModel_ComputePartitionFcn(int c)
{
    Q[c] = MathManager_Average( q[c][ns[c]-1] );
    if (!finite(Q[c])) programStatus = STOP;
}

//******************************************************************************
// ThermoModel_ComputePropagators
//******************************************************************************
void ThermoModel_ComputePropagator(int c)
{
    int m,s,x,y,z;
    MathManager_SetScalarField(q[c][0], 1.0);
    MathManager_SetScalarField(qdagger[c][ns[c]-1],1.0);
   
    double nratio=1.0*N[c]/minN;
    // Compute propagator q[c]
    for (s=1; s<ns[c]; s++)
    {
        int b = (s-1)/stepsPerBlock;
        double scaledS = deltaS[c][b]*nratio;
        double halfS=scaledS/2.0;

        MathManager_ClearScalarField(temp1);
        FOREACH(x,y,z) FOREACH_MONOMER(m,c) 
        temp1[x][y][z]+=omega[c][m][x][y][z]*lineDensity[c][m][b];
        FOREACH(x,y,z) temp1[x][y][z]=exp(-temp1[x][y][z]*halfS);
    
        FOREACH(x,y,z) temp2[x][y][z]= temp1[x][y][z]*q[c][s-1][x][y][z];
        ThermoModel_ApplyExponentialLaplacian(temp2,temp3,scaledS,c);
        FOREACH(x,y,z) q[c][s][x][y][z] = temp1[x][y][z]*temp3[x][y][z];
    }

    int symmetric=0;
    
    //TODO: generalize to all symmetric cases: triblocks, etc.
    if (numB[c]==1) symmetric=1;
    
    if (!symmetric)
    {
        // Compute copropagator q^dagger
        for (s=(ns[c]-1)-1; s>=0; s--)
        {
            int b = s/stepsPerBlock;
            double scaledS = deltaS[c][b]*nratio;
            double halfS=scaledS/2.0;

            MathManager_ClearScalarField(temp1);
            FOREACH(x,y,z) FOREACH_MONOMER(m,c) 
                temp1[x][y][z]+=omega[c][m][x][y][z]*lineDensity[c][m][b];
            FOREACH(x,y,z) temp1[x][y][z]=exp(-temp1[x][y][z]*halfS);
            
            FOREACH(x,y,z) temp2[x][y][z]= temp1[x][y][z]*qdagger[c][s+1][x][y][z];
            ThermoModel_ApplyExponentialLaplacian(temp2,temp3,scaledS,c);
            FOREACH(x,y,z) qdagger[c][s][x][y][z]= temp1[x][y][z]*temp3[x][y][z];
        }
    }
    
    if (symmetric)
    {
        for (s=0; s<ns[c]; s++) MathManager_CopyScalarField(q[c][s],qdagger[c][ns[c]-1-s]);
    }
}

//******************************************************************************
// ThermoModel_ComputeResidualError
//******************************************************************************
void ThermoModel_ComputeResidualErrors()
{
    scftErr=0;
    int c; FOREACH_POLYMER(c)
    {
        ThermoModel_ComputePropagator(c);
        ThermoModel_ComputePartitionFcn(c);
        ThermoModel_ComputeAuxilliaryField(c);
        
        int m,x,y,z; FOREACH_MONOMER(m,c) 
        {
            // Construct a bounded target concentration
            double PHI_MIN=1.0e-3;
            double PHI_MAX=1.0-PHI_MIN;
            MathManager_CopyScalarField(phi[c][m],temp1);
            FOREACH(x,y,z) 
            {
                temp1[x][y][z]=MIN(temp1[x][y][z],PHI_MAX);
                temp1[x][y][z]=MAX(temp1[x][y][z],PHI_MIN);
            }
            
            int x,y,z; FOREACH(x,y,z) 
            {
                omegaErr[c][m][x][y][z]=(phiAux[c][m][x][y][z] - temp1[x][y][z]);
            }
            rmsOmegaErr[c][m]= MathManager_RMSScalarField(omegaErr[c][m]);
            scftErr = MAX(scftErr,rmsOmegaErr[c][m]);
        }
    }
}

//******************************************************************************
// ThermoModel_DisplayStatus
//******************************************************************************
void ThermoModel_DisplayStatus()
{
    int c,m; FORALL_MONOMERS(m,c)
    {
        rmsOmega[c][m]= MathManager_RMSScalarField(omega[c][m]);
        rmsPhi[c][m]  = MathManager_RMSScalarField(phi[c][m]);
        rmsPhiAux[c][m]  = MathManager_RMSScalarField(phiAux[c][m]);
    }

    if ( ( scftStep % scftInterval == 0) && (processRank == 0) && (timeStep % displayInterval==0) )
        printf("\t%s%d %s=%5.2e %s=%5.2e %s=%5.2e %s=%f %s=%f %s=%f %s=%.2e\n",
                "scftStep#" ,scftStep,
                "omegaErrA" ,rmsOmegaErr[0][0], 
                "scftErr"   ,scftErr, 
                "rmsPhiA"   ,rmsPhi[0][0],
                "rmsAuxA"   ,rmsPhiAux[0][0],
                "stepLength",rmsStepLength,
                "omegaA"    ,rmsOmega[0][0],
                "Q[0]"      ,Q[0]);
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// ThermoModel_Finalize
//******************************************************************************
void ThermoModel_Finalize()
{
    MemoryManager_FreeVariableSizeTF(mu       ,numC,numM);
    MemoryManager_FreeVariableSizeTF(omega    ,numC,numM);
    MemoryManager_FreeVariableSizeTF(omegaErr ,numC,numM);
    MemoryManager_FreeVariableSizeTF(omegaErrOld,numC,numM);
    MemoryManager_FreeVariableSizeTF(phiAux   ,numC,numM);
    MemoryManager_FreeVariableSizeTF(q        ,numC,ns);
    MemoryManager_FreeVariableSizeTF(qdagger  ,numC,ns);
}

//******************************************************************************
// ThermoModel_Initialize
//******************************************************************************
void ThermoModel_Initialize()
{
    SimulationManager_Trace("ThermoModel_Initialize",2);

    mu          = MemoryManager_AllocateVariableSizeTF(numC,numM);
    omega       = MemoryManager_AllocateVariableSizeTF(numC,numM);
    omegaErr    = MemoryManager_AllocateVariableSizeTF(numC,numM);
    omegaErrOld = MemoryManager_AllocateVariableSizeTF(numC,numM);
    phiAux      = MemoryManager_AllocateVariableSizeTF(numC,numM);
    q           = MemoryManager_AllocateVariableSizeTF(numC,ns);
    qdagger     = MemoryManager_AllocateVariableSizeTF(numC,ns);
    
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// ThermoModel_ReduceResidualError
// Use line minimization to search for the saddle point
//******************************************************************************
void ThermoModel_ReduceResidualErrors()
{
    int c,m,x,y,z; 
    double trial = epsilon;
    // Take a small trial step
    FOREACH(x,y,z) FORALL_POLYMER_COMPONENTS(m,c) omega[c][m][x][y][z]+=trial*omegaErr[c][m][x][y][z];
    
    if (timeStep>1)
    {
        // Find the bottom of the parabola
        FORALL_POLYMER_COMPONENTS(m,c) MathManager_CopyScalarField(omegaErr[c][m],omegaErrOld[c][m]);
        ThermoModel_ComputeResidualErrors();
        
        FOREACH(x,y,z)
        {
            double r1r1 =0;
            double r2r1 =0;

            //Compute local step length
            FORALL_POLYMER_COMPONENTS(m,c) 
            {
                r1r1+=SQR(omegaErrOld[c][m][x][y][z]);
                r2r1+=omegaErr[c][m][x][y][z]*omegaErrOld[c][m][x][y][z];
            }
            double stepLength= fabs(1.0*epsilon*r1r1/NONZERO(r1r1-r2r1) -epsilon );
            temp1[x][y][z]=MIN(stepLength,maxThermoStep);
            
            // Step to the bottom of the parabola
            FORALL_POLYMER_COMPONENTS(m,c) omega[c][m][x][y][z]+=temp1[x][y][z]*omegaErr[c][m][x][y][z];
        }
        
        rmsStepLength = MathManager_RMSScalarField(temp1);
    }
}