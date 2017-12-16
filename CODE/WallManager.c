//******************************************************************************
//  WallManager.c
//******************************************************************************
#include <string.h>
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

#define GEOMETRY_DIR "GEOMETRY_FILES"
#define UNINITIALIZED   -1
//******************************************************************************
// WallManager_ComputePivotPoint 
// Compute the center of mass of a solid copolymer object
//******************************************************************************
void WallManager_ComputePivotPoint(int c)
{
    SimulationManager_Trace("WallManager_ComputeCenterOfMass",4);
    double total = 0;
    int i; FOREACH_DIM(i) pivot[c][i]=0.0;
    int m; FOREACH_MONOMER(m,c)
    {
        total+=MathManager_IntegrateScalarField(phi[c][m]);
        int x,y,z; FOREACH(x,y,z)
        {
            pivot[c][0]+=x*particle[c][x][y][z];
            pivot[c][1]+=y*particle[c][x][y][z];
            pivot[c][2]+=z*particle[c][x][y][z];
        }
        FOREACH_DIM(i) pivot[c][i]/=total;
        particleMass[c] = total*density[c][0];
        printf("Pivot = %.1f, %.1f, %.1f\n",pivot[c][0],pivot[c][1],pivot[c][2]);
        printf("Particle mass = %.1f\n",particleMass[c]);
    }
}

void WallManager_ComputeParticleMoment(int c)
{
    particleMoment[c]=0;
    // Sum the portion of the forces and torques which act on this object
    int x,y,z; FOREACH(x,y,z) if (particle[c][x][y][z]>0)
    {
        rprime[0]=fmod(L[0]*(x-pivot[c][0])/nx + L[0]*1.5 ,L[0]) - L[0]*0.5;
        rprime[1]=fmod(L[1]*(y-pivot[c][1])/ny + L[1]*1.5 ,L[1]) - L[1]*0.5;
        rprime[2]=fmod(L[2]*(z-pivot[c][2])/nz + L[2]*1.5 ,L[2]) - L[2]*0.5;
        double r2 = SQR(rprime[0]) + SQR(rprime[1]) + SQR(rprime[2]);
        particleMoment[c]+=density[c][0]*particle[c][x][y][z]*r2;
    }
    printf("Particle moment = %.1f\n",particleMoment[c]);

}

//******************************************************************************
// WallManager_ComputeExternalForce
// Construct an external force which induces flow in the system
//******************************************************************************
void WallManager_ComputeExternalForce()
{
    SimulationManager_Trace("WallManager_ComputeExternalForce",2);

    double unitVector[D];
    unitVector[0]=0.0; unitVector[1]=-1.0; unitVector[2]=0.0;

    int i,x,y,z; FOREACH(x,y,z) FOREACH_DIM(i)
    {
        externalForce[i][x][y][z]=efMagnitude*unitVector[i];
    }
}

//******************************************************************************
// WallManager_ComputeWallVelocityFields
// Compute externally driven velocity fields as specified by motion descriptors
//******************************************************************************
void WallManager_ComputeWallVelocityFields()
{
    SimulationManager_Trace("WallManager_ComputeWallVelocityFields",2);
    MathManager_ClearVectorField(vWall,D);

    int c,i,j,m,x,y,z;
    FORALL_SOLID_COMPONENTS(m,c)
    {
        MathManager_ClearVectorField(vPhi[c][m],D);

        FOREACH_PARTICLE(j,c)
        {
            // Add translational motion to solid velocity field
            WallManager_GetParticleField(c,j,temp1);
            FOREACH_DIM(i) FOREACH(x,y,z) if (temp1[x][y][z]!=0) vPhi[c][m][i][x][y][z]+=vCM[c][j][i];
        
            // Add rotational motion about z axis to solid velocity field
            FOREACH(x,y,z) if (temp1[x][y][z]>0)
            {
                // convert coord to a point withing L/2 of the pivot point
                rprime[0]=fmod(L[0]*x/nx - rCM[c][j][0] + L[0]*1.5 ,L[0]) - L[0]*0.5;
                rprime[1]=fmod(L[1]*y/ny - rCM[c][j][1] + L[1]*1.5 ,L[1]) - L[1]*0.5;
                rprime[2]=fmod(L[2]*z/nz - rCM[c][j][2] + L[2]*1.5 ,L[2]) - L[2]*0.5;

                vPhi[c][m][0][x][y][z]+=-avCM[c][j][2]*rprime[1];
                vPhi[c][m][1][x][y][z]+=+avCM[c][j][2]*rprime[0];
            }
            
            // Compute average wall velocity field
            FOREACH_DIM(i) FOREACH(x,y,z) 
                vWall[i][x][y][z]+=phi[c][m][x][y][z]*vPhi[c][m][i][x][y][z];
        }
    }
}

//******************************************************************************
// WallManager_GetParticleField
//******************************************************************************
void WallManager_GetParticleField(int c, int j, SF field)
{
    MathManager_RotateField(particle[c],temp3,aCM[c][j],pivot[c],bounds[c],&rbounds);
    MathManager_TranslateField(temp3,field,rCM[c][j],pivot[c],rbounds);
}

//******************************************************************************
// WallManager_GetParticlePEField
//******************************************************************************
void WallManager_GetParticlePEField(int c, int j, SF field)
{
    MathManager_RotateField(particlePE[c],temp3,aCM[c][j],pivot[c],peBounds[c],&rbounds);
    MathManager_TranslateField(temp3,field,rCM[c][j],pivot[c],rbounds);
}

//******************************************************************************
// WallManager_ConstructParticleField
//******************************************************************************
void WallManager_ConstructParticleField(int c)
{
    MathManager_ClearScalarField(phi[c][0]);
    int j; FOREACH_PARTICLE(j,c) if (rCM[c][j][0]!=UNINITIALIZED) 
    {
        WallManager_GetParticleField(c,j,temp1);
        int x,y,z; FOREACH(x,y,z) phi[c][0][x][y][z]+=temp1[x][y][z];
    }
}

//******************************************************************************
// WallManager_ComputeSolidField
//******************************************************************************
void WallManager_ComputeSolidField()
{
    int c,j,x,y,z; 
    MathManager_ClearScalarField(solid);
    MathManager_ClearScalarField(netParticlePE);
    
    FOREACH_SOLID(c) FOREACH_PARTICLE(j,c) if (rCM[c][j][0]!=UNINITIALIZED) 
    {
        WallManager_GetParticleField(c,j,temp1);
        WallManager_GetParticlePEField(c,j,temp2);
        FOREACH(x,y,z) 
        {
            netParticlePE[x][y][z]+=temp2[x][y][z];
            solid[x][y][z]+=temp1[x][y][z];
        }
    }

    fSolid = MathManager_Average(solid);
    //printf("fSolid=%f\n",fSolid);
}

//******************************************************************************
// WallManager_ComputeLiquidFraction
//******************************************************************************
void WallManager_ComputeLiquidFraction()
{
    SimulationManager_Trace("WallManager_ComputeLiquidFraction",2);
    int m,c,x,y,z;
    MathManager_ClearScalarField(liquid);
    FORALL_LIQUID_COMPONENTS(m,c) FOREACH(x,y,z) liquid[x][y][z]+=phi[c][m][x][y][z];
    WallManager_ComputeSolidField();
    fLiquid = MathManager_Average(liquid);
    //printf("fLiquid=%e\n",fLiquid);
}

//******************************************************************************
// WallManager_Initialize
//******************************************************************************
void WallManager_Initialize()
{
    SimulationManager_Trace("WallManager_Initialize",2);
    WallManager_AllocateMemory();
    WallManager_MakeGaussianFilter();
    WallManager_ConstructSingleParticlePotential();
    int i,j,c; FOREACH_MATERIAL(c) FOREACH_PARTICLE(j,c) FOREACH_DIM(i) vCM[c][j][i]=0;
}

//******************************************************************************
// WallManager_AllocateMemory
//******************************************************************************
void WallManager_AllocateMemory()
{
    filter  = MemoryManager_AllocateScalarField();
    liquid  = MemoryManager_AllocateScalarField();
    mask    = MemoryManager_AllocateScalarField();
    netParticlePE = MemoryManager_AllocateScalarField();
    potential = MemoryManager_AllocateScalarField();
    solid   = MemoryManager_AllocateScalarField();
    temp1   = MemoryManager_AllocateScalarField();
    temp2   = MemoryManager_AllocateScalarField();
    temp3   = MemoryManager_AllocateScalarField();
    total   = MemoryManager_AllocateScalarField();
    
    particle    = MemoryManager_AllocateVectorField(MAX_C);
    particlePE  = MemoryManager_AllocateVectorField(MAX_C);
}

//******************************************************************************
// WallManager_AllocateMemoryForParticle
//******************************************************************************
void WallManager_AllocateMemoryForParticle(int c,int numP)
{
    aCM[c]      = MemoryManager_AllocateT2(numP,D);
    avCM[c]     = MemoryManager_AllocateT2(numP,D);
    netTorque[c]= MemoryManager_AllocateT2(numP,D);
    netForce[c] = MemoryManager_AllocateT2(numP,D);
    rCM[c]      = MemoryManager_AllocateT2(numP,D);
    vCM[c]      = MemoryManager_AllocateT2(numP,D);
}

//******************************************************************************
// WallManager_FreeMemoryForParticle
//******************************************************************************
void WallManager_FreeMemoryForParticle(int c)
{
    MemoryManager_FreeT2(aCM[c],numParticles[c]);
    MemoryManager_FreeT2(avCM[c],numParticles[c]);
    MemoryManager_FreeT2(netForce[c],numParticles[c]);
    MemoryManager_FreeT2(netTorque[c],numParticles[c]);
    MemoryManager_FreeT2(rCM[c],numParticles[c]);
    MemoryManager_FreeT2(vCM[c],numParticles[c]);
}

//******************************************************************************
// WallManager_FreeMemory
//******************************************************************************
void WallManager_Finalize()
{
    MemoryManager_FreeScalarField(filter);
    MemoryManager_FreeScalarField(liquid);
    MemoryManager_FreeScalarField(mask);
    MemoryManager_FreeScalarField(netParticlePE);
    MemoryManager_FreeScalarField(potential);
    MemoryManager_FreeScalarField(solid);
    MemoryManager_FreeScalarField(temp1);
    MemoryManager_FreeScalarField(temp2);
    MemoryManager_FreeScalarField(temp3);
    MemoryManager_FreeScalarField(total);
    
    MemoryManager_FreeVectorField(particle,MAX_C);
    MemoryManager_FreeVectorField(particlePE,MAX_C);
    
    int c; FOREACH_SOLID(c) WallManager_FreeMemoryForParticle(c);
}

//******************************************************************************
// ConstructSingleParticlePotential
//******************************************************************************
void WallManager_ConstructSingleParticlePotential()
{
    int x,y,z;
    double r_c = 2.0*filterWidth*nx;
    double peMax = 0.1;
    FOREACH(x,y,z)
    {
        double xshift = fmod(x+nx/2.0,nx);
        double yshift = fmod(y+ny/2.0,ny);
        double zshift = fmod(z+nz/2.0,nz);
        
        double r2 = SQR(xshift-nx/2.0) + SQR(yshift-ny/2.0) + SQR(zshift-nz/2.0);
        double r = sqrt(r2);
        potential[x][y][z]=peMax*SQR(r/r_c - 1)*(r<r_c);
    }
}

//******************************************************************************
// WallManager_MakeGaussianFilter
// Construct a gaussian filter that can be used to smooth geometry masks
//******************************************************************************
void WallManager_MakeGaussianFilter()
{
    int x,y,z;
    double width2 = SQR(filterWidth*nx);
    
    FOREACH(x,y,z)
    {
        double xshift = fmod(x+nx/2.0,nx);
        double yshift = fmod(y+ny/2.0,ny);
        double zshift = fmod(z+nz/2.0,nz);
        
        double d2 = SQR(xshift-nx/2.0) + SQR(yshift-ny/2.0) + SQR(zshift-nz/2.0);
        filter[x][y][z]=exp(-d2/(2.0*width2));
    }
    double sum = MathManager_IntegrateScalarField(filter);
    FOREACH(x,y,z) filter[x][y][z]*=1.0/sum;
}

//******************************************************************************
// WallManager_ConstructSystemFromFile
//******************************************************************************
void WallManager_ConstructSystemFromFile()
{   
    int i,m,x,y,z;
    SimulationManager_Trace("WallManager_ConstructSystemFromFile",2);
    numLiquids=0;
    numSolids=0;
    
    FILE* file=0; if (processRank==0) file = fopen(systemFile,"r");
    printf("WallManager: Construct system from %s\n",systemFile);
    
    // Read in each liquid composition
    int numL; FileIO_ReadInts(buffer,&numL,file,1);
    for(i=0; i<numL; i++)
    {
        FileIO_ReadString(file, compositionFile);
        FileIO_ReadString(file, geometryFile);
        WAIT_FOR_ALL_PROCESSES;
        WallManager_FillGeometry(geometryFile,compositionFile);
        numLiquids+=copolymersInComp;
    }
    
    // Read in each solid composition
    int numW; FileIO_ReadInts(buffer,&numW,file,1);
    for(i=0; i<numW; i++)
    {
        FileIO_ReadString(file, compositionFile);
        FileIO_ReadString(file, geometryFile);
        FileIO_ReadString(file, motionDescriptorFile);
        FileIO_ReadInts(buffer,&numParticles[numC],file,1);
        
        //Dynamically allocate space for all particles of this type
        WallManager_AllocateMemoryForParticle(numC, numParticles[numC]);
        
        // Construct a single field for each solid object type
        WallManager_FillGeometry(geometryFile,compositionFile);
        MathManager_CopyScalarField(mask,particle[numC-1]);
        WallManager_ComputePivotPoint(numC-1);
        WallManager_ComputeParticleMoment(numC-1);
        WallManager_ConstructParticlePESurface(numC-1);

        numSolids+=copolymersInComp;
        
        int c; for(c=numC-copolymersInComp; c<numC; c++)
        {
            materialType[c]=SOLID;
            FileManager_ReadMotionDescriptor(c);
        }
        bounds[numC-1]=MathManager_ComputeBoundingBox(particle[numC-1]);
        
        WallManager_SetObjectLocations(numC-1);
        WallManager_ConstructParticleField(numC-1);
        
        FORALL_LIQUID_COMPONENTS(m,c)
        {
            FOREACH(x,y,z) phi[c][m][x][y][z]*=(1.0-phi[numC-1][0][x][y][z]);
        }
    }
    
    if (processRank ==0) fclose(file);
    FileManager_LookUpFloryParameters();
}

//******************************************************************************
// WallManager_MeasureForcesAndTorques
// Measure net force and toque on particle j of material c
//******************************************************************************
void WallManager_MeasureForcesAndTorques(int c,int j)
{
    int i,x,y,z;
    double torque[3];
    double force[3];
    
    FOREACH_DIM(i) netForce[c][j][i]=0.0;
    FOREACH_DIM(i) netTorque[c][j][i]=0.0;
    WallManager_GetParticleField(c,j,temp1);
    WallManager_GetParticlePEField(c,j,temp3);

    // Compute Solid-Solid chemical potential contribution
    MathManager_ClearScalarField(temp2);
    FOREACH(x,y,z) temp2[x][y][z] = netParticlePE[x][y][z] - temp3[x][y][z];
    MathManager_Gradient(temp2,vTemp2);
 
    // Compute particle-solid contribution to the osmotic pressure
    FOREACH_DIM(i) FOREACH(x,y,z) 
    {
        vTemp1[i][x][y][z]=-temp1[x][y][z]*vTemp2[i][x][y][z]/Ca;
    }
    
    // Sum the portion of the forces and torques which act on this object
    FOREACH(x,y,z) if (temp1[x][y][z]>0)
    {
        FOREACH_DIM(i) 
        {
            force[i] = fPsi[c][0][i][x][y][z] + vTemp1[i][x][y][z];
            netForce[c][j][i]+= force[i];
        }
        
        rprime[0]=fmod(L[0]*x/nx - rCM[c][j][0] + L[0]*1.5 ,L[0]) - L[0]*0.5;
        rprime[1]=fmod(L[1]*y/ny - rCM[c][j][1] + L[1]*1.5 ,L[1]) - L[1]*0.5;
        rprime[2]=fmod(L[2]*z/nz - rCM[c][j][2] + L[2]*1.5 ,L[2]) - L[2]*0.5;
        
        MathManager_CrossProduct(rprime,force,torque);
        FOREACH_DIM(i) netTorque[c][j][i]+=torque[i];
    }
    
   // FOREACH_DIM(i)
   //    {
   //        netTorque[c][j][i]*=volume*norm;
   //        netForce[c][j][i]*=volume*norm;
   //    }
}

//******************************************************************************
// WallManager_SetObjectLocations
// Choose a random location for the new object
//******************************************************************************
void WallManager_SetObjectLocations(int c)
{
    printf("Choosing Particle Locations\n");

    int j,x,y,z; 
    
    // Use x coord to indicate they should not be included in solid field.
    FOREACH_PARTICLE(j,c) rCM[c][j][0]=UNINITIALIZED;
    
    FOREACH_PARTICLE(j,c)
    {
        WallManager_ComputeSolidField();
        rCM[c][j][0]=L[0]*pivot[c][0]/gridSize[0];
        rCM[c][j][1]=L[1]*pivot[c][1]/gridSize[1];

        double collisions = 1;
        if (tConstrained[c] && aConstrained[c]) collisions=0;
        
        while (collisions>0)
        {
            if (numParticles[c]>1)
            {
                // Choose a random center of mass position
                rCM[c][j][0] = L[0]*(1.0*rand()/RAND_MAX);
                rCM[c][j][1] = L[1]*(1.0*rand()/RAND_MAX);
            }
            if (!aConstrained[c]) aCM[c][j][2] = TWO_PI*(1.0*rand()/RAND_MAX);

            // Compute collisions between this particle and the other solid objects
            collisions=0.0;
            WallManager_GetParticleField(c,j,temp1);
            FOREACH(x,y,z) collisions += temp1[x][y][z]*solid[x][y][z];
            printf("c=%d j=%d x=%.2f y=%.2f collisions = %f\n",c,j,rCM[c][j][0],rCM[c][j][1],collisions);
        }
    }
}

//******************************************************************************
// WallManager_ConstructMask
//******************************************************************************
void WallManager_ConstructMask(char* geometryFile)
{
    SimulationManager_Trace("WallManager_ConstructMask",4);
    int c,m,x,y,z;
    WallManager_MakeGeometryFromFile(mask,geometryFile);
    WallManager_SmoothEdges(mask);
    if(numParticles[numC]==0)
    FORALL_LIQUID_COMPONENTS(m,c) FOREACH(x,y,z) phi[c][m][x][y][z]*=(1.0-mask[x][y][z]);
}

//******************************************************************************
// WallManager_MakeNoisyData
//******************************************************************************
void WallManager_MakeNoisyData(SF field, double mean)
{
    SimulationManager_Trace("WallManager_MakeNoisyData",4);
    double amplitude = 1.0e-4;
    printf("target mean = %.3e\t amplitude=%.3e\n",mean,amplitude);

    MathManager_GaussianRandomField(field,mean,amplitude);
    //double order = 10.0;
    //MathManager_ApplyExponentialFilter(field,order);
    
    double result = MathManager_Average(field);
    printf("Resultant mean is:%e\n",result);
}

//******************************************************************************
// WallManager_FitDataInMask
//******************************************************************************
void WallManager_FitDataInMask(TF field,int numComponents, int total)
{
    SimulationManager_Trace("WallManager_FitDataInMask",4);

    
    int c,m,x,y,z; FOREACH(x,y,z)
    {
        double sum=0;
        for(c=total-numComponents; c<total; c++) FOREACH_MONOMER(m,c) 
            sum+=field[c][m][x][y][z];
        
        for(c=total-numComponents; c<total; c++) FOREACH_MONOMER(m,c) 
            field[c][m][x][y][z]*=mask[x][y][z]/sum;
    }
}

//******************************************************************************
// WallManager_ComputeMonomerFractionsOnCopolymer
//******************************************************************************
void WallManager_ComputeMonomerFractionsOnCopolymer(int c)
{
    SimulationManager_Trace("WallManager_ComputeMonomerFractionsOnCopolymer",4);
    
    //Rescale block lengths to add up to N[c]
    double blockSum=0;
    int b; 
    FOREACH_BLOCK(b,c) blockSum+=blockLength[c][b];
    FOREACH_BLOCK(b,c) blockLength[c][b]*=N[c]/blockSum;
    
    int m; FOREACH_MONOMER(m,c)
    {
        fM[c][m]=0.0;
        int b; FOREACH_BLOCK(b,c) fM[c][m]+=lineDensity[c][m][b]*blockLength[c][b];
        fM[c][m]*=lambdaC[c]/N[c];
        if (processRank == 0) printf("\t fM[%d][%d] \t=\t %f\n",c,m,fM[c][m]);
    }
}

//******************************************************************************
// WallManager_RescaleCopolymerFractions
// Rescale copolymer fractions to add up to one in this composition region
//******************************************************************************
void WallManager_RescaleCopolymerFractions()
{
    int c;
    double totalLambdaC=0;
    for (c=numC; c<numC+copolymersInComp; c++) totalLambdaC+=lambdaC[c];
    for (c=numC; c<numC+copolymersInComp; c++) lambdaC[c]/=totalLambdaC;
}


//******************************************************************************
// WallManager_FillGeometry
//******************************************************************************
void WallManager_FillGeometry(char* geometryFile,char* compositionFile)
{
    SimulationManager_Trace("WallManager_FillGeometry",3);
    
    int c,m;
    WallManager_ConstructMask(geometryFile);
    FileManager_ReadCompositionFile(compositionFile);
    WallManager_RescaleCopolymerFractions();
    
    // Fill the mask with the new composition
    int total = numC + copolymersInComp;
    for(c=numC; c<total; c++, numC++)
    {
        FileManager_ReadCopolymerFile(c);
        FileManager_ReadMaterialFiles(c);

        phi = MemoryManager_ExpandVariableSizeTF(phi,c+1,numM[c]);
        WallManager_ComputeMonomerFractionsOnCopolymer(c);
        FOREACH_MONOMER(m,c) WallManager_MakeNoisyData(phi[c][m],fM[c][m]);
    }
    WallManager_FitDataInMask(phi,copolymersInComp,numC);    
}

//******************************************************************************
// WallManager_ReadGeometryFile
//******************************************************************************
void WallManager_MakeGeometryFromFile(SF geom,char* fileName)
{
    SimulationManager_Trace("WallManager_MakeGeometryFromFile",5);
    char operation[MAX_STRING];
    char buffer[MAX_STRING];
    
    MathManager_ClearScalarField(geom);
    
    // Read geometry file
    printf("%d: reading %s\n",processRank,fileName);
    FILE* file=0; if (processRank ==0) file = fopen(fileName,"r");
    int numOps; FileIO_ReadInts(buffer,&numOps,file,1);
        
    //Perform geometry construction operations
    int i; for (i=0; i<numOps; i++)
    {
        double value; 
        FileIO_ReadDoubles(operation,&value,file,1);
        if (strcmp(operation,"BOX")==0)
        {
            FileIO_ReadDoubles(buffer,center,file,D);
            FileIO_ReadDoubles(buffer,width,file,D);
            WallManager_FillBox(geom,value,center,width);
        }
        if (strcmp(operation,"SPHERE")==0)
        {
            FileIO_ReadDoubles(buffer,center,file,D);
            FileIO_ReadDoubles(buffer,&radius,file,1);
            WallManager_FillSphere(geom,value,center,radius);
        }
        if (strcmp(operation,"CYLINDER")==0)
        {
            FileIO_ReadDoubles(buffer,center,file,D);
            FileIO_ReadDoubles(buffer,&radius,file,1);
            WallManager_FillCylinder(geom,value,center,radius);
        }        
        
        if (strcmp(operation,"TRIANGLE")==0)
        {
            FileIO_ReadDoubles(buffer,center,file,D);
            FileIO_ReadDoubles(buffer,width,file,D);
            WallManager_FillTriangle(geom,value,center,width);
        }
        
        if (strcmp(operation,"BITMAP")==0)
        {
            FileIO_ReadString(file,bitmapPath);
            FileIO_ReadDoubles(buffer,width,file,D);
            FileIO_ReadDoubles(buffer,center,file,D);
            WallManager_ImportBitmap();
        }
    }
    if (processRank ==0) fclose(file);
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// WallManager_ImportBitmap
//******************************************************************************
void WallManager_ImportBitmap()
{
    printf("Importing bitmapped text file:%s\n",bitmapPath);
    FILE* file = fopen(bitmapPath,"r");
    int xlow = (center[0]*nx)-width[0]/2;
    int ylow = (center[1]*ny)-width[1]/2;
    if (file!=NULL)  
    {
        int x,y; for (y=ylow+width[1]-1; y>=ylow; y--) for(x=xlow; x<xlow + width[0]; x++) 
        {
            fscanf(file,"%le", &(mask[x][y][0]) );
            mask[x][y][0]*=0.95/255;
        }
        fclose(file);
    }
    else 
    {
        printf("NULL file pointer.\n"); 
        exit(1);
    }
}


//******************************************************************************
// WallManager_ConstructParticlePESurface
//******************************************************************************
void WallManager_ConstructParticlePESurface(int c)
{
    int x,y,z;
    printf("%d: Construction particle-particle PE Surface %d\n",processRank,c);
        
    MathManager_CopyScalarField(particle[c],temp1);
    MathManager_CopyScalarField(potential,temp2);
    MathManager_Convolve(temp1,temp2,particlePE[c]);
    FOREACH(x,y,z) particlePE[c][x][y][z]*= particlePE[c][x][y][z]>1.0e-4;

    peBounds[c]=MathManager_ComputeBoundingBox(particlePE[c]);

    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// WallManager_SmoothEdges
// Apply a Gaussian filter to smooth the object, then rescale it to its maximum
//******************************************************************************
void WallManager_SmoothEdges(SF field)
{
    int x,y,z;
    printf("%d: Smoothing Edges. Radius = %f\n",processRank,filterWidth);

    double oldMax = MathManager_Maximum(field);
    
    MathManager_CopyScalarField(field,temp1);
    MathManager_CopyScalarField(filter,temp2);
    MathManager_Convolve(temp1,temp2,field);
    FOREACH(x,y,z) field[x][y][z]*= field[x][y][z]>0.02;
    
    double newMax = MathManager_Maximum(field);
    if (newMax!=0) FOREACH(x,y,z) field[x][y][z]*= oldMax/newMax;
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// WallManager_FillSphere
// Fill sphere, making it smooth by super-sampling...
//******************************************************************************
void WallManager_FillSphere(SF field,double val, double* center, double radius)
{
    double c[3];
    int i,j,k; FOREACH_DIM(i) c[i] = center[i]*gridSize[i];
    double r2 = SQR(1.0*radius/L[0]*gridSize[0]);
    
    int n=3; //number of supersample points
    double step = 1.0/n;
    double norm = 1.0/CUBE((n+1));
    
    //For each pixel
     FOREACH(i,j,k)
     {
         double frac=0.0;
         
         //Supersample the edges
         int is,js,ks; 
         for(is=0; is<=n; is++) for(js=0; js<=n; js++) for(ks=0; ks<=n; ks++)
         {
             double x = i + is*step;
             double y = j + js*step;
             double z = k + ks*step;
             double d2 = SQR(x-c[0]) + SQR(y-c[1]) + SQR(z-c[2]);
             frac+= (d2<=r2);
         }
        frac*=norm;
         //and assign the average value to that pixel.
         if (frac>0)
         {
             if (val>0) field[i][j][k] = val*frac;
             else field[i][j][k]*=(1.0-frac);
         }
    }
}

//******************************************************************************
// WallManager_FillCylinder
//******************************************************************************
void WallManager_FillCylinder(SF field,double val, double* center, double radius)
{
    double c[3];
    int i,j,k; FOREACH_DIM(i) c[i] = center[i]*gridSize[i];
    double r2 = SQR(1.0*radius/L[0]*gridSize[0]);
    
    int n=3; //number of supersample points
    double step = 1.0/n;
    double norm = 1.0/CUBE((n+1));
    
    //For each pixel
    FOREACH(i,j,k)
    {
        double frac=0.0;
        
        //Supersample the edges
        int is,js,ks; 
        for(is=0; is<=n; is++) for(js=0; js<=n; js++) for(ks=0; ks<=n; ks++)
        {
            double x = i + is*step;
            double y = j + js*step;
            //double z = k + ks*step;
            double d2 = SQR(x-c[0]) + SQR(y-c[1]);
            frac+= (d2<=r2);
        }
            frac*=norm;
        //and assign the average value to that pixel.
        if (frac>0)
        {
            if (val>0) field[i][j][k] = val*frac;
            else field[i][j][k]*=(1.0-frac);
        }
    }
}

//******************************************************************************
// WallManager_FillTriangle
//******************************************************************************
void WallManager_FillTriangle(SF field,double val, double* center,double* width)
{
    int low[3],high[3];
    int i,j,k;
    FOREACH_DIM(i)
    {
        low[i]  = (center[i] - width[i]/2.0)*gridSize[i];
        high[i] = (center[i] + width[i]/2.0)*gridSize[i];
    }
    
    int n=3; //number of supersample points
    double step = 1.0/n;
    double norm = 1.0/CUBE((n+1));    
    double m = 0.5*width[0]/width[1];

    //For each pixel
    FOREACH(i,j,k)
    {
        double frac=0.0;
        double xlow = low[0] + m*(j-low[1]);
        double xhigh= high[0]- m*(j-low[1]);
        
        //Supersample the edges
        int is,js,ks; 
        for(is=0; is<=n; is++) for(js=0; js<=n; js++) for(ks=0; ks<=n; ks++)
        {
            double x = i + is*step;
            double y = j + js*step;
            double z = k + ks*step;
            frac+= ( (xlow<=x) && (x<=xhigh) && 
                     (low[1]<=y) && (y<=high[1]) && 
                     (low[2]<=z) && (z<=high[2]) );
        }
            frac*=norm;
        //and assign the average value to that pixel.
        if (frac>0)
        {
            if (val>0) field[i][j][k] = val*frac;
            else field[i][j][k]*=(1.0-frac);
        }
    }
    
}

//******************************************************************************
// WallManager_FillBox
//******************************************************************************
void WallManager_FillBox(SF field,double val, double* center,double* width)
{
    int low[3],high[3];
    int i; FOREACH_DIM(i)
    {
        low[i]  = (center[i] - width[i]/2.0)*gridSize[i];
        high[i] = (center[i] + width[i]/2.0)*gridSize[i];
    }
    
    int x,y,z; FOREACH(x,y,z)
    {
        if ( (low[0]<=x) && (x<=high[0]) && 
             (low[1]<=y) && (y<=high[1]) && 
             (low[2]<=z) && (z<=high[2]) ) 
            field[x][y][z]=val;
    }
}
