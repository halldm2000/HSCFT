//******************************************************************************
//  FileManager.c
//******************************************************************************
#include <string.h>
#include "Constants.h"
#include "FileIO.h"
#include "FileManager.h"
#include "FFTWManager.h"
#include "HydroModel.h"
#include "MemoryManager.h"
#include "SimulationManager.h"
#include "StressModel.h"
#include "ThermoModel.h"
#include "WallManager.h"

#define DOMAIN_SIZE_NAME        "DOMAIN_SIZE"
#define DIMENSIONAL_PARAMS_NAME "DIMENSIONAL_PARAMETERS"

//******************************************************************************
// FileManager_LookUpFloryParameters
//******************************************************************************
void FileManager_LookUpFloryParameters()
{
    SimulationManager_Trace("FileManager_LookUpFloryParameters",3);
    printf("Reading flory params from %s\n",floryParamsFile);
        
    int c1,c2,m1,m2; FORALL_MONOMERS(m1,c1) FORALL_MONOMERS(m2,c2)
    {
        if (strcmp(materialName[c1][m1],materialName[c2][m2])==0)
        {
            printf("\tchi[%d][%d][%d][%d]=0 \n",c1,m1,c2,m2);
            chi[c1][m1][c2][m2]=0.0;
        }
        else
        {
            char targetM1M2[MAX_STRING];
            char targetM2M1[MAX_STRING];
            sprintf(targetM1M2,"%s,%s",materialName[c1][m1],materialName[c2][m2]);
            sprintf(targetM2M1,"%s,%s",materialName[c2][m2],materialName[c1][m1]);
            
            FILE* file=0;
            int found=0;
            double value=0;
            if (processRank ==0) file = fopen(floryParamsFile,"r");
            while(!found)
            {
                FileIO_ReadDoubles(buffer,&value,file,1);
                if (strcmp(buffer,targetM1M2)==0) found=1;
                if (strcmp(buffer,targetM2M1)==0) found=1;
                if (found) 
                {
                    chi[c1][m1][c2][m2]=value;
                    printf("\tchi[%d][%d][%d][%d] = %f\n",
                           c1,m1,c2,m2,
                           chi[c1][m1][c2][m2]);
                }
            }
            if (processRank ==0) fclose(file);  
        }
    }
}

//******************************************************************************
// FileManager_WriteDimensionalParameters
//******************************************************************************
void FileManager_WriteDimensionalParameters()
{
    SimulationManager_Trace("FileManager_WriteDimensionalParameters",1);
    
    char path[MAX_STRING];
    sprintf(path,"%s/%s.txt",startupDest,DIMENSIONAL_PARAMS_NAME);
    
    FILE* file = fopen(path,"w");
    FileIO_WriteAsDouble(file,Ca        ,"CAPILLARY_NUMBER");
    FileIO_WriteAsDouble(file,vC        ,"CONVECTION_SPEED");
    FileIO_WriteAsDouble(file,tC        ,"CONVECTION_TIME");
    FileIO_WriteAsDouble(file,De        ,"DOBORAH_NUMBER");
    FileIO_WriteAsDouble(file,DENSITY   ,"DENSITY");
    FileIO_WriteAsDouble(file,wC        ,"DIFFUSION_SPEED");
    FileIO_WriteAsDouble(file,tauD      ,"DIFFUSION_TIME");
    FileIO_WriteAsDouble(file,Gc        ,"ELASTIC_MODULUS");
    FileIO_WriteAsDouble(file,Gamma     ,"FRICTION_FACTOR");
    FileIO_WriteAsDouble(file,lC        ,"LENGTH");
    FileIO_WriteAsDouble(file,NDENSITY  ,"NUMBER_DENSITY");
    FileIO_WriteAsDouble(file,zeta0C    ,"MONOMER_FRICTION");
    FileIO_WriteAsDouble(file,muC       ,"POTENTIAL_ENERGY_DENSITY");
    FileIO_WriteAsDouble(file,kT        ,"THERMAL_ENERGY");
    FileIO_WriteAsDouble(file,Gc*tauD   ,"VISCOSITY");
    fclose(file);
}

//******************************************************************************
// FileManager_Initialize
//******************************************************************************
void FileManager_Initialize()
{
    SimulationManager_Trace("FileManager_Initialize",2);
    resampleInput=0;
    printf("%d: restart=%d, newExperiment=%d\n",processRank,restart,newExperiment);
    
    FileManager_ReadExperimentFile();
    
    if (!restart) FileIO_GenerateUniqueRunPath();
    else strcpy(runDir,restartDir);
    FileIO_CreateDirectories();
    
    if (!restart) FileIO_CopyToRunDir("*");
    FileManager_ReadRunParameters(runParamsFile);
    FileManager_ReadOutputPrefs();
    
}

//******************************************************************************
// FileManager_Finalize
//******************************************************************************
void FileManager_Finalize()
{
    if (resampleInput) MemoryManager_FreeScalarField(resampleBuffer);
}

//******************************************************************************
// FileManager_ReadOutputPrefs
//******************************************************************************
void FileManager_ReadOutputPrefs()
{   
    FILE* file=0;
    if (processRank ==0) 
    {
        printf("0: Reading output prefs from %s\n",outputPrefsFile);
        file = fopen(outputPrefsFile,"r");
    }
    FileIO_ReadInts(buffer,&writeMu,file,1);
    FileIO_ReadInts(buffer,&writeOmegas,file,1);
    FileIO_ReadInts(buffer,&writeP,file,1);
    FileIO_ReadInts(buffer,&writePhi,file,1);
    FileIO_ReadInts(buffer,&writeProps,file,1);
    FileIO_ReadInts(buffer,&writeV,file,1);
    FileIO_ReadInts(buffer,&writeVDiv,file,1);
    FileIO_ReadInts(buffer,&writeVM,file,1);
    FileIO_ReadInts(buffer,&writeVSubsteps,file,1);
    FileIO_ReadInts(buffer,&writeVTube,file,1);
    FileIO_ReadInts(buffer,&writeW,file,1);
    FileIO_ReadInts(buffer,&writeForces,file,1);
    FileIO_ReadInts(buffer,&appendEnergy,file,1);
    FileIO_ReadInts(buffer,&appendEVF,file,1);
    FileIO_ReadInts(buffer,&appendForces,file,1);
    FileIO_ReadInts(buffer,&appendMu,file,1);
    FileIO_ReadInts(buffer,&appendOmegas,file,1);
    FileIO_ReadInts(buffer,&appendPhi,file,1);
    FileIO_ReadInts(buffer,&appendPhiGrads,file,1);
    FileIO_ReadInts(buffer,&appendResiduals,file,1);
    FileIO_ReadInts(buffer,&appendVels,file,1);
    FileIO_ReadInts(buffer,&traceLevel,file,1);
    if (processRank ==0) fclose(file);
    
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileManager_ReadExperimentFile
//******************************************************************************
void FileManager_ReadExperimentFile()
{   
    SimulationManager_DisplayMessage("Reading EXPERIMENT file.");
    FILE* file=0;
    
    if (processRank ==0) file = fopen(experimentFile,"r");
    FileIO_ReadString(file, description);
    FileIO_ReadString(file, runParamsFile);
    FileIO_ReadString(file, systemFile);
    FileIO_ReadString(file, floryParamsFile);
    FileIO_ReadString(file, outputPrefsFile);
    if (processRank ==0) fclose(file);
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileManager_ReadCompositionFile
//******************************************************************************
void FileManager_ReadCompositionFile(char* fileName)
{
    SimulationManager_Trace("FileManager_ReadCompositionFile",4);
    printf("COMPOSITION FILE=%s\n",fileName);
    
    FILE* file=0; if (processRank ==0) file = fopen(fileName,"r");
    FileIO_ReadInts(buffer,&copolymersInComp,file,1); 
    int c; for (c=numC; c<numC+copolymersInComp; c++) FileIO_ReadString(file, copolymerPath[c]);
    FileIO_ReadDoubles(buffer,&lambdaC[numC],file,copolymersInComp);
    if (processRank ==0) fclose(file); 
    
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileManager_ReadCopolymerFile
//******************************************************************************
void FileManager_ReadCopolymerFile(int c)
{
    SimulationManager_Trace("FileManager_ReadCopolymerFile",4);
    int m; 
    FILE* file=0;
    printf("MOLECULE FILE=%s\n",copolymerPath[c]);

    if (processRank ==0)file = fopen(copolymerPath[c],"r");
    FileIO_ReadInts(buffer,&numB[c],file,1);
    FileIO_ReadInts(buffer,&numM[c],file,1);
    FileIO_ReadDoubles(buffer,blockLength[c],file,numB[c]);
    FOREACH_MONOMER(m,c) FileIO_ReadDoubles(buffer,lineDensity[c][m],file,numB[c]);
    FileIO_ReadDoubles(buffer,&N[c],file,1);
    FileIO_ReadDoubles(buffer,&Ne[c],file,1);
    FileIO_ReadDoubles(buffer,&tauNe[c],file,1);
    FOREACH_MONOMER(m,c) FileIO_ReadString(file, materialPath[c][m]);
    if (processRank ==0) fclose(file);  
    
    if(N[c]>1) 
    {
        materialType[c]=POLYMER;
        printf("\tMaterial IS a polymer\n");
    }
    else 
    {
        materialType[c]=SIMPLE_LIQUID;
        printf("\tMaterial is NOT a polymer\n");

    }
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileManager_ReadMaterialFiles
//******************************************************************************
void FileManager_ReadMaterialFiles(int c)
{
    SimulationManager_DisplayMessage("Reading MATERIAL files.");
    FILE* file=0;

    int m; FOREACH_MONOMER(m,c)
    {
        if (processRank ==0) file = fopen(materialPath[c][m],"r");
        printf("Reading %s\n",materialPath[c][m]);
        FileIO_ReadString(file, materialName[c][m]);
        FileIO_ReadDoubles(buffer,&G0[c][m],file,1);
        FileIO_ReadDoubles(buffer,&K0[c][m],file,1);
        FileIO_ReadDoubles(buffer,&zeta0[c][m],file,1);
        FileIO_ReadDoubles(buffer,&density[c][m],file,1);
        FileIO_ReadDoubles(buffer,&mViscosity[c][m],file,1);
        if (processRank ==0) fclose(file);  
    }
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileManager_ReadMotionDescriptor
//******************************************************************************
void FileManager_ReadMotionDescriptor(int c)
{
    FILE* file=0;
    if (processRank ==0) 
    {
        file = fopen(motionDescriptorFile,"r");
        printf("0: Reading motion descriptor %s\n",motionDescriptorFile);
    }
    FileIO_ReadDoubles(buffer,&aConstrained[c],file,1);
    FileIO_ReadDoubles(buffer,&angSpeedMax[c],file,1);
    FileIO_ReadDoubles(buffer,&angFreq[c],file,1);
    FileIO_ReadDoubles(buffer,&angPhase[c],file,1);
    FileIO_ReadDoubles(buffer,&tConstrained[c],file,1);
    FileIO_ReadDoubles(buffer,transVMax[c],file,D);
    FileIO_ReadDoubles(buffer,transFreq[c],file,D);
    FileIO_ReadDoubles(buffer,transPhase[c],file,D);
    if (processRank ==0) fclose(file);
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileManager_WriteDomainFile
//******************************************************************************
void FileManager_WriteDomainFile()
{
    SimulationManager_Trace("FileManager_WriteDomainFile",1);
    char path[MAX_STRING];
    sprintf(path,"%s/%s_p%d.txt",startupDest,DOMAIN_SIZE_NAME,processRank);
    printf("Domain file path:%s\n",path);
    
    FILE* file = fopen(path,"w");
    FileIO_WriteAsDouble(file,nxLocal,"nxLocal");
    FileIO_WriteAsDouble(file,ny,"nyLocal");
    FileIO_WriteAsDouble(file,nz,"nzLocal");
    FileIO_WriteAsDouble(file,startX,"startX");
    FileIO_WriteAsDouble(file,startY,"startY");
    FileIO_WriteAsDouble(file,startZ,"startZ");
    FileIO_WriteAsDouble(file,nyTrans,"nyTrans");
    FileIO_WriteAsDouble(file,totalSize ,"totalSize");
    FileIO_WriteAsDouble(file,numProcesses ,"numProcesses");
    FileIO_WriteAsDouble(file,numBlocks,"numBlocks");
    FileIO_WriteAsDouble(file,numLiquids,"numLiquids");
    FileIO_WriteAsDouble(file,numSolids,"numSolids");
    FileIO_WriteAsDouble(file,numC,"numC");
    fclose(file);
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileManger_ReadRunParameters
// Read parameters from file on master node and broadcast to other nodes
//******************************************************************************
void FileManager_ReadRunParameters(char* path)
{
    FILE* file=0;
    if (processRank ==0)
    {
        printf("0: Reading run params from %s\n",path);
        file=fopen(path,"r");
    }
    FileIO_ReadInts     (buffer,&antialias,file,1);
    FileIO_ReadInts     (buffer,&appendInterval,file,1);
    FileIO_ReadDoubles  (buffer,&CahnNumber,file,1);
    FileIO_ReadInts     (buffer,&displayInterval,file,1);
    FileIO_ReadDoubles  (buffer,&efMagnitude,file,1);
    FileIO_ReadDoubles  (buffer,&filterWidth,file,1);
    FileIO_ReadInts     (buffer,gridSize,file,3);
    FileIO_ReadDoubles  (buffer,&epsilon,file,1);
    FileIO_ReadDoubles  (buffer,L,file,3);
    FileIO_ReadDoubles  (buffer,&maxSCFTResidual,file,1);
    FileIO_ReadInts     (buffer,&maxSCFTSteps,file,1);
    FileIO_ReadDoubles  (buffer,&maxThermoStep,file,1);
    FileIO_ReadDoubles  (buffer,&noise,file,1);
    FileIO_ReadInts     (buffer,&stepsPerBlock,file,1);
    FileIO_ReadDoubles  (buffer,&runEndTime,file,1);
    FileIO_ReadInts     (buffer,&scftInterval,file,1);
    FileIO_ReadDoubles  (buffer,&deltaT,file,1);
    FileIO_ReadDoubles  (buffer,&deltaT_v,file,1);
    FileIO_ReadDoubles  (buffer,&RE,file,1);
    FileIO_ReadInts     (buffer,&updateStresses,file,1);
    FileIO_ReadInts     (buffer,&updateThermo,file,1);
    FileIO_ReadInts     (buffer,&useUpperConvected,file,1);
    FileIO_ReadInts     (buffer,&vDisplayInterval,file,1);
    FileIO_ReadDoubles  (buffer,&vThreshold,file,1);
    FileIO_ReadDoubles  (buffer,&writeInterval,file,1);
    FileIO_ReadDoubles  (buffer,&restartInterval,file,1);
    if (processRank ==0) fclose(file);
    
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileManager_WriteRestartFiles
//******************************************************************************
void FileManager_WriteRestartFiles()
{
    int restartStep = runTime/restartInterval;
    if (restartStep!=lastRestartStep)
    { 
        lastRestartStep = restartStep;
        FileIO_AppendScalar("restart_runTime",runTime);
        FileIO_AppendScalar("restart_timeStep",timeStep);
        FileIO_AppendScalar("restart_outputStep",outputStep);
        
        FileIO_WriteVariableSizeTF("omega",omega,numC,numM,outputStep);
        FileIO_WriteVariableSizeTF("phi",phi,numC,numM,outputStep);
        FileIO_WriteVectorField("v" ,v,D,outputStep);
        FileIO_WriteScalarField("p" ,P,outputStep);
        
        if (updateStresses) FileIO_WriteVariableSizeTF4("sigma",sigma,numC,numM,D,D,outputStep);
        
        int c; FOREACH_SOLID(c)
        {
            FileIO_WriteParticleData("rCM",rCM,outputStep);
            FileIO_WriteParticleData("vCM",vCM,outputStep);
            FileIO_WriteParticleData("aCM",aCM,outputStep);
            FileIO_WriteParticleData("avCM",avCM,outputStep);
        }
    }
}

//******************************************************************************
// FileManager_ReadRestartFiles
//******************************************************************************
void FileManager_ReadRestartFiles()
{
    
    if (restartTime==0)
    {
        timeStep    = FileIO_ReadLastScalar("restart_timeStep");
        runTime     = FileIO_ReadLastScalar("restart_runTime");
        outputStep  = FileIO_ReadLastScalar("restart_outputStep");
    }
    else
    {
        int index = FileIO_SeekScalarInFile("restart_runTime",restartTime);

        timeStep    = FileIO_ReadNthScalar("restart_timeStep",index);
        runTime     = FileIO_ReadNthScalar("restart_runTime",index);
        outputStep  = FileIO_ReadNthScalar("restart_outputStep",index);
        restartTimestep=timeStep;
    }
    
    FileIO_ReadVariableSizeTF("omega",omega,numC,numM);
    FileIO_ReadVariableSizeTF("phi",phi,numC,numM);
    FileIO_ReadVectorField("v" ,v,D);
    FileIO_ReadScalarField("p" ,P);
    if (updateStresses) FileIO_ReadVariableSizeTF4("sigma",sigma,numC,numM,D,D);
    
    int c; FOREACH_SOLID(c)
    {
        FileIO_ReadParticleData("rCM",rCM);
        FileIO_ReadParticleData("vCM",vCM);
        FileIO_ReadParticleData("aCM",aCM);
        FileIO_ReadParticleData("avCM",avCM);
    }
}

//******************************************************************************
// FileManager_WriteFieldsToFile
//******************************************************************************
void FileManager_WriteFieldsToFile()
{
    SimulationManager_Trace("FileManager_WriteDataToFile",1);
    
    outputStep = runTime/writeInterval;
    if (outputStep!=lastOutputStep)
    { 
        lastOutputStep = outputStep;
        FileIO_AppendVector("L", L, 3);        
        
        FileIO_WriteScalarField("vMag" ,vMag,outputStep);
        FileIO_WriteScalarField("solid" ,solid,outputStep);
      
        if (outputStep==0) 
        {
         //   FileIO_WriteVectorField("particle" ,particle,MAX_C,outputStep);
         //   FileIO_WriteScalarField("potential" ,potential,outputStep);
         //   FileIO_WriteVectorField("particlePE" ,particlePE,MAX_C,outputStep);
        }
        
        if(writeForces)
        {
            FileIO_WriteVectorField("wallForces", wallForces,D,outputStep);
            FileIO_WriteVectorField("elasticForce", elasticForceSum,D,outputStep);
            FileIO_WriteVectorField("gradPi",gradPi,D,outputStep);
            FileIO_WriteVectorField("gradP",gradP,D,outputStep);
            FileIO_WriteVariableSizeTF3("fPsi",fPsi,numC,numM,D,outputStep);
        }
        
        if(writeMu) FileIO_WriteVariableSizeTF("mu",mu,numC,numM,outputStep);
        
        if(writeOmegas) FileIO_WriteVariableSizeTF("omega",omega,numC,numM,outputStep);

        if(writeP)  FileIO_WriteScalarField("p" ,P,outputStep);
                
        if(writePhi) FileIO_WriteVariableSizeTF("phi",phi,numC,numM,outputStep);
        
        if(writeProps)
        {
            FileIO_WriteVariableSizeTF("q",       q, numC,ns,outputStep);
            FileIO_WriteVariableSizeTF("qdagger", qdagger,numC,ns,outputStep);
        }
        
        if (writeVDiv) FileIO_WriteScalarField("vDiv" ,vDiv,outputStep);
        
        if(writeV)  
        {
            FileIO_WriteVectorField("v" ,v,D,outputStep);
            //FileIO_WriteVectorField("vWall" ,vWall,D,outputStep);
        }
        
        if(writeVTube) FileIO_WriteVectorField("vTube",vTube,D,outputStep);
        if(writeVM) FileIO_WriteVariableSizeTF3("vPhi",vPhi,numC,numM,D,outputStep);
        if(writeW) FileIO_WriteVariableSizeTF3("w",w,numC,numM,D,outputStep);
    }
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileManager_AppendScalarsToFile
//******************************************************************************
void FileManager_AppendScalarsToFile()
{
    SimulationManager_Trace("FileManager_AppendScalarsToFile",1);

    if ((appendInterval>0) && (timeStep % appendInterval ==0))
    {
        SimulationManager_MeasureSystemProperties();
                
        //    FileIO_AppendParticleData("netForce",netForce,numC,numParticles,D);
        //    FileIO_AppendParticleData("netTorque",netTorque,numC,numParticles,D);
        //    FileIO_AppendParticleData("rCM",rCM,numC,numParticles,D);
        //    FileIO_AppendParticleData("vCM",vCM,numC,numParticles,D);
        //    FileIO_AppendParticleData("aCM",aCM,numC,numParticles,D);
        //    FileIO_AppendParticleData("avCM",avCM,numC,numParticles,D);

        FileIO_AppendScalar("fSolid",fSolid);
        FileIO_AppendScalar("runTime",runTime);

        if (appendEVF)
        {
            FileIO_AppendTensor("evf",evf);
            FileIO_AppendTensor("avf",avf);
        }
        
        if(appendResiduals)
        {
            FileIO_AppendScalar("scftErr",scftErr);
            FileIO_AppendTensor("omegaErr",rmsOmegaErr);
            FileIO_AppendScalar("rmsVDiv",rmsVDiv);
            FileIO_AppendScalar("rmsVErr",rmsVErr);
            FileIO_AppendScalar("rmsdVdL",rmsdVdL);
        }
        
        if(appendPhi)
        {
            FileIO_AppendTensor("maxPhi",maxPhi);
            FileIO_AppendTensor("meanPhi",meanPhi);
            FileIO_AppendTensor("minPhi",minPhi);
        }
        
        if (appendPhiGrads) FileIO_AppendTensor("rmsGradPhi",rmsGradPhi);

        if (appendEnergy) {}
        
        if (appendMu)FileIO_AppendTensor("rmsMu",rmsMu);

        if (appendOmegas)
        {
            FileIO_AppendTensor("meanOmega",meanOmega);
            FileIO_AppendTensor("rmsOmega",rmsOmega);
        }
        
        if (appendVels)
        {
            FileIO_AppendScalar("rmsV",rmsV);
            FileIO_AppendScalar("meanV",meanV);
            FileIO_AppendTensor("rmsVM",rmsVM);
            FileIO_AppendScalar("rmsVTube",rmsVTube);
            FileIO_AppendTensor("rmsW",rmsW);
        }
        
        if (appendForces)
        {
            FileIO_AppendScalar("elasticForceMax",maxElasticForce);
            FileIO_AppendScalar("osmoticForceMax",maxOsmoticForce);
            FileIO_AppendScalar("viscousForceMax",maxViscousForce);
            FileIO_AppendScalar("wallForceMax",maxWallForce);
        }
    }
    WAIT_FOR_ALL_PROCESSES;
}
