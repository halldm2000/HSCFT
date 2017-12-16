//******************************************************************************
//  SimulationManager.c
//  Hydrodynamic Self Consistent Field Theory Simulation
//  Created     David M. Hall   Mar 06 2005
//******************************************************************************
#include "Constants.h"
#include "FileManager.h"
#include "FileIO.h"
#include "FFTWManager.h"
#include "HydroModel.h"
#include "MathManager.h"
#include "MemoryManager.h"
#include "SimulationManager.h"
#include "StressModel.h"
#include "ThermoModel.h"
#include "WallManager.h"

//******************************************************************************
// SimulationManager_ComputeDimensionlessGroups
//******************************************************************************
void SimulationManager_ComputeDimensionlessGroups()
{
    SimulationManager_Trace("SimulationManager_ComputeDimensionlessGroups",2);

    int c,m;    

    // Calculate Characteristic Scales
    lC      = Rg[0];                FOREACH_MATERIAL(c) lC = MIN(lC,Rg[c]);
    tauD    = tauRep[0];            FOREACH_MATERIAL(c) tauD = MAX(tauD,tauRep[c]);
    zeta0C  = zeta0[0][0];          FORALL_MONOMERS(m,c) zeta0C = MIN(zeta0C, zeta0[c][m]);
    muC     = kT*NDENSITY;  // erg/cm3
    wC      = lC/tauD;      // cm/sec
    vC      = wC;           // cm/sec
    tC      = lC/vC;        // sec
    zetaC   = zeta0C*NDENSITY;

    // Compute Dimensionless Nunmbers
    De      = tauD/tC;
    Ca      = 5.0e-3;
    Gamma   = 1.0;
    
    // Display scaling factors
    printf("\t lC   = %.3e cm\n",lC);
    printf("\t tC   = %.3e sec\n",tC);
    printf("\t De   = %.3e \n",De);
    printf("\t Ca   = %.3e \n",Ca);
    printf("\t Gamma= %.3e \n",Gamma);
    printf("\t zetac= %.3e \n",zeta0C);
}

//******************************************************************************
// SimulationManager_InitializeSystemSize
//******************************************************************************
void SimulationManager_InitializeSystemSize()
{
    SimulationManager_Trace("SimulationManager_InitializeSystemSize",2);
    nx = gridSize[0]; ny = gridSize[1]; nz = gridSize[2];
    nmin=MIN(nx,ny);
    
    volume = L[0]*L[1]*L[2];
    norm = 1.0/(nx*ny*nz);
    
    GET_RUN_DIMENSION(dim);
}

//******************************************************************************
// SimulationManager_ComputeAverageN
// Compute the average degree of polymerization across all copolymers
//******************************************************************************
void SimulationManager_ComputeAverageN()
{
    meanN   = 0.0;
    minN    = N[0];
    
    int c; FOREACH_POLYMER(c)
    {
        printf("\t N[%d] = %f\n",c,N[c]);
        meanN+=N[c]*lambda[c];
        minN=MIN(minN,N[c]);
    }
    printf("\t meanN = %f\n",meanN);
}

//******************************************************************************
// SimulationManager_ComputeTotalNumBlocks
//******************************************************************************
void SimulationManager_ComputeTotalNumBlocks()
{
    int c; 
    numBlocks=0; FOREACH_MATERIAL(c) numBlocks+=numB[c];
    printf("Number of blocks in system = %d\n",numBlocks);
}

//******************************************************************************
// SimulationManager_ComputeStepsInPropagators
//******************************************************************************
void SimulationManager_ComputeStepsInPropagators()
{
    int b,c; 
    FOREACH_MATERIAL(c) ns[c]= stepsPerBlock*numB[c] + 1;//*N[c]/minN
    FORALL_BLOCKS(b,c) 
    {
        int steps = stepsPerBlock;//*N[c]/minN;
        deltaS[c][b]= (blockLength[c][b]/N[c])/steps;
    }
}

//******************************************************************************
// SimulationManager_MeasureMonomerVolumeFractions
//******************************************************************************
void SimulationManager_MeasureMonomerVolumeFractions()
{
    SimulationManager_Trace("SimulationManager_MeasureMonomerVolumeFractions",3);

    int m,c; FORALL_MONOMERS(m,c) 
    {
        f[c][m]=MathManager_Average(phi[c][m]);
        printf("\t f[%d][%d] \t=\t %f\n",c,m,f[c][m]);
    }
}

//******************************************************************************
// SimulationManager_MeasureCopolymerVolumeFractions
// Compute the total volume fraction (lambda) occupied by a copolymer by summing
// the volume frations of all monomers on that copolymer
//******************************************************************************
void SimulationManager_MeasureCopolymerVolumeFractions()
{    
    SimulationManager_Trace("SimulationManager_MeasureCopolymerVolumeFractions",3);

    int c; FOREACH_MATERIAL(c)
    {
        lambda[c]=0;
        int m; FOREACH_MONOMER(m,c) lambda[c]+=f[c][m];
        printf("\t lambda[%d]\t=\t %f\n",c,lambda[c]);
    }
}

//******************************************************************************
// SimulationManager_ComputeDerivedQuantities
//******************************************************************************
void SimulationManager_ComputeDerivedQuantities()
{   
    SimulationManager_Trace("SimulationManager_ComputeDerivedQuantities",2);
    int c; FOREACH_MATERIAL(c) tauRep[c]=tauNe[c]*pow(N[c]/Ne[c],3.4);

    SimulationManager_MeasureMonomerVolumeFractions();
    SimulationManager_MeasureCopolymerVolumeFractions();
    SimulationManager_ComputeAverageN();
    SimulationManager_ComputeTotalNumBlocks();
    SimulationManager_ComputeStepsInPropagators();
    FOREACH_MATERIAL(c) Rg[c]=MONOMER_SIZE*sqrt(N[c])/(2*dim);
    kT  = Kb*T;
}

//****************************************************************************
// SimulationManager_Trace
// Display trace information to track progress of the code
//****************************************************************************
void SimulationManager_Trace(char* message,int level)
{
    if (traceLevel>=level) 
    {
        printf("%d: ",processRank);
        int i; for (i=0; i<level-1; i++) printf(" ");
        printf("Trace %s\n",message);
    }
}

//****************************************************************************
// SimulationManager_DisplayMessage
// Display a message to standard out. 
//****************************************************************************
void SimulationManager_DisplayMessage(char* message)
{
    printf("%d: %s\n", processRank, message);
}

//******************************************************************************
// DisplayStatus
//******************************************************************************
void SimulationManager_DisplayStatus()
{
    SimulationManager_Trace("SimulationManager_DisplayStatus",1);
    
    if (timeStep % displayInterval == 0)
    {
        if (processRank == 0)
        {
            GET_TIME(timeEnd);
            printf("\n%s%d %s%.3f %s%.2e %s%d %s%.3f %s%.3f %s%.2e %s%f %s%.2e\n",   
                    "Step#",    timeStep,       "runTime=", runTime,
                    "deltaT=",  deltaT,   "#SCFT=",   scftStep, 
                    "maxPhiA=", maxPhi[0][0],   "minPhiA=", minPhi[0][0],
                    "scftErr=", scftErr,        "etime=",   timeEnd-timeStart ,
                    "scftErr1=",scftErr1);
            
            printf("%s%.1f %s%.2e %s%.2e %s%.6f %s%.6f \n",
                   "stepL=" ,rmsStepLength, 
                   "v="     ,rmsV,              "vTube=" ,rmsVTube,
                   "auxPhiA=", rmsPhiAux[0][0], "rmsPhiA=", rmsPhi[0][0]);
            
            printf("%s%.6f %s%.6f %s%.4f %s%.4f\n",
                   "meanPhiA=",meanPhi[0][0],   "meanAuxA=",meanPhiAux[0][0],
                   "rmsAlphaA=",rmsAlpha[0][0], "chi[0][0][0][1]=",chi[0][0][0][1]);
        }
    }
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// SimulationManager_InitializeData
//******************************************************************************
void SimulationManager_InitializeData()
{
    SimulationManager_Trace("SimulationManager_InitializeData",2);
    int c,m,x,y,z; FORALL_MONOMERS(m,c) FOREACH(x,y,z)
    {
        omega[c][m][x][y][z]=0;
    }    
    
    // If the run is a restart, override initial values with the last saved state
    if (restart) FileManager_ReadRestartFiles();
}

//******************************************************************************
// SimulationManager_PrintUsage
//******************************************************************************
void SimulationManager_PrintUsage()
{
    printf("Usage:\n");
    printf("hscft -x experiment_file\n");
    printf("hscft -r run_directory\n");
    printf("hscft -x experiment_file -r run_directory \n");
    printf("hscft -x experiment_file -t restart_time -r run_directory \n");
    exit(0);
}

//******************************************************************************
// SimulationManager_ReadCommandLineParams
//******************************************************************************
void SimulationManager_ReadCommandLineParams(int argc,char** argv)
{
    SimulationManager_Trace("SimulationManager_ReadCommandLineParams",2);
    
    restart=0;
    restartTime=0;
    restartTimestep=0;
    newExperiment=0;
    
    int i; for(i=1; i<argc; i++)
    {
        printf("%d: arg[%d] = %s\n",processRank,i,argv[i]);
        if(argv[i][0]=='-')
        {
            char c=argv[i][1];
            switch(c)
            {
                case 'r'://Restart finished (or broken) run
                    restart=1;
                    BROADCAST_DOUBLES(&restart,1);
                    i++;
                    strcpy(restartDir,argv[i]);
                    BROADCAST_STRING(restartDir);
                    printf("%d: runDir = %s\n",processRank,restartDir);
                    break;
                    
                case 'x'://Restart finished (or broken) run
                    newExperiment=1;
                    BROADCAST_DOUBLES(&newExperiment,1);
                    
                    i++;
                    strcpy(experimentFile,argv[i]);
                    BROADCAST_STRING(experimentFile);
                    printf("%d: experimentFile = %s\n",processRank,experimentFile);
                    break;
                    
                case 'd'://set distination dir
                    i++;
                    strcpy(dataDir,argv[i]);
                    BROADCAST_STRING(dataDir);
                    printf("%d: data_dir = %s\n",processRank,dataDir);
                    break;
                    
                case 't'://Choose target restart time
                    i++;
                    sscanf(argv[i],"%le",&restartTime);
                    BROADCAST_DOUBLES(&restartTime,1);
                    printf("%d: restart time = %f\n",processRank,restartTime);
                    break;
                    
                case 'h':
                    SimulationManager_PrintUsage();
                    break;
                    
                default:
                    printf("Unknown option %c\n",c);
                    SimulationManager_PrintUsage();
                    break;
            }
        }
        else SimulationManager_PrintUsage();
    }
    
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// SimulationManager_SeedRandomNumberGenerator
//******************************************************************************
void SimulationManager_SeedRandomNumberGenerator()
{
    SimulationManager_Trace("SimulationManager_SeedRandomNumberGenerator",2);
    unsigned seed = time(NULL) * (1+processRank);
    srand(seed);
    printf("%d: seed=%u\n",processRank,seed);
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// SimulationManager_SetStartTime
//******************************************************************************
void SimulationManager_SetStartTime()
{
    timeStep        = 0;
    runTime         = 0;
    lastRestartStep =-1;
    lastOutputStep  =-1;    
}

//******************************************************************************
// SimulationManager_InitializeModel
//******************************************************************************
void SimulationManager_InitializeModel(int argc,char** argv)
{
    SimulationManager_Trace("SimulationManager_InitializeModel",1);
        
    INITIALIZE_ALL_PROCESSES(argc,argv);
    GET_PROCESS_RANK(processRank);
    GET_NUM_PROCESSES(numProcesses);
    GET_TIME(timeStart);
    SimulationManager_SeedRandomNumberGenerator();
    
    SimulationManager_ReadCommandLineParams(argc,argv);
    FileManager_Initialize();
    SimulationManager_InitializeSystemSize();
    FFTWManager_Initialize();
    FFTWManager_MakeWaveTable();
    WallManager_Initialize();
    
    WallManager_ConstructSystemFromFile();
    SimulationManager_ComputeDerivedQuantities();
    SimulationManager_ComputeDimensionlessGroups();
    
    HydroModel_Initialize();
    ThermoModel_Initialize();
    StressModel_Initialize();
    SimulationManager_SetStartTime();
    SimulationManager_InitializeData();
}

//******************************************************************************
// SimulationManager_FinalizeModel
//******************************************************************************
void SimulationManager_FinalizeModel()
{
    SimulationManager_Trace("SimulationManager_FinalizeModel",1);
    FileManager_Finalize();
    HydroModel_Finalize();
    StressModel_Finalize();
    ThermoModel_Finalize();
    FFTWManager_Finalize();
    WallManager_Finalize();
    WAIT_FOR_ALL_PROCESSES;
    FINALIZE_ALL_PROCESSES;
}

//******************************************************************************
// SimulationManager_AntialiasFields
//******************************************************************************
void SimulationManager_AntialiasFields()
{
    int c,m,i,j;
    SimulationManager_Trace("SimulationManager_AntialiasFields",1);
    
    double order = 10.0;    
    FOREACH_DIM(i) MathManager_ApplyExponentialFilter(v[i],order);
    if (updateStresses) FORALL_MONOMERS(m,c) FOREACH_DIM(i) FOREACH_DIM(j)
    {
        MathManager_ApplyExponentialFilter(sigma[c][m][i][j],order);
    }
}

//******************************************************************************
// SimulationManager_MeasureSystemProperties
//******************************************************************************
void SimulationManager_MeasureSystemProperties()
{
    HydroModel_MeasureSystemProperties();
    StressModel_MeasureSystemProperties();
    ThermoModel_MeasureSystemProperties();
}

//******************************************************************************
// SimulationManager_IncrementTime
//******************************************************************************
void SimulationManager_IncrementTime()
{
    SimulationManager_Trace("SimulationManager_IncrementTime",1);
    timeStep++;
    runTime+=deltaT;
    if (runTime>=runEndTime) programStatus=STOP;
}

//******************************************************************************
// SimulationManager_AdvanceFields
//******************************************************************************
void SimulationManager_AdvanceFields()
{
    SimulationManager_Trace("SimulationManager_AdvanceFields",1);
    
    WallManager_ComputeLiquidFraction();
    if (updateThermo) ThermoModel_SeekLocalEquilibrium();
    if (updateThermo) ThermoModel_ComputeChemicalPotentials();
    WallManager_ComputeExternalForce();
    HydroModel_SeekSteadyState();
    HydroModel_AdvancePhi();
}

//******************************************************************************
// main
//******************************************************************************
int main(int argc, char** argv)
{   
    traceLevel= DEFAULT_TRACE_LEVEL; 
    SimulationManager_Trace("main",1);
    SimulationManager_InitializeModel(argc,argv);
    WallManager_ComputeLiquidFraction();
    HydroModel_RescalePolymerFields();
    FileManager_WriteDomainFile();
    FileManager_WriteDimensionalParameters();
    
    // Determined if it is okay to suppress mean field drift.
    suppressVDrift=0;
    numConstrained=0;
    int c; FOREACH_SOLID(c) if (tConstrained[c] || aConstrained[c]) numConstrained++;
    if (efMagnitude==0 && numConstrained==0) suppressVDrift=1;    
    
    programStatus = RUN;
    while(programStatus==RUN)
    {
        SimulationManager_DisplayStatus();
        FileManager_WriteFieldsToFile();
        FileManager_WriteRestartFiles();
        FileManager_AppendScalarsToFile();
        SimulationManager_IncrementTime();
        if (antialias) SimulationManager_AntialiasFields();
        SimulationManager_AdvanceFields();
    }
    SimulationManager_FinalizeModel();
    return 0;
}



