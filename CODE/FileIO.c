//******************************************************************************
//  FileIO.c
//******************************************************************************
#include "FFTWManager.h"
#include "FileIO.h"
#include "FileManager.h"
#include "HydroModel.h"
#include "MathManager.h"
#include "SimulationManager.h"
#include "StressModel.h"
#include "ThermoModel.h"
#include "WallManager.h"

#define APPENDS_DIRECTORY   "SCALARS"
#define STARTUP_DIRECTORY   "STARTUP_FILES"
#define WRITE_FLOAT_FORMAT  "%s\t%e\n"
#define DATE_FORMAT         "%b%d%Y"
#define TIME_FORMAT         "%H%M%S"
#define RUN_NAME_BASE       "run"
#define RUN_PARAMS_FILENAME "RunParameters.txt"
#define RUN_STATS_FILENAME  "RunStatstics.txt"
#define STRING_FORMAT       "%s\t%s\n"

//******************************************************************************
// FileIO_CopyToRunDir
//******************************************************************************
void FileIO_CopyToRunDir(char* descriptor)
{    
    if (processRank == 0)
    {
        char command[MAX_STRING];
        sprintf(command,"/bin/cp ./%s %s", descriptor, startupDest);
        printf("0: %s\n",command);
        system(command);
    }
}

//******************************************************************************
// FileIO_VerifyDir
//******************************************************************************
void FileIO_VerifyDir(char* path)
{
    if ((processRank == 0) && (stat(path,&FileIO_statBuffer)<0))
    {
        char command[MAX_STRING];
        sprintf(command,"/bin/mkdir %s", path);
        printf("%d: %s\n",processRank,command);
        system(command);
    }
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileIO_GetStartupFile
//******************************************************************************
void FileIO_GetStartupFile(char* fileName,char* pathToFile)
{
    sprintf(pathToFile,"%s/%s/%s", restartDir,STARTUP_DIRECTORY,fileName);
    printf("RESTART: %s = %s\n",fileName,pathToFile); 
}

//******************************************************************************
// FileIO_CreateRunDir
//******************************************************************************
void FileIO_CreateDirectories()
{
    FileIO_VerifyDir(dataDir);
    FileIO_VerifyDir(runDir);
    BROADCAST_STRING(runDir);
    
    if (processRank == 0) sprintf(appendDir,"%s/%s",runDir,APPENDS_DIRECTORY);
    FileIO_VerifyDir(appendDir);
    BROADCAST_STRING(appendDir);
    
    if (processRank == 0) sprintf(startupDest,"%s/%s",runDir,STARTUP_DIRECTORY);
    
    printf("startupDest=%s\n",startupDest);
    FileIO_VerifyDir(startupDest);
    BROADCAST_STRING(startupDest);

    sprintf(runStatsPath ,"%s/%s",startupDest,RUN_STATS_FILENAME );
    
    char command[MAX_STRING];
    sprintf(command,"echo -en \"\\e]0;%s\a\"",runName);
    system(command);
}

//******************************************************************************
// FileIO_AppendValue
//******************************************************************************
double FileIO_ReadLastScalar(char* fieldName)
{
    double result=-1;
    int incount=0;

    if (processRank == 0)
    {
        char path[MAX_STRING];
        sprintf(path,"%s/%s/%s.txt",restartDir,APPENDS_DIRECTORY,fieldName);
        printf("READ LAST SCALAR FROM %s\n",path);

        FILE* file = fopen(path,"r");
        do incount= fscanf(file,"%le", &result );       
        while(incount!=EOF);
        fclose(file);
        printf("result=%.3f\n",result);
    }
    WAIT_FOR_ALL_PROCESSES;
    return result;
}

//******************************************************************************
// FileIO_SeekScalarInFile
//******************************************************************************
int FileIO_SeekScalarInFile(char* fileName,double targetValue)
{
    int position=0;
    int incount=0;
    double result=0;
    
    if (processRank == 0)
    {
        char path[MAX_STRING];
        sprintf(path,"%s/%s/%s.txt",restartDir,APPENDS_DIRECTORY,fileName);
        printf("SEEK %f in %s\n",targetValue,path);
        
        FILE* file = fopen(path,"r");
        do 
        {
            incount+= fscanf(file,"%le", &result ); 
            position++;
        }
        while(incount!=EOF && result<targetValue);
        fclose(file);
        
        printf("result value=%.3f\n",result);
        printf("result position=%d\n",position);
    }
    WAIT_FOR_ALL_PROCESSES;
    return position;
}

//******************************************************************************
// FileIO_SeekScalarInFile
//******************************************************************************
double FileIO_ReadNthScalar(char* fileName,int index)
{
    int incount=0;
    int position=0;
    double result=-1;
    
    if (processRank == 0)
    {
        char path[MAX_STRING];
        sprintf(path,"%s/%s/%s.txt",restartDir,APPENDS_DIRECTORY,fileName);
        printf("READ scalar at position %d in %s\n",index,path);
        
        FILE* file = fopen(path,"r");
        do 
        {
            incount+= fscanf(file,"%le", &result ); 
            position++;
        }
        while(incount!=EOF && position<index);
        fclose(file);
        
        printf("result value=%.3f\n",result);
    }
    WAIT_FOR_ALL_PROCESSES;
    return result;
}

//******************************************************************************
// FileIO_AppendValue
//******************************************************************************
void FileIO_AppendScalar(char* fieldName,double value)
{
   
    if (processRank == 0)
    {
        char path[MAX_STRING];
        sprintf(path,"%s/%s.txt",appendDir,fieldName);
        FILE* file = fopen(path,"a");
        fprintf(file,"%e\n",value);
        fclose(file);
    }
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileIO_AppendVector
//******************************************************************************
void FileIO_AppendVector(char* fieldName,VEC vec,int dim)
{
    printf("Append vector %s\n",fieldName);
    if (processRank == 0)
    {
        char path[MAX_STRING];
        sprintf(path,"%s/%s.txt",appendDir,fieldName);
        FILE* file = fopen(path,"a");
        int i; for (i=0; i<dim; i++) fprintf(file,"%e ",vec[i]);
        fprintf(file,"\n");
        fclose(file);
    }
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileIO_WriteVector
//******************************************************************************
void FileIO_WriteVector(char* fieldName,VEC vec,int dim,int step)
{
    char dir[MAX_STRING];
    char path[MAX_STRING];

    //TODO: For MPI Case send the missing bits of data across
    if (processRank == 0)
    {
        sprintf(dir,"%s/%s",runDir,fieldName);
        FileIO_VerifyDir(dir);
        sprintf(path,"%s/%s_%d.txt",dir,fieldName,step);
        FILE* file = fopen(path,"w");
        
        int i; for (i=0; i<dim; i++) fprintf(file,"%e ",vec[i]);
        fprintf(file,"\n");
        fclose(file);
    }
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FilIO_AppendTensor
//******************************************************************************
void FileIO_AppendTensor(char* fieldName,double vals[MAX_C][MAX_M])
{
    int count=0;
    int m,c;
    FORALL_MONOMERS(m,c)
    {
        char name[MAX_STRING];
        sprintf(name,"%s%d",fieldName,count);
        FileIO_AppendScalar(name,vals[c][m]);
        count++;
    }
}

//******************************************************************************
// FileIO_GenerateUniqueRunName
//******************************************************************************
void FileIO_GenerateUniqueRunPath()
{
    char timeString[MAX_STRING];
    char dateString[MAX_STRING];
    char buffer[MAX_STRING];

    if (processRank==0)
    { 
        time_t calendarTime;
        time(&calendarTime);
        strftime(timeString, MAX_STRING, TIME_FORMAT, localtime(&calendarTime));
        strftime(dateString, MAX_STRING, DATE_FORMAT, localtime(&calendarTime));
        
        strcpy(buffer,description);
        sprintf(runName,"%s_%s_%s",dateString,timeString,description);  
    }
    BROADCAST_STRING(runName);
    
    if (processRank == 0) sprintf(runDir,"%s/%s",dataDir,runName);
}

//******************************************************************************
// FileIO_ReadVectorField
//******************************************************************************
void FileIO_ReadVectorField(char* field,VF vf,int numDims)
{
    char fieldName[MAX_STRING];
    int i=0; for (i=0; i<numDims; i++)
    {
        sprintf(fieldName,"%s%d",field,i);
        FileIO_ReadScalarField(fieldName,vf[i]);
    }
}

//******************************************************************************
// FileIO_ReadScalarField
//******************************************************************************
void FileIO_ReadScalarField(char* fieldName,SF sf)
{
    int x,y,z;
    char path[MAX_STRING];
    
    if (newExperiment)  
        sprintf(path,"%s/%s/%s_p%d_%d.txt",
                restartDir,fieldName,
                fieldName,processRank,outputStep);
    else
        sprintf(path,"%s/%s/%s_p%d_%d.txt",
                runDir,fieldName,
                fieldName,processRank,outputStep);
    
    if (resampleInput)
    {   
        printf("resampling file:%s\n",path);
        
        FILE* file = fopen(path,"r");
        if (file!=NULL)
        {
            for (x=0; x<oldGridSize[0]; x++) 
                for (y=0; y<oldGridSize[1]; y++) 
                    for (z=0; z<oldGridSize[2]; z++)
                        fscanf(file,"%le", &(resampleBuffer[x][y][z]) );
            fclose(file);
        }
        else printf("NULL file pointer.\n");
        MathManager_ResampleData(resampleBuffer,sf);
    }
    else 
    {
        printf("reading file:%s\n",path);
        FILE* file = fopen(path,"r");
        if (file!=NULL)  
        {
            FOREACH(x,y,z) fscanf(file,"%le", &(sf[x][y][z]) );
            fclose(file);
        }
        else printf("NULL file pointer.\n");
    }
}

//******************************************************************************
// FileIO_ReadDoubles
//******************************************************************************
void FileIO_ReadDoubles(char* buffer,double* result, FILE* file,int num)
{
    if (processRank ==0)
    {
        fscanf(file, "%s\t", buffer);
        printf("\t%s\t\t=", buffer);
        
        int i; for (i=0; i<num-1; i++)
        {
            fscanf(file, "%le\t", &(result[i]));
            printf("%.2e\t",result[i]);
        }
        fscanf(file, "%le\n", &(result[i]));
        printf("%.2e\n",result[i]);
    }
    BROADCAST_DOUBLES(result,num);
    BROADCAST_STRING(buffer);

}

//******************************************************************************
// FileIO_ReadInts
//******************************************************************************
void FileIO_ReadInts(char* buffer,int* result, FILE* file,int num)
{
    if (processRank ==0)
    {
        double temp;
        fscanf(file, "%s\t", buffer);
        printf("\t%s\t\t=", buffer);
        
        int i; for (i=0; i<num-1; i++)
        {
            fscanf(file, "%le\t", &temp);
            result[i] = (int)temp;
            printf("%d\t",result[i]);
        }
        fscanf(file, "%le\n", &temp);
        result[i] = (int)temp;
        printf("%d\n",result[i]);
    }
    BROADCAST_INTEGERS(result,num);
    BROADCAST_STRING(buffer);
}

//******************************************************************************
// FileIO_ReadString
//******************************************************************************
void FileIO_ReadString(FILE* file, char* result)
{
    if (processRank==0)
    {
        fscanf(file, STRING_FORMAT, buffer , result);
        printf("\t %s \t\t= %s\n", buffer, result);
    }
    BROADCAST_STRING(result);
}

//******************************************************************************
// FileIO_VerifyFile
//******************************************************************************
int FileIO_VerifyFile(char* fileName)
{
    int result = stat(fileName, &FileIO_statBuffer);
    return result;
}

//******************************************************************************
// FileIO_WriteParticleData
//******************************************************************************
void FileIO_WriteParticleData(char* name, T2 tensor[MAX_C],int step)
{
    char fileName[MAX_STRING];
    int c; FOREACH_SOLID(c)
    {
        sprintf(fileName,"%s%d",name,c);
        FileIO_WriteT2(fileName,tensor[c],numParticles[c],D,step);
    }
}

//******************************************************************************
// FileIO_ReadParticleData
//******************************************************************************
void FileIO_ReadParticleData(char* name, T2 tensor[MAX_C])
{
    char fileName[MAX_STRING];
    int c; FOREACH_SOLID(c)
    {
        sprintf(fileName,"%s%d",name,c);
        FileIO_ReadT2(fileName,tensor[c],numParticles[c],D);
    }
}

//******************************************************************************
// FileIO_WriteT2
//******************************************************************************
void FileIO_WriteT2(char* fieldName,T2 tensor,int dim1, int dim2, int step)
{
    char dir[MAX_STRING];
    char path[MAX_STRING];
    
    sprintf(dir,"%s/%s",runDir,fieldName);
    FileIO_VerifyDir(dir);
    sprintf(path,"%s/%s_p%d_%d.txt",dir,fieldName,processRank,step);
    FILE* file = fopen(path,"w");
    
    int i,j; for(i=0; i<dim1; i++) for (j=0; j<dim2; j++)
        fprintf(file,"%10e\n",tensor[i][j]);
    fclose(file);
    WAIT_FOR_ALL_PROCESSES;
    
}

//******************************************************************************
// FileIO_ReadT2
//******************************************************************************
void FileIO_ReadT2(char* fieldName,T2 tensor, int dim1, int dim2)
{
    int i,j;
    char path[MAX_STRING];
    
    sprintf(path,"%s/%s/%s_p%d_%d.txt",
            runDir,fieldName,
            fieldName,processRank,outputStep);
    
    printf("reading file:%s\n",path);
    FILE* file = fopen(path,"r");
    if (file!=NULL)  
    {
        for(i=0; i<dim1; i++) for (j=0;j<dim2; j++)
            fscanf(file,"%le", &(tensor[i][j]) );
        fclose(file);
    }
    else printf("NULL file pointer.\n");
}

//******************************************************************************
// FileIO_WriteScalarField
//******************************************************************************
void FileIO_WriteScalarField(char* fieldName,SF sf,int step)
{
    char dir[MAX_STRING];
    char path[MAX_STRING];
    
    sprintf(dir,"%s/%s",runDir,fieldName);
    FileIO_VerifyDir(dir);
    sprintf(path,"%s/%s_p%d_%d.txt",dir,fieldName,processRank,step);
    FILE* file = fopen(path,"w");

    int x,y,z; FOREACH(x,y,z) fprintf(file,"%10e\n",sf[x][y][z]);
    fclose(file);
    WAIT_FOR_ALL_PROCESSES;
}

//******************************************************************************
// FileIO_WriteVectorField
//******************************************************************************
void FileIO_WriteVectorField(char* field,VF vf,int numDims,int step)
{
    char fieldName[MAX_STRING];
    int i=0; for (i=0; i<numDims; i++)
    {
        sprintf(fieldName,"%s%d",field,i);
        FileIO_WriteScalarField(fieldName,vf[i],step);
    }
}

//******************************************************************************
// FileIO_WriteTensorField
//******************************************************************************
void FileIO_WriteTensorField(char* field,TF tf, int td2, int td1,int step)
{
    char fieldName[MAX_STRING];
    int i,j; for(i=0; i<td2; i++) for (j=0; j<td1; j++)
    {
        sprintf(fieldName,"%s%d%d",field,i,j);
        FileIO_WriteScalarField(fieldName,tf[i][j],step);
    }
}

//******************************************************************************
// FileIO_WriteTensorField
//******************************************************************************
void FileIO_ReadTensorField(char* field,TF tf, int td2, int td1)
{
    char fieldName[MAX_STRING];
    int i,j; for(i=0; i<td2; i++) for (j=0; j<td1; j++)
    {
        sprintf(fieldName,"%s%d%d",field,i,j);
        FileIO_ReadScalarField(fieldName,tf[i][j]);
    }
}

//******************************************************************************
// FileIO_WriteVariablSizeTF
//******************************************************************************
void FileIO_WriteVariableSizeTF(char* field,TF tf, int td2, int* td1,int step)
{
    int block =0;
    
    char fieldName[MAX_STRING];
    int i,j; for(i=0; i<td2; i++) for (j=0; j<td1[i]; j++)
    {
        sprintf(fieldName,"%s%d",field,block);
        FileIO_WriteScalarField(fieldName,tf[i][j],step);
        block++;
    }
}

//******************************************************************************
// FileIO_ReadVariableSizeTF
//******************************************************************************
void FileIO_ReadVariableSizeTF(char* field,TF tf, int td2, int* td1)
{
    int block =0;
    
    char fieldName[MAX_STRING];
    int i,j; for(i=0; i<td2; i++) for (j=0; j<td1[i]; j++)
    {
        sprintf(fieldName,"%s%d",field,block);
        FileIO_ReadScalarField(fieldName,tf[i][j]);
        block++;
    }
}

//******************************************************************************
// FileIO_ReadVariableSizeTF
//******************************************************************************
void FileIO_ReadVariableSizeTF4(char* field,TF4 tf, int td4, int* td3, int td2, int td1)
{
    int block =0;
    
    char fieldName[MAX_STRING];
    int i,j,k,l; for(i=0; i<td4; i++) for (j=0; j<td3[i]; j++)
    {
        int component=0;
        for (k=0; k<td2; k++) for (l=0; l<td1; l++)
        {
            sprintf(fieldName,"%s%d_%d",field,block,component);
            FileIO_ReadScalarField(fieldName,tf[i][j][k][l]);
            component++;
        }
        block++;
    }
}


//******************************************************************************
// FileIO_AppendParticleData
//******************************************************************************
void FileIO_AppendParticleData(char* name,T2 val[MAX_C])
{
    int c,j; 
    char fieldName[MAX_STRING];
    
    FOREACH_SOLID(c) FOREACH_PARTICLE(j,c) 
    {
        sprintf(fieldName,"%s%d",name,c);
        FileIO_AppendVector(fieldName,val[c][j],D);
    }
}

//******************************************************************************
// FileIO_AppendVariableSizeT3
//******************************************************************************
void FileIO_AppendVariableSizeT3(char* field,T3 t3, int td3, int* td2, int td1)
{
    int i,j; 
    char fieldName[MAX_STRING];
    int block=0;
    for(i=0; i<td3; i++) for (j=0; j<td2[i]; j++) 
    {
        sprintf(fieldName,"%s%d",field,block);
        FileIO_AppendVector(fieldName,t3[i][j],td1);
        block++;
    }
}

//******************************************************************************
// FileIO_WriteVariablSizeTF3
//******************************************************************************
void FileIO_WriteVariableSizeTF3(char* field,TF3 t3, int td3, int* td2, int td1,int step)
{
    int i,j,k; 
    char fieldName[MAX_STRING];
    int block=0;
    for(i=0; i<td3; i++) for (j=0; j<td2[i]; j++) 
    {
        for (k=0; k<td1; k++)
        {
            sprintf(fieldName,"%s%d_%d",field,block,k);
            FileIO_WriteScalarField(fieldName,t3[i][j][k],step);
        }
        block++;
    }
}

//******************************************************************************
// FileIO_WriteVariableSizeTF4
//******************************************************************************
void FileIO_WriteVariableSizeTF4(char* field,TF4 TF4, int td4, int* td3, int td2,int td1, int step)
{
    int i,j,k,l; 
    char fieldName[MAX_STRING];
    int block=0;
    for(i=0; i<td4; i++) for (j=0; j<td3[i]; j++) 
    {
        int component =0;
        for (k=0; k<td2; k++) for (l=0; l<td1; l++)
        {
            sprintf(fieldName,"%s%d_%d",field,block,component);
            FileIO_WriteScalarField(fieldName,TF4[i][j][k][l],step);
            component++;
        }
        block++;
    }
}

//******************************************************************************
// FileIO_WriteAsDouble
//******************************************************************************
void FileIO_WriteAsDouble(FILE* file, double value, char* label )
{
    fprintf(file, WRITE_FLOAT_FORMAT, label,(double)value);
}

//******************************************************************************
// FileIO_WriteString
//******************************************************************************
void FileIO_WriteString(FILE* file, char* value, char* label)
{
    fprintf(file, STRING_FORMAT, value , label);
}

