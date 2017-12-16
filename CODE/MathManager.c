//******************************************************************************
//  MathManager.c
//  Created by halldm on Fri Jul 23 2004.
//******************************************************************************
#include <float.h>
#include "MathManager.h"
#include "FFTWManager.h"
#include "SimulationManager.h"
#include "StressModel.h"
#include "HydroModel.h"

//******************************************************************************
// MathManager_MinDoubleInVector
//******************************************************************************
double MathManager_MinDoubleInVector(VEC vector,int dim)
{
    double min=vector[0];
    int i; for(i=0; i<dim; i++) min=MIN(min,vector[i]);
    return min;
}

//******************************************************************************
// MathManager_MinIntegerInVector
//******************************************************************************
int MathManager_MinIntegerInVector(int* vector,int dim)
{
    int min=vector[0];
    int i; for(i=0; i<dim; i++) min=MIN(min,vector[i]);
    return min;
}

//******************************************************************************
// MathManager_Divergence
//******************************************************************************
void MathManager_Divergence(VF vf,SF sf)
{
    int freqs[D],index,i;
    double k[D];
    
    MathManager_ClearScalarField(sf);
    fftw_complex* sfTilda = (fftw_complex*)sf[0][0];
    
    FOREACH_DIM(i)
    {
        MathManager_CopyScalarField(vf[i],work2);
        fftw_complex* vfTilda = (fftw_complex*)work2[0][0];
        FOURIER_TRANSFORM(work2,work1);
        
        FOREACHFREQ(freqs)
        {
            index = FINDEX(freqs);
            GETWAVENUMS(freqs,k);
        
            sfTilda[index].re += -k[i]*vfTilda[index].im*norm;
            sfTilda[index].im += +k[i]*vfTilda[index].re*norm;
        }
    }
    INVERSE_TRANSFORM(sf,work1);
}

//******************************************************************************
// MathManager_DivTensor
//******************************************************************************
void MathManager_DivTensorField(TF tf, VF vf)
{
    int r,c;
    FOREACH_DIM(c)
    {
        FOREACH_DIM(r) MathManager_CopyScalarField(tf[r][c],workVector[r]);
        MathManager_Divergence(workVector,vf[c]);
    }
}

//******************************************************************************
// MathManager_VectorDotVector
//******************************************************************************
void MathManager_VectorDotVector(VF vf1, VF vf2, SF vf1_DOT_vf2)
{
    int x,y,z,i;
    MathManager_ClearScalarField(vf1_DOT_vf2);
    FOREACH_VINDEX(i,x,y,z)
    {
        vf1_DOT_vf2[x][y][z]+= vf1[i][x][y][z]*vf2[i][x][y][z];
    }
}

//******************************************************************************
// MathManager_Gradient
// ...with linear strainYX
//******************************************************************************
void MathManager_Gradient(SF sf,VF vf)
{
    int freqs[D],index,i;
    double k[D];
        
    MathManager_CopyScalarField(sf,work2);
    fftw_complex* sfTilda = (fftw_complex*)work2[0][0];    
    
    fftw_complex* deriv[D];
    FOREACH_DIM(i) deriv[i] = (fftw_complex*)vf[i][0][0];
        
    FOURIER_TRANSFORM(work2,work1);
    FOREACHFREQ(freqs)
    {
        index = FINDEX(freqs);
        GETWAVENUMS(freqs,k);        
        
        FOREACH_DIM(i)
        {
            double temp=k[i]*norm;
            deriv[i][index].re = -temp*sfTilda[index].im;
            deriv[i][index].im = +temp*sfTilda[index].re;
        }
    }
    FOREACH_DIM(i) INVERSE_TRANSFORM(vf[i],work1);
}

//******************************************************************************
// MathManager_Laplacian
//******************************************************************************
void MathManager_Laplacian(SF sf,SF laplacian)
{
    int freqs[D],index;
    double k[D],k2;
    
    MathManager_CopyScalarField(sf,work2);
    
    fftw_complex* sfTilda = (fftw_complex*)work2[0][0];
    fftw_complex* lTilda  = (fftw_complex*)laplacian[0][0];
    
    FOURIER_TRANSFORM(work2,work1);    
    FOREACHFREQ(freqs)
    {
        index = FINDEX(freqs);
        GETWAVENUMS(freqs,k); 	 
        k2 = SQUARE_MAG(k);
        double laplacian = -( SQR(k[0]) + SQR(k[1]) + SQR(k[2]) );
        lTilda[index].re = laplacian*sfTilda[index].re*norm;
        lTilda[index].im = laplacian*sfTilda[index].im*norm;
    }
    INVERSE_TRANSFORM(laplacian,work1);
}

//******************************************************************************
// MathManager_MeasureLengthScale
//******************************************************************************
double  MathManager_MeasureMeanLengthScale(SF sf)
{
    int freqs[D],index;
    double k[D],k2,kmag,s;
    MathManager_CopyScalarField(sf,work2);
    fftw_complex* sfTilda = (fftw_complex*)work2[0][0];    

    double kSum=0;
    double sSum=0;
    FOURIER_TRANSFORM(work2,work1);    
    FOREACHFREQ(freqs)
    {
        index = FINDEX(freqs);
        GETWAVENUMS(freqs,k); 
        k2 = SQUARE_MAG(k);
        kmag = sqrt(k2);
        if (kmag!=0)
        {
            s = SQR(sfTilda[index].re) + SQR(sfTilda[index].im);
            kSum+=s*kmag/TWO_PI;
            sSum+=s;
        }
    }

    double kTotal=0 ;FIND_GLOBAL_SUM(kSum,kTotal);
    double sTotal=0; FIND_GLOBAL_SUM(sSum,sTotal);
    double lengthScale = sTotal/kTotal;
    return lengthScale;
}

//******************************************************************************
// MathManager_SolvePoissonEquation
//******************************************************************************
void MathManager_SolvePoissonEquation(SF sf,SF result)
{
    int freqs[D],index;
    double k[D],k2;
    
    MathManager_CopyScalarField(sf,work2);
    
    fftw_complex* sfTilda = (fftw_complex*)work2[0][0];
    fftw_complex* rTilda  = (fftw_complex*)result[0][0];
    
    FOURIER_TRANSFORM(work2,work1);
    FOREACHFREQ(freqs)
    {
        index = FINDEX(freqs);
        GETWAVENUMS(freqs,k); 	 
        k2 = SQUARE_MAG(k);
        double laplacian = -(SQR(k[0]) + SQR(k[1]) + SQR(k[2]));
        rTilda[index].re = sfTilda[index].re*norm/laplacian;
        rTilda[index].im = sfTilda[index].im*norm/laplacian;
    }
    rTilda[0].re=0.0;
    rTilda[0].im=0.0;
    INVERSE_TRANSFORM(result,work1);
}

//******************************************************************************
// MathManager_DownsampleData
//******************************************************************************
void MathManager_DownsampleData(SF oldField,SF newField)
{
    int rx = oldGridSize[0]/gridSize[0];
    int ry = oldGridSize[1]/gridSize[1];
    int rz = oldGridSize[2]/gridSize[2];
    
    printf("rx=%d, ry=%d rz=%d ",rx,ry,rz);
    
    int x,y,z; FOREACH(x,y,z)
    {
        int count=0;
        double value = 0.0;
        int i,j,k; for (i=-rx+1; i<rx; i++) for (j=-ry+1; j<ry; j++) for (k=-rz+1; k<rz; k++)
        {
            int x2 = (x*rx + i + oldGridSize[0])%oldGridSize[0];
            int y2 = (y*ry + j + oldGridSize[1])%oldGridSize[1];
            int z2 = (z*rz + k + oldGridSize[2])%oldGridSize[2];
            value+= oldField[x2][y2][z2];
            count++;
        }
            value/=count;
        newField[x][y][z]= value;
    }
}

//******************************************************************************
// MathManager_UpsampleData
//******************************************************************************
void MathManager_UpsampleData(SF oldField,SF newField)
{
    int rx = gridSize[0]/oldGridSize[0];
    int ry = gridSize[1]/oldGridSize[1];
    int rz = gridSize[2]/oldGridSize[2];
    
    printf("rx=%d, ry=%d rz=%d\n",rx,ry,rz);
    
    /*int i,j,k; FOREACH(i,j,k)
    {
        double value = 0.0;
        
        //use trilinear interpolation to find the value at a higher resolution
        int x = (i/rx); double fx = i/rx-x;
        int y = (j/ry); double fy = j/ry-y;
        int z = (k/rz); double fz = k/rz-z;
        
        //interp along z;
        double a[2][2];
        int h,v; 
        for (h=0; h<2; h++) for (v=0; v<2; v++)
        {
            int x0 = (x+h)%oldGridSize[0];
            int y0 = (y+v)%oldGridSize[1];
            int z0 = (z)%oldGridSize[2];
            int z1 = (z+1)%oldGridSize[2];

            a[h][v] = oldField[x0][y0][z0]*(1-fz) + oldField[x0][y0][z1]*fz;
        }
        //intep along y;
        double b[2];
        b[0] = a[0][0]*(1-fy) + a[0][1]*fy;
        b[1] = a[1][0]*(1-fy) + a[1][1]*fy;

        //interp along x;
        value = b[0]*(1-fx) + b[1]*fx;

        newField[i][j][k]= value;
    }*/
    
    //interpolate to find the value at a higher resolution
    int i,j,k; FOREACH(i,j,k)
    {        
        int x = (i/rx); double fx = i/rx-x;
        int y = (j/ry); double fy = j/ry-y;
        int z = (k/rz); double fz = k/rz-z;
                
        //linear interp along z
        double data2d[4][4];
        int h,v; for (h=0; h<4; h++) for (v=0; v<4; v++)
        {
            int x0 = (x+h-1+oldGridSize[0])%oldGridSize[0];
            int y0 = (y+v-1+oldGridSize[1])%oldGridSize[1];
            int z0 = (z)%oldGridSize[2];
            int z1 = (z+1)%oldGridSize[2];
            
            data2d[h][v] = oldField[x0][y0][z0]*(1-fz) + oldField[x0][y0][z1]*fz;
        }
        
        //cubic intep along each y;
        double data1d[4];
        double fy2=fy*fy;
        double fy3=fy2*fy;

        for (h=0; h<4; h++)
        {
            double f0 = data2d[h][0];
            double f1 = data2d[h][1];
            double f2 = data2d[h][2];
            double f3 = data2d[h][3];

            double a0 = f3 - f2 - f0 +f1;
            double a1 = f0 - f1 - a0;
            double a2 = f2 - f0;
            double a3 = f1;
            
            data1d[h]=a0*fy3 + a1*fy2 + a2*fy + a3;
        }
        
        //interp along x;
        double fx2=fx*fx;
        double fx3=fx2*fx;

        double f0 = data1d[0];
        double f1 = data1d[1];
        double f2 = data1d[2];
        double f3 = data1d[3];
        
        double a0 = f3 - f2 - f0 + f1;
        double a1 = f0 - f1 - a0;
        double a2 = f2 - f0;
        double a3 = f1;
        
        newField[i][j][k]=a0*fx3+a1*fx2+a2*fx+a3;
    }
}

//******************************************************************************
// MathManager_Resample
//******************************************************************************
void MathManager_ResampleData(SF oldField,SF newField)
{
    double rx = oldGridSize[0]/gridSize[0];
    double ry = oldGridSize[1]/gridSize[1];
    double rz = oldGridSize[2]/gridSize[2];

    printf("rx=%f, ry=%f rz=%f ",rx,ry,rz);

    if (rx>1) MathManager_DownsampleData(oldField,newField);
    if (rx<1) MathManager_UpsampleData(oldField,newField);
}

//******************************************************************************
// MathManager_Maximum
//******************************************************************************
double MathManager_Maximum(SF sf)
{
    double localMax = sf[0][0][0];
    int x,y,z; FOREACH(x,y,z) localMax=MAX(fabs(localMax),sf[x][y][z]);
    
    double max=localMax;
    FIND_GLOBAL_MAX(localMax,max);
    
    return max;
}

//******************************************************************************
// MathManager_Maximum
//******************************************************************************
double MathManager_Minimum(SF sf)
{
    double localMin = sf[0][0][0];
    int x,y,z; FOREACH(x,y,z) localMin=MIN(localMin,sf[x][y][z]);
    
    double min=localMin;
    FIND_GLOBAL_MIN(localMin,min);    
    return min;
}


//******************************************************************************
// MathManager_SumOfSquares
//******************************************************************************
double  MathManager_SumOfSquares(SF sf)
{
    double localSum=0.0;
    int x,y,z; FOREACH(x,y,z) localSum+=SQR(sf[x][y][z]);

    double sum=0.0;
    FIND_GLOBAL_SUM(localSum,sum);
    return sum;
}

//******************************************************************************
// MathManager_Average
// Compute the Average value across all processors
//******************************************************************************
double MathManager_Average(SF sf)
{
    double localSum = 0.0;
    int x,y,z; FOREACH(x,y,z) localSum+=sf[x][y][z];
    
    double sum=0.0;
    FIND_GLOBAL_SUM(localSum,sum);
    return (sum*norm);
}

//******************************************************************************
// MathManager_CrossProduct
// Compute the cross product of two 3D vectors
//******************************************************************************
void MathManager_CrossProduct(V3 a,V3 b,V3 result)
{
    result[0]=a[1]*b[2]-a[2]*b[1];
    result[1]=a[2]*b[0]-a[0]*b[2];
    result[2]=a[0]*b[1]-a[1]*b[0];
}

//******************************************************************************
// MathManager_Sum
// Compute the Average value across all processors
//******************************************************************************
double MathManager_Sum(SF sf)
{
    double localSum = 0.0;
    int x,y,z; FOREACH(x,y,z) localSum+=sf[x][y][z];
    
    double sum=0.0;
    FIND_GLOBAL_SUM(localSum,sum);
    return sum;
}

//******************************************************************************
// MathManager_LocalMagnitude
//******************************************************************************
double MathManager_LocalMagnitude(VF vf,int x,int y,int z)
{
    double sum = 0;
    int i; FOREACH_DIM(i) sum+=SQR(vf[i][x][y][z]);
    return sqrt(sum);
}

//******************************************************************************
// MathManager_MaximumMagnitude
//******************************************************************************
double MathManager_MaximumMagnitude(VF vf)
{
    double localMax = 0.0;
    int x,y,z; FOREACH(x,y,z)
    {
        localMax = MAX(localMax,MathManager_LocalMagnitude(vf,x,y,z));
    }
    double max=0;
    FIND_GLOBAL_MAX(localMax,max);    
    return max;
}

//******************************************************************************
// MathManager_Average
// Compute the Average value across all processors
//******************************************************************************
double MathManager_AverageMagnitude(VF vf)
{
    double localSum = 0.0;
    int x,y,z; FOREACH(x,y,z) localSum+=MathManager_LocalMagnitude(vf,x,y,z);;
    
    double sum=0;
    FIND_GLOBAL_SUM(localSum,sum);    
    return sum*norm;
}

//******************************************************************************
// MathManager_RMSVectorField
// Compute the Average value across all processors
//******************************************************************************
double MathManager_RMSVectorField(VF vf)
{
    int x,y,z; 
    double localSum = 0.0;
    FOREACH(x,y,z) 
    {
        double magSquared = 0;
        int i; FOREACH_DIM(i) magSquared+=SQR(vf[i][x][y][z]);
        localSum+=magSquared;
    }
    
    double sum=0;
    FIND_GLOBAL_SUM(localSum,sum);    
    return sqrt(sum*norm);
}
//******************************************************************************
// MathManager_RMSScalarField
// Compute the Average value across all processors
//******************************************************************************
double MathManager_RMSScalarField(SF sf)
{
    int x,y,z; 
    double localSum = 0.0;
    FOREACH(x,y,z) localSum+=SQR(sf[x][y][z]);
    
    double sum=0;
    FIND_GLOBAL_SUM(localSum,sum);    
    return sqrt(sum*norm);
}

//******************************************************************************
// MathManager_RandomizeValues
//******************************************************************************
void MathManager_RandomizeValues(SF sf, double mean, double range)
{
    int x,y,z; FOREACH(x,y,z) sf[x][y][z] = range*(1.0*rand()/RAND_MAX - 0.5);
    double average = MathManager_Average(sf);
    FOREACH(x,y,z) sf[x][y][z]+= mean - average;
}

//******************************************************************************
// MathManager_ClearScalarField
//******************************************************************************
void MathManager_ClearScalarField(SF sf)
{
    memset(sf[0][0],0,totalSize*sizeof(double));
}

//******************************************************************************
// MathManager_Scale
//******************************************************************************
void MathManager_Scale(SF sf,double val)
{
    int x,y,z; FOREACH(x,y,z) sf[x][y][z]*=val;
}

//******************************************************************************
// MathManager_SetScalarField
//******************************************************************************
void MathManager_SetScalarField(SF sf,double val)
{
    int x,y,z; FOREACH(x,y,z) sf[x][y][z]=val;
}

//******************************************************************************
// MathManager_ClearVectorField
//******************************************************************************
void MathManager_ClearVectorField(VF vf, int d)
{
    int i; for (i=0; i<d; i++) MathManager_ClearScalarField(vf[i]);
}

//******************************************************************************
// MathManager_ClearTensorField
//******************************************************************************
void MathManager_ClearTensorField(TF tf,int td2,int td1)
{
    int i,j;
    for(i=0; i<td2; i++) for(j=0; j<td1; j++)
    {
        MathManager_ClearScalarField(tf[i][j]);
    }
}

//******************************************************************************
// MathManager_ScalarField
//******************************************************************************
void MathManager_CopyScalarField(SF source,SF dest)
{
    memcpy(dest[0][0],source[0][0],totalSize*sizeof(double));
}

//******************************************************************************
// MathManager_CopyVectorField
//******************************************************************************
void MathManager_CopyVectorField(VF source,VF destination, int d)
{
    int i; for (i=0; i<d; i++) MathManager_CopyScalarField(source[i],destination[i]);
}

//******************************************************************************
// MathManager_CopyTensorField
//******************************************************************************
void MathManager_CopyTensorField(TF source,TF destination,int td2, int td1)
{
    int i,j; FOREACH_DIM(i) FOREACH_DIM(j)
        MathManager_CopyScalarField(source[i][j],destination[i][j]);
}

//******************************************************************************
// MathManager_SetVectorField
//******************************************************************************
void MathManager_SetVectorField(VF vf,int d,double val)
{
    int i; for (i=0; i<d; i++)  MathManager_SetScalarField(vf[i],val);
}

//******************************************************************************
// MathManager_InnerProdVectorTensor
//******************************************************************************
void MathManager_InnerProductVT(VF vf, TF tf, VF v_dot_t)
{
    int c,r,x,y,z;
    
    FOREACH_DIM(c) FOREACH(x,y,z)
    {
        v_dot_t[c][x][y][x]=0;
        FOREACH_DIM(r) v_dot_t[c][x][y][z]+=vf[r][x][y][z]*tf[r][c][x][y][z];
    }
}

//******************************************************************************
// MathManager_ComputeBoundingBox
//******************************************************************************
IBOX MathManager_ComputeBoundingBox(SF field)
{
    GP min, max;
    int i; FOREACH_DIM(i) 
    {
        min[i] = gridSize[i];
        max[i] = 0;
    }
    
    GP r; FOREACH_GRIDPOINT(r) if (field[r[0]][r[1]][r[2]]!=0)
    {
        FOREACH_DIM(i) min[i]=MIN(min[i],r[i]);
        FOREACH_DIM(i) max[i]=MAX(max[i],r[i]);
    }
    
    IBOX bounds;
    bounds.ll[0]=min[0]; bounds.ll[1]=min[1]; 
    bounds.lr[0]=max[0]; bounds.lr[1]=min[1]; 
    bounds.ul[0]=min[0]; bounds.ul[1]=max[1]; 
    bounds.ur[0]=max[0]; bounds.ur[1]=max[1]; 
    printf("bounds: x=[%d %d] y=[%d %d] \n",min[0],max[0],min[1],max[1]);
    return bounds;
}

//******************************************************************************
// MathManager_CenterOfMass
//******************************************************************************
void MathManager_CenterOfMass(SF sf, V3 cm)
{
    int i; 
    FOREACH_DIM(i) cm[i]=0.0;
    double total = MathManager_IntegrateScalarField(sf);
    int x,y,z; FOREACH(x,y,z)
    {
        cm[0]+=x*sf[x][y][z];
        cm[1]+=y*sf[x][y][z];
        cm[2]+=z*sf[x][y][z];
    }
    FOREACH_DIM(i)cm[i]/=total;
}

//******************************************************************************
// MathManager_LinearInterpXYZ
//******************************************************************************
double MathManager_LinearInterpXYZ(SF sf,double x,double y,double z)
{
    int x0=floor(x); int y0=floor(y); int z0=floor(z);
    int x1=x0+1; int y1=y0+1; int z1=z0+1;
    double xf = x-x0;
    double yf = y-y0;
    double zf = z-z0;
    x0=IMOD(x0,nx); y0=IMOD(y0,ny); z0=IMOD(z0,nz);
    x1=IMOD(x1,nx); y1=IMOD(y1,ny); z1=IMOD(z1,nz);

    //interp in x dir
    double a1 = sf[x0][y0][z0]*(1.0-xf) + sf[x1][y0][z0]*xf;
    double a2 = sf[x0][y1][z0]*(1.0-xf) + sf[x1][y1][z0]*xf;
    double a3 = sf[x0][y0][z1]*(1.0-xf) + sf[x1][y0][z1]*xf;
    double a4 = sf[x0][y1][z1]*(1.0-xf) + sf[x1][y1][z1]*xf;

    // interp in y dir
    double b1 = a1*(1.0-yf) + a2*yf;
    double b2 = a3*(1.0-yf) + a4*yf;

    //interp in z dir
    double c1 = b1*(1.0-zf) + b2*zf;
   // if (c1>1.0) printf ("warning c1=%.3f\n>1. xf=%.2f, yf=%.2f zf=%.2f",c1,xf,yf,zf);
    return c1;
}


void MathManager_RotateGridpoint(GP gp,V3 origin,V3 result, double a)
{
    double dx= gp[0]-origin[0];
    double dy= gp[1]-origin[1];
    double dx1=cos(a)*dx - sin(a)*dy;
    double dy1=sin(a)*dx + cos(a)*dy;
    result[0] = dx1 + origin[0];
    result[1] = dy1 + origin[1];
}

//******************************************************************************
// MathManager_RotateField
//******************************************************************************
void MathManager_RotateField(SF sf1,SF sf2,V3 angle,V3 pivot,IBOX b,IBOX* rb)
{
    //rotate the particle bounding box
    BOX b1;
    MathManager_RotateGridpoint(b.ll,pivot,b1.ll,angle[2]);
    MathManager_RotateGridpoint(b.lr,pivot,b1.lr,angle[2]);
    MathManager_RotateGridpoint(b.ul,pivot,b1.ul,angle[2]);
    MathManager_RotateGridpoint(b.ur,pivot,b1.ur,angle[2]);
    
    rb->ll[0] = MIN4(b1.ll[0],b1.lr[0],b1.ul[0],b1.ur[0]);
    rb->ll[1] = MIN4(b1.ll[1],b1.lr[1],b1.ul[1],b1.ur[1]);
    
    rb->lr[0] = MAX4(b1.ll[0],b1.lr[0],b1.ul[0],b1.ur[0])+1;
    rb->lr[1] = MIN4(b1.ll[1],b1.lr[1],b1.ul[1],b1.ur[1]);
    
    rb->ul[0] = MIN4(b1.ll[0],b1.lr[0],b1.ul[0],b1.ur[0]);
    rb->ul[1] = MAX4(b1.ll[1],b1.lr[1],b1.ul[1],b1.ur[1])+1;

    rb->ur[0] = MAX4(b1.ll[0],b1.lr[0],b1.ul[0],b1.ur[0])+1;
    rb->ur[1] = MAX4(b1.ll[1],b1.lr[1],b1.ul[1],b1.ur[1])+1;
    
    double c= cos(angle[2]);
    double s= sin(angle[2]);
    MathManager_ClearScalarField(sf2);
    int x,y,z; for (x=rb->ll[0]; x<=rb->ur[0]; x++) for (y=rb->ll[1]; y<=rb->ur[1]; y++) FORZ(z)
    {
        double dx= x-pivot[0];
        double dy= y-pivot[1];
        double x1 = c*dx + s*dy + pivot[0];
        double y1 =-s*dx + c*dy + pivot[1];
        double z1 = z;
        sf2[IMOD(x,nx)][IMOD(y,ny)][IMOD(z,nz)]= MathManager_LinearInterpXYZ(sf1,x1,y1,z1);
    }
}

//******************************************************************************
// MathManager_TranslateField
//******************************************************************************
void MathManager_TranslateField(SF sf1,SF sf2, V3 r, V3 pivot,IBOX b)
{
    double dnx= r[0]*nx/L[0] - pivot[0];
    double dny= r[1]*ny/L[1] - pivot[1];
    double dnz= r[2]*nz/L[2] - pivot[2];
    
    //translate the extreme points of the particle bounding box
    IBOX tb;
    tb.ll[0]= b.ll[0] + dnx;
    tb.ll[1]= b.ll[1] + dny;
    tb.ur[0]= b.ur[0] + dnx + 1;
    tb.ur[1]= b.ur[1] + dny + 1;
    
   // printf("bounds: x=[%d %d] y=[%d %d] \n",b.ll[0],b.ur[0],b.ll[1],b.ur[1]);
    //printf("tbounds: x=[%d %d] y=[%d %d] \n",tb.ll[0],tb.ur[0],tb.ll[1],tb.ur[1]);

    MathManager_ClearScalarField(sf2);
    int x,y,z; for (x=tb.ll[0]; x<=tb.ur[0]; x++) for (y=tb.ll[1]; y<=tb.ur[1]; y++) FORZ(z)
    {
        sf2[IMOD(x,nx)][IMOD(y,ny)][IMOD(z,nz)]= MathManager_LinearInterpXYZ(sf1,x-dnx,y-dny,z-dnz);
    }
}

//******************************************************************************
// MathManager_VectorFieldMagnitude
//******************************************************************************
void MathManager_Convolve(SF sf1, SF sf2, SF result)
{
    int freqs[D],i;
    double k[D];
    
    MathManager_CopyScalarField(sf1,work2);
    MathManager_CopyScalarField(sf2,work3);
    
    fftw_complex* sf1Tilda = (fftw_complex*)sf1[0][0];
    fftw_complex* sf2Tilda = (fftw_complex*)sf2[0][0];
    fftw_complex* rTilda = (fftw_complex*)result[0][0];

    FOURIER_TRANSFORM(sf1,work1);
    FOURIER_TRANSFORM(sf2,work1);

    FOREACHFREQ(freqs)
    {
        i = FINDEX(freqs);
        GETWAVENUMS(freqs,k); 	 
        
        rTilda[i].re=norm*(sf1Tilda[i].re*sf2Tilda[i].re - sf1Tilda[i].im*sf2Tilda[i].im);
        rTilda[i].im=norm*(sf1Tilda[i].re*sf2Tilda[i].im + sf1Tilda[i].im*sf2Tilda[i].re);
    }
    INVERSE_TRANSFORM(result,work1);
}

//******************************************************************************
// MathManager_IntegralScalarField
// Compute the Average value across all processors
//******************************************************************************
double MathManager_IntegrateScalarField(SF sf)
{
    int x,y,z; 
    double localSum = 0.0;
    FOREACH(x,y,z) localSum+=sf[x][y][z];
    
    double sum=0;
    FIND_GLOBAL_SUM(localSum,sum);
    return sum;
}

//******************************************************************************
// MathManager_InnerProdTensorVector
//******************************************************************************
void MathManager_InnerProductTV(TF tf, VF vf, VF t_dot_v)
{
    int c,r,x,y,z;
    
    FOREACH_DIM(r) FOREACH(x,y,z)
    {
        t_dot_v[r][x][y][x]=0;
        FOREACH_DIM(c) t_dot_v[r][x][y][z]+=vf[c][x][y][z]*tf[r][c][x][y][z];
    }
}

//******************************************************************************
// MathManager_InnerProdTensorVector
//******************************************************************************
void MathManager_InnerProductTT(TF tf1, TF tf2, TF t1_dot_t2)
{
    int i,j,k;
    int x,y,z;
    FOREACH_TINDEX(i,k,x,y,z)
    {
        t1_dot_t2[i][k][x][y][x]=0;
        
        FOREACH_DIM(j) 
            t1_dot_t2[i][k][x][y][z]+=
            tf1[i][j][x][y][z]*tf2[j][k][x][y][z];
    }
}

//******************************************************************************
// MathManager_SumScalarFields
//******************************************************************************
void MathManager_SumScalarFields(SF sf1, double a, SF sf2, SF result)
{
    int x,y,z;
    FOREACH(x,y,z) result[x][y][z]=sf1[x][y][z] + a*sf2[x][y][z];
}

//******************************************************************************
// MathManager_SumVectorFields
//******************************************************************************
void MathManager_SumVectorFields(VF vf1, double a, VF vf2, VF result)
{
    int i,x,y,z; FOREACH_VINDEX(i,x,y,z) 
        result[i][x][y][z]=vf1[i][x][y][z] + a*vf2[i][x][y][z];
}

//******************************************************************************
// MathManager_SumTensorFields
//******************************************************************************
void MathManager_SumTensorFields(TF tf1, double a, TF tf2, TF result)
{
    int i,j,x,y,z; FOREACH_TINDEX(i,j,x,y,z) 
        result[i][j][x][y][z]= tf1[i][j][x][y][z] + a*tf2[i][j][x][y][z];
}

//******************************************************************************
// MathManager_Transpose
//******************************************************************************
void MathManager_Transpose(TF tf)
{
    int i,j;
    SF swap;
    
    FOREACH_LOWER(i,j)
    {
        swap=tf[i][j]; 
        tf[i][j]=tf[j][i];
        tf[j][i]=swap;
    }
}

//******************************************************************************
// MathManager_GradVector
//******************************************************************************
void   MathManager_GradVector(VF vf,TF tf)
{
    int i; FOREACH_DIM(i) MathManager_Gradient(vf[i],tf[i]);
    MathManager_Transpose(tf);
}

//******************************************************************************
// MathManager_ApplyExponentialFilter
//******************************************************************************
void MathManager_ApplyExponentialFilter(SF sf,double order)
{
    int freqs[D],index;
    double k[D],k2;

    fftw_complex* sfTilda = (fftw_complex*)sf[0][0];
    FOURIER_TRANSFORM(sf,work1);
    FOREACHFREQ(freqs)
    {
        index = FINDEX(freqs);
        GETWAVENUMS(freqs,k); 	 
        k2 = SQUARE_MAG(k);
        double power = DBL_MIN_EXP*pow(k2/k2Max,order);
        double temp=norm*exp(power);
        sfTilda[index].re*=temp;
        sfTilda[index].im*=temp;
    }
    INVERSE_TRANSFORM(sf,work1);
}

//******************************************************************************
// MathManager_VectorFieldMagnitude
//******************************************************************************
void   MathManager_VectorFieldMagnitude(VF vf,SF sf)
{
    int x,y,z;
    FOREACH(x,y,z) sf[x][y][z]=MathManager_LocalMagnitude(vf,x,y,z);
}

//******************************************************************************
// MathManager_GaussianRandomField
//******************************************************************************
void MathManager_GaussianRandomField(SF sf, double targetMean, double amplitude)
{
    int x,y,z; FOREACH(x,y,z)
    {
        sf[x][y][z] = amplitude*MathManager_GaussianRandomNum();
    }
    double mean = MathManager_Average(sf);
    FOREACH(x,y,z)
    {
        sf[x][y][z]+= targetMean - mean;
    }
}

//******************************************************************************
// MathManager_GaussianNoise
// Generate 2 gaussian random numbers using the polar Box-Muller algorithm
//******************************************************************************
double MathManager_GaussianRandomNum()
{
    double u1,u2,s;
    do
    {
        u1= 2.0*rand()/RAND_MAX - 1.0;
        u2 = 2.0*rand()/RAND_MAX - 1.0;
        s = u1*u1 + u2*u2;
    }
    while(s >= 1.0 );
    
    double r1 = sqrt( -2.0*log(s)/s)* u1;
    //double r2 = sqrt( -2.0*log(s)/s)* u2;
    
    return r1;
}

//******************************************************************************
// MathManager_VecDotGradVector
//******************************************************************************
void MathManager_VectorDotGradVector(VF v1,VF v2,VF result)
{
    int i;
    FOREACH_DIM(i)
    {
        MathManager_Gradient(v2[i],workVector);
        MathManager_VectorDotVector(v1,workVector,result[i]);
    }
}

//******************************************************************************
// MathManager_VecDotGradTensor
//******************************************************************************
void MathManager_VectorDotGradTensor(VF v1,TF t1,TF result)
{
    int i,j;
    FOREACH_DIM(i) FOREACH_DIM(j)
    {
        MathManager_Gradient(t1[i][j],workVector);
        MathManager_VectorDotVector(v1,workVector,result[i][j]);
    }
}

