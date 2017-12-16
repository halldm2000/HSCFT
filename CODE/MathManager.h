//******************************************************************************
//  MathManager.h
//******************************************************************************
#ifndef MathManager_H
#define MathManager_H
#include "Constants.h"

void    MathManager_ApplyExponentialFilter(SF sf,double order);
double  MathManager_Average               (SF sf);
double  MathManager_AverageMagnitude      (VF vf);
void    MathManager_CenterOfMass          (SF sf, V3 cm);
void    MathManager_ClearScalarField      (SF sf);
void    MathManager_ClearVectorField      (VF vf,int d);
void    MathManager_ClearTensorField      (TF tf,int td2,int td1);
IBOX    MathManager_ComputeBoundingBox    (SF field);
void    MathManager_Convolve              (SF field, SF kernel, SF result);
void    MathManager_CopyScalarField       (SF source,SF dest);
void    MathManager_CopyVectorField       (VF source,VF dest,int d);
void    MathManager_CopyTensorField       (TF source,TF dest,int td2, int td1);
void    MathManager_CrossProduct          (V3 a,V3 b,V3 result);
void    MathManager_Divergence            (VF vf,SF sf);
void    MathManager_DivTensorField        (TF tf,VF vf);
void    MathManager_DownsampleData        (SF oldSF,SF newSF);
void    MathManager_VectorDotVector       (VF vf1, VF vf2, SF result);
void    MathManager_VectorDotTensor       (VF vf, TF tf, VF result);
void    MathManager_Gradient              (SF sf,VF gradSF);
void    MathManager_GradVector            (VF vf,TF tf);
double  MathManager_GaussianRandomNum();
void    MathManager_GaussianRandomField   (SF sf, double targetMean, double amplitude);
void    MathManager_InnerProductVT        (VF vf, TF tf, VF v_dot_t);
void    MathManager_InnerProductTV        (TF tf, VF vf, VF t_dot_v);
void    MathManager_InnerProductTT        (TF tf1, TF tf2, TF t1_dot_t2);
double  MathManager_IntegrateScalarField  (SF sf);
void    MathManager_Laplacian             (SF sf,SF result);
double  MathManager_LinearInterpXYZ       (SF sf,double x,double y,double z);
double  MathManager_LocalMagnitude        (VF vf,int x,int y,int z);
double  MathManager_Maximum               (SF sf);
double  MathManager_MaximumMagnitude      (VF vf);
double  MathManager_MeasureMeanLengthScale(SF sf);
double  MathManager_Minimum               (SF sf);
int     MathManager_MinIntegerInVector    (int* vector,int dim);
double  MathManager_MinDoubleInVector     (VEC vector,int dim);
double  MathManager_RMSScalarField        (SF sf);
double  MathManager_RMSVectorField        (VF vf);
void    MathManager_RandomizeValues       (SF sf, double mean, double range);
void    MathManager_ResampleData          (SF oldField,SF newField);
void    MathManager_RotateGridpoint       (GP gp,V3 origin,V3 result, double a);
void    MathManager_RotateField           (SF sf1,SF sf2,V3 angle,V3 pivot,IBOX b,IBOX* rb);
void    MathManager_Scale                 (SF sf, double val);
void    MathManager_SetScalarField        (SF sf, double val);
void    MathManager_SetVectorField        (VF vf, int d, double val);
void    MathManager_SolvePoissonEquation  (SF sf,SF result);
double  MathManager_SumOfSquares          (SF sf);
double  MathManager_Sum                   (SF sf);
void    MathManager_SumScalarFields       (SF sf1, double a, SF sf2, SF result);
void    MathManager_SumVectorFields       (VF vf1, double a, VF vf2, VF result);
void    MathManager_SumTensorFields       (TF tf1, double a, TF tf2, TF result);
void    MathManager_TranslateField        (SF sf1,SF sf2, V3 r, V3 pivot,IBOX b);
void    MathManager_UpsampleData          (SF oldSF,SF newSF);
void    MathManager_VectorFieldMagnitude  (VF vf, SF sf);
void    MathManager_VectorDotGradVector   (VF v1,VF v2,VF result);
void    MathManager_VectorDotGradTensor   (VF v1,TF t1,TF result);

#endif
