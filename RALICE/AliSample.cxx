/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.2  1999/09/29 09:24:28  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////
// Class AliSample
// Perform statistics on various multi-dimensional data samples
// A data sample can be filled using the "Enter" and/or "Remove" functions,
// whereas the "Reset" function resets the complete sample to 'empty'.
// The info which can be extracted from a certain data sample are the
// sum, mean, variance, sigma, covariance and correlation.
// The "Info" function provides all statistics data for a certain sample.
// The variables for which these stat. parameters have to be calculated
// are indicated by the index of the variable which is passed as an
// argument to the various member functions.
// The index convention for a data point (x,y) is : x=1  y=2
//
// Example :
// ---------
// For an AliSample s a data point (x,y) can be entered as s.Enter(x,y) and
// the mean_x can be obtained as s.GetMean(1) whereas the mean_y is obtained
// via s.GetMean(2).
// The correlation between x and y is available via s.GetCor(1,2).
// The x-statistics are obtained via s.Info(1), y-statistics via s.Info(2),
// and the covariance and correlation between x and y via s.Info(1,2).
// All statistics of a sample are obtained via s.Info().
//
//--- Author: Nick van Eijndhoven 30-mar-1996 CERN Geneva
///////////////////////////////////////////////////////////////////////////

#include "AliSample.h"
 
AliSample::AliSample()
{
// Creation of an Aliample object and resetting the statistics values
// The dimension is initialised to maximum
 fDim=fMaxdim;
 fNames[0]='X';
 fNames[1]='Y';
 fNames[2]='Z';
 fN=0;
 Reset();
}
///////////////////////////////////////////////////////////////////////////
AliSample::~AliSample()
{
// Default destructor
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Reset()
{
// Resetting the statistics values for a certain Sample object
// Dimension is NOT changed
 fN=0;
 for (Int_t i=0; i<fDim; i++)
 {
  fSum[i]=0.;
  fMean[i]=0.;
  fVar[i]=0.;
  fSigma[i]=0.;
  for (Int_t j=0; j<fDim; j++)
  {
   fSum2[i][j]=0.;
   fCov[i][j]=0.;
   fCor[i][j]=0.;
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Enter(Float_t x)
{
// Entering a value into a 1-dim. sample
// In case of first entry the dimension is set to 1
 if (fN == 0)
 {
  fDim=1;
  fNames[0]='X';
  fNames[1]='-';
  fNames[2]='-';
 }
 if (fDim != 1)
 {
  cout << " *AliSample::enter* Error : Not a 1-dim sample." << endl;
 }
 else
 {
  fN+=1;
  fSum[0]+=x;
  fSum2[0][0]+=x*x;
  Compute();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Remove(Float_t x)
{
// Removing a value from a 1-dim. sample
 if (fDim != 1)
 {
  cout << " *AliSample::remove* Error : Not a 1-dim sample." << endl;
 }
 else
 {
  fN-=1;
  fSum[0]-=x;
  fSum2[0][0]-=x*x;
  Compute();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Enter(Float_t x,Float_t y)
{
// Entering a pair (x,y) into a 2-dim. sample
// In case of first entry the dimension is set to 2
 if (fN == 0)
 {
  fDim=2;
  fNames[0]='X';
  fNames[1]='Y';
  fNames[2]='-';
 }
 if (fDim != 2)
 {
  cout << " *AliSample::enter* Error : Not a 2-dim sample." << endl;
 }
 else
 {
  fN+=1;
  fSum[0]+=x;
  fSum[1]+=y;
  fSum2[0][0]+=x*x;
  fSum2[0][1]+=x*y;
  fSum2[1][0]+=y*x;
  fSum2[1][1]+=y*y;
  Compute();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Remove(Float_t x,Float_t y)
{
// Removing a pair (x,y) from a 2-dim. sample
 if (fDim != 2)
 {
  cout << " *AliSample::remove* Error : Not a 2-dim sample." << endl;
 }
 else
 {
  fN-=1;
  fSum[0]-=x;
  fSum[1]-=y;
  fSum2[0][0]-=x*x;
  fSum2[0][1]-=x*y;
  fSum2[1][0]-=y*x;
  fSum2[1][1]-=y*y;
  Compute();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Enter(Float_t x,Float_t y,Float_t z)
{
// Entering a set (x,y,z) into a 3-dim. sample
// In case of first entry the dimension is set to 3
 if (fN == 0)
 {
  fDim=3;
  fNames[0]='X';
  fNames[1]='Y';
  fNames[2]='Z';
 }
 if (fDim != 3)
 {
  cout << " *AliSample::enter* Error : Not a 3-dim sample." << endl;
 }
 else
 {
  fN+=1;
  fSum[0]+=x;
  fSum[1]+=y;
  fSum[2]+=z;
  fSum2[0][0]+=x*x;
  fSum2[0][1]+=x*y;
  fSum2[0][2]+=x*z;
  fSum2[1][0]+=y*x;
  fSum2[1][1]+=y*y;
  fSum2[1][2]+=y*z;
  fSum2[2][0]+=z*x;
  fSum2[2][1]+=z*y;
  fSum2[2][2]+=z*z;
  Compute();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Remove(Float_t x,Float_t y,Float_t z)
{
// Removing a set (x,y,z) from a 3-dim. sample
 if (fDim != 3)
 {
  cout << " *AliSample::remove* Error : Not a 3-dim sample." << endl;
 }
 else
 {
  fN-=1;
  fSum[0]-=x;
  fSum[1]-=y;
  fSum[2]-=z;
  fSum2[0][0]-=x*x;
  fSum2[0][1]-=x*y;
  fSum2[0][2]-=x*z;
  fSum2[1][0]-=y*x;
  fSum2[1][1]-=y*y;
  fSum2[1][2]-=y*z;
  fSum2[2][0]-=z*x;
  fSum2[2][1]-=z*y;
  fSum2[2][2]-=z*z;
  Compute();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Compute()
{
// Computation of the various statistical values
// after each entering or removing action on a certain sample
 Float_t rn=fN;
 for (Int_t k=0; k<fDim; k++)
 {
  fMean[k]=fSum[k]/rn;
  fVar[k]=(fSum2[k][k]/rn)-(fMean[k]*fMean[k]);
  if (fVar[k] < 0.) fVar[k]=0.;
  fSigma[k]=sqrt(fVar[k]);
 }
 for (Int_t i=0; i<fDim; i++)
 {
  for (Int_t j=0; j<fDim; j++)
  {
   fCov[i][j]=(fSum2[i][j]/rn)-(fMean[i]*fMean[j]);
   Float_t sigij=fSigma[i]*fSigma[j];
   if (sigij != 0.) fCor[i][j]=fCov[i][j]/sigij;
  }
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSample::GetDimension()
{
// Provide the dimension of a certain sample
 return fDim;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSample::GetN()
{
// Provide the number of entries of a certain sample
 return fN;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSample::GetSum(Int_t i)
{
// Provide the sum of a certain variable
 if (fDim < i)
 {
  cout << " *AliSample::sum* Error : Dimension less than " << i << endl;
  return 0.;
 }
 else
 {
 return fSum[i-1];
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSample::GetMean(Int_t i)
{
// Provide the mean of a certain variable
 if (fDim < i)
 {
  cout << " *AliSample::mean* Error : Dimension less than " << i << endl;
  return 0.;
 }
 else
 {
 return fMean[i-1];
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSample::GetVar(Int_t i)
{
// Provide the variance of a certain variable
 if (fDim < i)
 {
  cout << " *AliSample::var* Error : Dimension less than " << i << endl;
  return 0.;
 }
 else
 {
 return fVar[i-1];
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSample::GetSigma(Int_t i)
{
// Provide the standard deviation of a certain variable
 if (fDim < i)
 {
  cout << " *AliSample::sigma* Error : Dimension less than " << i << endl;
  return 0.;
 }
 else
 {
 return fSigma[i-1];
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSample::GetCov(Int_t i,Int_t j)
{
// Provide the covariance between variables i and j
 if ((fDim < i) || (fDim < j))
 {
  Int_t k=i;
  if (j > i) k=j;
  cout << " *AliSample::cov* Error : Dimension less than " << k << endl;
  return 0.;
 }
 else
 {
 return fCov[i-1][j-1];
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSample::GetCor(Int_t i,Int_t j)
{
// Provide the correlation between variables i and j
 if ((fDim < i) || (fDim < j))
 {
  Int_t k=i;
  if (j > i) k=j;
  cout << " *AliSample::cor* Error : Dimension less than " << k << endl;
  return 0.;
 }
 else
 {
 return fCor[i-1][j-1];
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Info()
{
// Printing of statistics of all variables
 for (Int_t i=0; i<fDim; i++)
 {
 cout << " " << fNames[i] << " : N = " << fN;
 cout << " Sum = " << fSum[i] << " Mean = " << fMean[i];
 cout << " Var = " << fVar[i] << " Sigma = " << fSigma[i] << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Info(Int_t i)
{
// Printing of statistics of ith variable
 if (fDim < i)
 {
  cout << " *AliSample::Info(i)* Error : Dimension less than " << i << endl;
 }
 else
 {
  cout << " " << fNames[i-1] << " : N = " << fN;
  cout << " Sum = " << fSum[i-1] << " Mean = " << fMean[i-1];
  cout << " Var = " << fVar[i-1] << " Sigma = " << fSigma[i-1] << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Info(Int_t i,Int_t j)
{
// Printing of covariance and correlation between variables i and j
 if ((fDim < i) || (fDim < j))
 {
  Int_t k=i;
  if (j > i) k=j;
  cout << " *AliSample::Info(i,j)* Error : Dimension less than " << k << endl;
 }
 else
 {
  cout << " " << fNames[i-1] << "-" << fNames[j-1] << " correlation :";
  cout << " Cov. = " << fCov[i-1][j-1] << " Cor. = " << fCor[i-1][j-1] << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
