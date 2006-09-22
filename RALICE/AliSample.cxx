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

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class AliSample
// Perform statistics on various multi-dimensional data samples.
// A data sample can be filled using the "Enter" and/or "Remove" functions,
// whereas the "Reset" function resets the complete sample to 'empty'.
// The info which can be extracted from a certain data sample are the
// sum, mean, variance, sigma, covariance and correlation.
// The "Data" function provides all statistics data for a certain sample.
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
// The x-statistics are obtained via s.Data(1), y-statistics via s.Data(2),
// and the covariance and correlation between x and y via s.Data(1,2).
// All statistics of a sample are obtained via s.Data().
//
//--- Author: Nick van Eijndhoven 30-mar-1996 CERN Geneva
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliSample.h"
#include "Riostream.h"
 
ClassImp(AliSample) // Class implementation to enable ROOT I/O
 
AliSample::AliSample()
{
// Creation of an Aliample object and resetting the statistics values
// The dimension is initialised to maximum
 fDim=fMaxdim;
 fNames[0]='X';
 fNames[1]='Y';
 fNames[2]='Z';
 fN=0;
 fStore=0;
 fX=0;
 fY=0;
 fZ=0;
 fArr=0;
 Reset();
}
///////////////////////////////////////////////////////////////////////////
AliSample::~AliSample()
{
// Default destructor
 if (fX)
 {
  delete fX;
  fX=0;
 }
 if (fY)
 {
  delete fY;
  fY=0;
 }
 if (fZ)
 {
  delete fZ;
  fZ=0;
 }
 if (fArr)
 {
  delete fArr;
  fArr=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Reset()
{
// Resetting the statistics values for a certain Sample object
// Dimension and storage flag are NOT changed
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

 // Set storage arrays to initial size
 if (fX) fX->Set(10);
 if (fY) fY->Set(10);
 if (fZ) fZ->Set(10);
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

  // Store all entered data when storage mode has been selected 
  if (fStore)
  {
   if (!fX) fX=new TArrayF(10);
   if (fX->GetSize() < fN) fX->Set(fN+10);
   fX->AddAt(x,fN-1);
  }

  // Compute the various statistics
  Compute();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Remove(Float_t x)
{
// Removing a value from a 1-dim. sample

 if (!fN) return;

 if (fDim != 1)
 {
  cout << " *AliSample::remove* Error : Not a 1-dim sample." << endl;
 }
 else
 {
  fN-=1;
  fSum[0]-=x;
  fSum2[0][0]-=x*x;

  // Remove data entry from the storage
  if (fStore && fX)
  {
   Float_t val=0;
   for (Int_t i=0; i<=fN; i++)
   {
    val=fX->At(i);
    if (fabs(x-val)>1.e-10) continue;

    for (Int_t j=i+1; j<=fN; j++)
    {
     val=fX->At(j);
     fX->AddAt(val,j-1);
    }
   }
  }

  // Compute the various statistics
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

  // Store all entered data when storage mode has been selected 
  if (fStore)
  {
   if (!fX) fX=new TArrayF(10);
   if (fX->GetSize() < fN) fX->Set(fN+10);
   fX->AddAt(x,fN-1);
   if (!fY) fY=new TArrayF(10);
   if (fY->GetSize() < fN) fY->Set(fN+10);
   fY->AddAt(y,fN-1);
  }

  // Compute the various statistics
  Compute();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Remove(Float_t x,Float_t y)
{
// Removing a pair (x,y) from a 2-dim. sample

 if (!fN) return;

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

  // Remove data entry from the storage
  if (fStore && fX && fY)
  {
   Float_t val=0;
   for (Int_t i=0; i<=fN; i++)
   {
    val=fX->At(i);
    if (fabs(x-val)>1.e-10) continue;
    val=fY->At(i);
    if (fabs(y-val)>1.e-10) continue;

    for (Int_t j=i+1; j<=fN; j++)
    {
     val=fX->At(j);
     fX->AddAt(val,j-1);
     val=fY->At(j);
     fY->AddAt(val,j-1);
    }
   }
  }

  // Compute the various statistics
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

  // Store all entered data when storage mode has been selected 
  if (fStore)
  {
   if (!fX) fX=new TArrayF(10);
   if (fX->GetSize() < fN) fX->Set(fN+10);
   fX->AddAt(x,fN-1);
   if (!fY) fY=new TArrayF(10);
   if (fY->GetSize() < fN) fY->Set(fN+10);
   fY->AddAt(y,fN-1);
   if (!fZ) fZ=new TArrayF(10);
   if (fZ->GetSize() < fN) fZ->Set(fN+10);
   fZ->AddAt(z,fN-1);
  }

  // Compute the various statistics
  Compute();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Remove(Float_t x,Float_t y,Float_t z)
{
// Removing a set (x,y,z) from a 3-dim. sample

 if (!fN) return;

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

  // Remove data entry from the storage
  if (fStore && fX && fY && fZ)
  {
   Float_t val=0;
   for (Int_t i=0; i<=fN; i++)
   {
    val=fX->At(i);
    if (fabs(x-val)>1.e-10) continue;
    val=fY->At(i);
    if (fabs(y-val)>1.e-10) continue;
    val=fZ->At(i);
    if (fabs(z-val)>1.e-10) continue;

    for (Int_t j=i+1; j<=fN; j++)
    {
     val=fX->At(j);
     fX->AddAt(val,j-1);
     val=fY->At(j);
     fY->AddAt(val,j-1);
     val=fZ->At(j);
     fZ->AddAt(val,j-1);
    }
   }
  }

  // Compute the various statistics
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
Int_t AliSample::GetDimension() const
{
// Provide the dimension of a certain sample
 return fDim;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSample::GetN() const
{
// Provide the number of entries of a certain sample
 return fN;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSample::GetSum(Int_t i) const
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
Float_t AliSample::GetMean(Int_t i) const
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
Float_t AliSample::GetVar(Int_t i) const
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
Float_t AliSample::GetSigma(Int_t i) const
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
Float_t AliSample::GetCov(Int_t i,Int_t j) const
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
Float_t AliSample::GetCor(Int_t i,Int_t j) const
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
void AliSample::Data()
{
// Printing of statistics of all variables
 for (Int_t i=0; i<fDim; i++)
 {
 cout << " " << fNames[i] << " : N = " << fN;
 cout << " Sum = " << fSum[i] << " Mean = " << fMean[i];
 cout << " Var = " << fVar[i] << " Sigma = " << fSigma[i];
 if (fStore)
 {
  cout << endl;
  cout << "     Minimum = " << GetMinimum(i+1);
  cout << " Maximum = " << GetMaximum(i+1);
  cout << " Median = " << GetMedian(i+1);
 }
 cout << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Data(Int_t i)
{
// Printing of statistics of ith variable
 if (fDim < i)
 {
  cout << " *AliSample::Data(i)* Error : Dimension less than " << i << endl;
 }
 else
 {
  cout << " " << fNames[i-1] << " : N = " << fN;
  cout << " Sum = " << fSum[i-1] << " Mean = " << fMean[i-1];
  cout << " Var = " << fVar[i-1] << " Sigma = " << fSigma[i-1];
  if (fStore)
  {
   cout << endl;
   cout << "     Minimum = " << GetMinimum(i);
   cout << " Maximum = " << GetMaximum(i);
   cout << " Median = " << GetMedian(i);
  }
  cout << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::Data(Int_t i,Int_t j) const
{
// Printing of covariance and correlation between variables i and j
 if ((fDim < i) || (fDim < j))
 {
  Int_t k=i;
  if (j > i) k=j;
  cout << " *AliSample::Data(i,j)* Error : Dimension less than " << k << endl;
 }
 else
 {
  cout << " " << fNames[i-1] << "-" << fNames[j-1] << " correlation :";
  cout << " Cov. = " << fCov[i-1][j-1] << " Cor. = " << fCor[i-1][j-1] << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSample::SetStoreMode(Int_t mode)
{
// Set storage mode for all entered data.
//
// mode = 0 : Entered data will not be stored
//        1 : All data will be stored as entered
//
// By default the storage mode is set to 0 in the constructor of this class.
// The default at invokation of this memberfunction is mode=1.
//
// For normal statistics evaluation (e.g. mean, sigma, covariance etc...)
// storage of entered data is not needed. This is the default mode
// of operation and is the most efficient w.r.t. cpu time and memory.
// However, when calculation of a median, minimum or maximum is required,
// then the data storage mode has be activated.
//
// Note : Activation of storage mode can only be performed before the
//        first data item is entered. 

 if (fN)
 {
  cout << " *AliSample::SetStore* Storage mode can only be set before first data." << endl;
 }
 else
 {
  if (mode==0 || mode==1) fStore=mode;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSample::GetStoreMode() const
{
// Provide the storage mode
 return fStore;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSample::GetMedian(Int_t i)
{
// Provide the median of a certain variable.
// For this functionality the storage mode has to be activated.

 if (fDim < i)
 {
  cout << " *AliSample::GetMedian* Error : Dimension less than " << i << endl;
  return 0;
 }

 if (!fStore)
 {
  cout << " *AliSample::GetMedian* Error : Storage of data entries was not activated." << endl;
  return 0;
 }

 if (fN<=0) return 0;

 Float_t median=0;

 if (fN==1)
 {
  if (i==1) median=fX->At(0);
  if (i==2) median=fY->At(0);
  if (i==3) median=fZ->At(0);
  return median;
 }

 // Prepare temp. array to hold the ordered values
 if (!fArr)
 {
  fArr=new TArrayF(fN);
 }
 else
 {
  if (fArr->GetSize() < fN) fArr->Set(fN);
 }

 // Order the values of the specified variable
 Float_t val=0;
 Int_t iadd=0;
 for (Int_t j=0; j<fN; j++)
 {
  if (i==1) val=fX->At(j);
  if (i==2) val=fY->At(j);
  if (i==3) val=fZ->At(j);

  iadd=0;
  if (j==0)
  {
   fArr->AddAt(val,j);
   iadd=1;
  }
  else
  {
   for (Int_t k=0; k<j; k++)
   {
    if (val>=fArr->At(k)) continue;
    // Put value in between the existing ones
    for (Int_t m=j-1; m>=k; m--)
    {
     fArr->AddAt(fArr->At(m),m+1);
    }
    fArr->AddAt(val,k);
    iadd=1;
    break;
   }

   if (!iadd)
   {
    fArr->AddAt(val,j);
   }
  }
 }

 median=0;
 Int_t index=fN/2;
 if (fN%2) // Odd number of entries
 {
  median=fArr->At(index);
 }
 else // Even number of entries
 {
  median=(fArr->At(index-1)+fArr->At(index))/2.;
 }
 return median;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSample::GetMinimum(Int_t i) const
{
// Provide the minimum value of a certain variable.
// For this functionality the storage mode has to be activated.

 if (fDim < i)
 {
  cout << " *AliSample::GetMinimum* Error : Dimension less than " << i << endl;
  return 0;
 }

 if (!fStore)
 {
  cout << " *AliSample::GetMinimum* Error : Storage of data entries was not activated." << endl;
  return 0;
 }

 if (fN<=0) return 0;

 Float_t min=0;

 if (i==1) min=fX->At(0);
 if (i==2) min=fY->At(0);
 if (i==3) min=fZ->At(0);

 for (Int_t k=1; k<fN; k++)
 {
  if (i==1 && fX->At(k)<min) min=fX->At(k);
  if (i==2 && fY->At(k)<min) min=fY->At(k);
  if (i==3 && fZ->At(k)<min) min=fZ->At(k);
 }

 return min;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSample::GetMaximum(Int_t i) const
{
// Provide the maxmum value of a certain variable.
// For this functionality the storage mode has to be activated.

 if (fDim < i)
 {
  cout << " *AliSample::GetMaximum* Error : Dimension less than " << i << endl;
  return 0;
 }

 if (!fStore)
 {
  cout << " *AliSample::GetMaximum* Error : Storage of data entries was not activated." << endl;
  return 0;
 }

 if (fN<=0) return 0;

 Float_t max=0;

 if (i==1) max=fX->At(0);
 if (i==2) max=fY->At(0);
 if (i==3) max=fZ->At(0);

 for (Int_t k=1; k<fN; k++)
 {
  if (i==1 && fX->At(k)>max) max=fX->At(k);
  if (i==2 && fY->At(k)>max) max=fY->At(k);
  if (i==3 && fZ->At(k)>max) max=fZ->At(k);
 }

 return max;
}
///////////////////////////////////////////////////////////////////////////
