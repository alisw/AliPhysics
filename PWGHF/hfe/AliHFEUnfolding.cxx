#include "AliHFEUnfolding.h"
#include <iostream>
// Response matrix is R(true value, measured value)

ClassImp(AliHFEUnfolding);

AliHFEUnfolding::AliHFEUnfolding()
{
  // Constructor
  fDataIsSet = kFALSE;
  fResponseIsSet = kFALSE;
  fLowestBin = 1;
  fHighestBin = 0;
  fDataVector = new TVectorD(1);
  fResponseMatrix = new TMatrixD(1,1);
  fCovarianceMatrix = new TMatrixD(1,1);
}

AliHFEUnfolding::~AliHFEUnfolding(void)
{
  // Desctructor
}

void AliHFEUnfolding::SetData(TH1D * DataHistogram)
{
  // Read in data histogram.
  fDataHistogram=DataHistogram;
  fDataIsSet=kTRUE;  
  if(fHighestBin==0) fHighestBin = fDataHistogram->GetNbinsX();
  fDataVector->ResizeTo(fHighestBin-fLowestBin+1);
  for(int i=fLowestBin;i<=fHighestBin;i++)
  { // Fill data Vector
    (*fDataVector)(i-fLowestBin) = fDataHistogram->GetBinContent(i);
  }
}

void AliHFEUnfolding::SetData(TH1D * DataHistogram, Int_t LowestBin, Int_t HighestBin) // Set Data Historgam and the range of bins considered in unfolding.
{
  // Read in data histogram. Unfolding is done from first bin to last bin
  fLowestBin=LowestBin;
  fHighestBin=HighestBin;
  SetData(DataHistogram);
}

void AliHFEUnfolding::SetResponseMatrix(TH2D * ResponseMatrix)
{
  // Response matrix must be normalized correctly beforehand.
  // It must also be invertible
  if(!fDataIsSet){
    printf("Set Data first, then response matrix!\n");
    return;
  }
  fResponseMatrixHistogram=ResponseMatrix;
  fResponseIsSet=kTRUE;
  fResponseMatrix->ResizeTo(fHighestBin-fLowestBin+1, fHighestBin-fLowestBin+1);
  fCovarianceMatrix->ResizeTo(fHighestBin-fLowestBin+1, fHighestBin-fLowestBin+1);
  for(int i=fLowestBin;i<=fHighestBin;i++)
    for(int j=fLowestBin;j<=fHighestBin;j++)
      { // Fill Response Matrix
	(*fResponseMatrix)(j-fLowestBin,i-fLowestBin) = fResponseMatrixHistogram->GetBinContent(i,j);
      }
}


TH1D * AliHFEUnfolding::Unfold(void)
{
  // Starts unfolding, returns unfolded histogram with errors from diagonal elements of covaiance matrix
  // Creates Covariance matrix which has to be read out using the GetCovarianceMatrix-function
  if(!fDataIsSet)
    {printf("No Data Histogram is set. Unfolding cancelled.\n"); return 0;}
  if(!fResponseIsSet)
    {printf("No Response Matrix Histogram is set. Unfolding cancelled.\n"); return 0;}
   
  // Invert ResponseMatrix
  TMatrixD InverseResponseMatrix(fHighestBin-fLowestBin+1, fHighestBin-fLowestBin+1);
  TVectorD UnfoldedSpectrum(fHighestBin-fLowestBin+1);
  InverseResponseMatrix=*fResponseMatrix;
  InverseResponseMatrix.Invert();
  
  UnfoldedSpectrum = InverseResponseMatrix * (*fDataVector);
  
  TH1D * UnfoldedSpectrumHistogram = (TH1D*) fDataHistogram->Clone(); UnfoldedSpectrumHistogram->SetName("UnfoldedSpectrumHistogram");
  for(int i=fLowestBin;i<=fHighestBin;i++)
  { // Save unfolded Spectrum to output Historam
    UnfoldedSpectrumHistogram->SetBinContent(i, UnfoldedSpectrum(i-fLowestBin));
  }
  
  double tempCov = 0.;
  for(int i=0;i<=fHighestBin-fLowestBin;i++)
    for(int j=0;j<=fHighestBin-fLowestBin;j++)
    {
      tempCov = 0;
      for(int k=0;k<=fHighestBin-fLowestBin;k++)
      {
	tempCov += InverseResponseMatrix(i,k) * InverseResponseMatrix(j,k) * fDataHistogram->GetBinError(k+fLowestBin) * fDataHistogram->GetBinError(k+fLowestBin);
      }
      (*fCovarianceMatrix)(i,j) = tempCov;
    }
  for(int i=fLowestBin;i<=fHighestBin;i++)
  {
    UnfoldedSpectrumHistogram->SetBinError(i, TMath::Sqrt((*fCovarianceMatrix)(i-fLowestBin, i-fLowestBin)));
  }
  
  return UnfoldedSpectrumHistogram;
}
