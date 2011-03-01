/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//
// Dielectron helper functions wrapped in a namespace
// 
//
// Authors: 
//   Jens Wiechula <Jens.Wiechula@cern.ch> 
// 




#include <TError.h>
#include <TMath.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TVectorD.h>

#include <AliVEvent.h>
#include <AliKFParticle.h>

#include "AliDielectronHelper.h"

//_____________________________________________________________________________
TVectorD* AliDielectronHelper::MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make logarithmic binning
  // the user has to delete the array afterwards!!!
  //
  
  //check limits
  if (xmin<1e-20 || xmax<1e-20){
    Error("AliDielectronHelper::MakeLogBinning","For Log binning xmin and xmax must be > 1e-20. Using linear binning instead!");
    return AliDielectronHelper::MakeLinBinning(nbinsX, xmin, xmax);
  }
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t expMax=TMath::Log(last/first);
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first*TMath::Exp(expMax/nbinsX*(Double_t)i);
  }
  return binLim;
}

//_____________________________________________________________________________
TVectorD* AliDielectronHelper::MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make linear binning
  // the user has to delete the array afterwards!!!
  //
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t binWidth=(last-first)/nbinsX;
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first+binWidth*(Double_t)i;
  }
  return binLim;
}

//_____________________________________________________________________________
TVectorD* AliDielectronHelper::MakeArbitraryBinning(const char* bins)
{
  //
  // Make arbitrary binning, bins separated by a ','
  //
  TString limits(bins);
  if (limits.IsNull()){
    Error("AliDielectronHelper::MakeArbitraryBinning","Bin Limit string is empty, cannot add the variable");
    return 0x0;
  }
  
  TObjArray *arr=limits.Tokenize(",");
  Int_t nLimits=arr->GetEntries();
  if (nLimits<2){
    Error("AliDielectronHelper::MakeArbitraryBinning","Need at leas 2 bin limits, cannot add the variable");
    delete arr;
    return 0x0;
  }
  
  TVectorD *binLimits=new TVectorD(nLimits);
  for (Int_t iLim=0; iLim<nLimits; ++iLim){
    (*binLimits)[iLim]=(static_cast<TObjString*>(arr->At(iLim)))->GetString().Atof();
  }
  
  delete arr;
  return binLimits;
}

//_____________________________________________________________________________
void AliDielectronHelper::RotateKFParticle(AliKFParticle * kfParticle,Double_t angle, const AliVEvent * const ev){
  // Before rotate needs to be moved to position 0,0,0, ; move back after rotation
  if (!kfParticle) return;
  Double_t dx = 0.;
  Double_t dy = 0.;
  Double_t dz = 0.;

  if (ev){
    dx = ev->GetPrimaryVertex()->GetX()-0.;
    dy = ev->GetPrimaryVertex()->GetY()-0.;
    dz = ev->GetPrimaryVertex()->GetZ()-0.;
  }
  
  kfParticle->X() = kfParticle->GetX() - dx;
  kfParticle->Y() = kfParticle->GetY() - dy;
  kfParticle->Z() = kfParticle->GetZ() - dz;
  
  
  // Rotate the kf particle
  Double_t c = cos(angle);
  Double_t s = sin(angle);
  
  Double_t mA[8][ 8];
  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<8; j++){
      mA[i][j] = 0;
    }
  }
  for( int i=0; i<8; i++ ){
    mA[i][i] = 1;
  }
  mA[0][0] =  c;  mA[0][1] = s;
  mA[1][0] = -s;  mA[1][1] = c;
  mA[3][3] =  c;  mA[3][4] = s;
  mA[4][3] = -s;  mA[4][4] = c;
  
  Double_t mAC[8][8];
  Double_t mAp[8];
  
  for( Int_t i=0; i<8; i++ ){
    mAp[i] = 0;
    for( Int_t k=0; k<8; k++){
      mAp[i]+=mA[i][k] * kfParticle->GetParameter(k);
    }
  }
  
  for( Int_t i=0; i<8; i++){
    kfParticle->Parameter(i) = mAp[i];
  }
  
  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<8; j++ ){
      mAC[i][j] = 0;
      for( Int_t k=0; k<8; k++ ){
        mAC[i][j]+= mA[i][k] * kfParticle->GetCovariance(k,j);
      }
    }
  }
  
  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<=i; j++ ){
      Double_t xx = 0;
      for( Int_t k=0; k<8; k++){
        xx+= mAC[i][k]*mA[j][k];
      }
      kfParticle->Covariance(i,j) = xx;
    }
  }
  
  kfParticle->X() = kfParticle->GetX() + dx;
  kfParticle->Y() = kfParticle->GetY() + dy;
  kfParticle->Z() = kfParticle->GetZ() + dz;
  
}

