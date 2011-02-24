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

/* $Id$ */

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

