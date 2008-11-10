/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
 * Class does fitting on an histogram
 *
 * @file   AliHLTPHOSPhysicsAnalyzerPeakFitter.cxx
 * @author Oystein Djuvsland
 * @date
 * @brief  Fitter for PHOS HLT
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <Riostream.h>
#include "AliHLTPHOSPhysicsAnalyzerPeakFitter.h"
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"


ClassImp(AliHLTPHOSPhysicsAnalyzerPeakFitter);

AliHLTPHOSPhysicsAnalyzerPeakFitter::AliHLTPHOSPhysicsAnalyzerPeakFitter() : fGainLow(80), fGainHigh(5),
									     fRootHistPtr(0)
{
  //Constructor
}

AliHLTPHOSPhysicsAnalyzerPeakFitter::~AliHLTPHOSPhysicsAnalyzerPeakFitter()
{
  //Destructor
}

AliHLTPHOSPhysicsAnalyzerPeakFitter::AliHLTPHOSPhysicsAnalyzerPeakFitter(const AliHLTPHOSPhysicsAnalyzerPeakFitter&): fGainLow(80), fGainHigh(5),
														      fRootHistPtr(0)
{
  //Copy constructor
}


Int_t
AliHLTPHOSPhysicsAnalyzerPeakFitter::FitGaussian()
{
  //FitGaussian

  Int_t maxBin = fRootHistPtr->GetMaximumBin();
  Float_t binWidth = fRootHistPtr->GetBinWidth(maxBin);
  Float_t maxBinValue = (Float_t)(maxBin * binWidth);
  Float_t lowRange = maxBinValue - 0.05;
  Float_t highRange = maxBinValue + 0.05;

  TF1* gaussian = new TF1("gaussian", "gaus", 0.1, 0.2);
    
  fRootHistPtr->Fit(gaussian->GetName(), "", "",lowRange, highRange);
  
  return 0;

}

Int_t
AliHLTPHOSPhysicsAnalyzerPeakFitter::FitLorentzian()
{
  //FitLorentzia
  Int_t maxBin = fRootHistPtr->GetMaximumBin();
  Float_t binWidth = fRootHistPtr->GetBinWidth(maxBin);
  Float_t maxBinValue = (Float_t)(maxBin * binWidth);
  Double_t lowRange = maxBinValue - 0.03;
  Double_t highRange = maxBinValue + 0.03;

  //char* name = "lorentzian";
  
  TF1* lorentzian = new TF1("lorentzian", "([0]*1/TMath::Pi())*[1]/((x[0]-[2])*(x[0]-[2])+[1]*[1])", lowRange, highRange);

  Double_t params[3] = {fRootHistPtr->GetBinContent(maxBin)/20, 0.01, 0.135};
  lorentzian->SetParameters(params);

  fRootHistPtr->Fit(lorentzian->GetName(), "", "", lowRange, highRange);

  lorentzian->GetParameters(params);

//   TFile *outfile = new TFile(,"recreate");  
//   fRootHistPtr->Write();
//   outfile->Close();

  return 0;
}

