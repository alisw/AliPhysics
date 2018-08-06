/**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Rudiger Haake <ruediger.haake@cern.ch>, Yale          *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************/


// AliRoot/EMCal system
#include "AliCaloRawAnalyzerKStandardFast.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliLog.h"

// Standard libraries
#include <stdexcept>
#include <Math/MinimizerOptions.h>


// Root system
#include "TF1.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TMath.h"

#include "AliEMCALRawResponse.h"

using namespace std;

/// \cond CLASSIMP
ClassImp( AliCaloRawAnalyzerKStandardFast ) ;
/// \endcond

///
/// Constructor
//_______________________________________________________________________
AliCaloRawAnalyzerKStandardFast::AliCaloRawAnalyzerKStandardFast() : AliCaloRawAnalyzerFitter("Chi Square ( kStandardFast )", "KStandardFast")
{
  fAlgo = Algo::kStandardFast;

  // Define fitter function
  fTf1 = new TF1("signalFit", RawResponseFunction, 0, TIMEBINS , 2);
  fTf1->SetParameters(10.,5.); //set all defaults once, just to be safe
  fTf1->SetParNames("amp","t0");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);

  fSignal =  new TGraph(TIMEBINS); 
}

///
/// Destructor
//_______________________________________________________________________
AliCaloRawAnalyzerKStandardFast::~AliCaloRawAnalyzerKStandardFast()
{
  //delete fTf1; // already done in base class
  if(fSignal)
    delete fSignal;
}

///
/// Evaluation Amplitude and TOF
//_______________________________________________________________________
AliCaloFitResults
AliCaloRawAnalyzerKStandardFast::Evaluate( const vector<AliCaloBunchInfo>  &bunchlist,
                                       UInt_t altrocfg1, UInt_t altrocfg2 )
{
  Float_t pedEstimate  = 0;
  short maxADC = 0;
  Int_t first = 0;
  Int_t last = 0;
  Int_t bunchIndex = 0;
  Float_t ampEstimate = 0;
  short timeEstimate  = 0;
  Float_t time = 0;
  Float_t amp=0;
  Float_t chi2 = 0;
  Int_t ndf = 0;
  Bool_t fitDone = kFALSE;

  int nsamples = PreFitEvaluateSamples( bunchlist, altrocfg1, altrocfg2, bunchIndex, ampEstimate, 
                                       maxADC, timeEstimate, pedEstimate, first, last,   (int)fAmpCut ); 
  
  
  if (ampEstimate >= fAmpCut  ) 
  { 
    time = timeEstimate; 
    Int_t timebinOffset = bunchlist.at(bunchIndex).GetStartBin() - (bunchlist.at(bunchIndex).GetLength()-1); 
    amp = ampEstimate; 
    
    if ( nsamples > 1 && maxADC< OVERFLOWCUT ) 
    { 
      FitRaw(first, last, amp, time, chi2, fitDone);
      time += timebinOffset;
      timeEstimate += timebinOffset;
      ndf = nsamples - 2;
    }
  } 
  if ( fitDone ) 
  { 
    Float_t ampAsymm = (amp - ampEstimate)/(amp + ampEstimate);
    Float_t timeDiff = time - timeEstimate;
    
    if ( (TMath::Abs(ampAsymm) > 0.1) || (TMath::Abs(timeDiff) > 2) ) 
    {
      amp     = ampEstimate;
      time    = timeEstimate; 
      fitDone = kFALSE;
    } 
  }  

  if (amp >= fAmpCut ) 
  { 
    if ( ! fitDone) 
    { 
      amp += (0.5 - gRandom->Rndm()); 
    }
    time = time * TIMEBINWITH; 
    time -= fL1Phase;
    
    return AliCaloFitResults( -99, -99, fAlgo , amp, time,
                             (int)time, chi2, ndf, Ret::kDummy );
  }
  return AliCaloFitResults( Ret::kInvalid, Ret::kInvalid );
}


///
/// Fits the raw signal time distribution
//_______________________________________________________________________
void
 AliCaloRawAnalyzerKStandardFast::FitRaw( Int_t firstTimeBin, Int_t lastTimeBin,
                                      Float_t & amp, Float_t & time, Float_t & chi2, Bool_t & fitDone) const
{ 
  // Changes wrt to kStandard fitter:
  // * Fit function not created in each FitRaw() call
  // * Signal histogram fSignal not created in each FitRaw() call
  // * Adjusted fit range to datapoints
  // * Fit function simplification
  // * Use fit option N0
  // * Minimizer strategy adjusted

  Int_t nsamples = lastTimeBin - firstTimeBin + 1;
  fitDone = kFALSE;
  if ( nsamples < 3 )  return;  

  fSignal->Set(nsamples);
  for (int i=0; i<nsamples; i++) 
  {
    Int_t timebin = firstTimeBin + i;    
    fSignal->SetPoint(i, timebin, GetReversed(timebin)); 
  }

  // Initial fit function
  fTf1->SetRange(lastTimeBin, firstTimeBin);
  fTf1->SetParameter(1, time);
  fTf1->SetParameter(0, amp);
  fTf1->SetParLimits(0, 0.5*amp, 2*amp );
  fTf1->SetParLimits(1, time - 4, time + 4); 
  
  try 
  {
    fSignal->Fit(fTf1, "QROW N0"); // Note option 'W': equal errors on all points, 'R': Use range we set earlier, N0: Do not plot result, Q: quite mode
    amp  = fTf1->GetParameter(0);
    time = fTf1->GetParameter(1);
    chi2 = fTf1->GetChisquare();

    fitDone = kTRUE;
  }
  catch (const std::exception & e) 
  {
    AliError( Form("TH1 Fit exception %s", e.what()) ); 
    // stay with default amp and time in case of exception, i.e. no special action required
    fitDone = kFALSE;
  }

  return;
}

///
/// Approximate response function of the EMCal electronics.
/// Simplified version, adapted from AliEMCALRawResponse
//_______________________________________________________________________
Double_t AliCaloRawAnalyzerKStandardFast::RawResponseFunction(Double_t *x, Double_t *par)
{
  Double_t signal = 0.;
  Double_t xx     = ( x[0] - par[1] + TAU ) / TAU;
  
  if(xx > 0)
    signal = par[0] * TMath::Power(xx , ORDER) * TMath::Exp(ORDER * (1 - xx )) ;

  return signal;
}

