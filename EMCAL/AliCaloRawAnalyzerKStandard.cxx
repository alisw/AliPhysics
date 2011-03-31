/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille <p.t.hille@fys.uio.no>                *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to p.t.hille@fys.uio.no                             *
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
// Extraction of amplitude and peak position
// FRom CALO raw data using
// least square fit for the
// Moment assuming identical and 
// independent errors (equivalent with chi square)
// 

#include "AliCaloRawAnalyzerKStandard.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliLog.h"
#include "TMath.h"
#include <stdexcept>
#include <iostream>
#include "TF1.h"
#include "TGraph.h"
#include "TRandom.h"

#include "AliEMCALRawResponse.h"


using namespace std;

ClassImp( AliCaloRawAnalyzerKStandard )


AliCaloRawAnalyzerKStandard::AliCaloRawAnalyzerKStandard() : AliCaloRawAnalyzerFitter("Chi Square ( kStandard )", "KStandard")
{
  
  fAlgo = Algo::kStandard;
}


AliCaloRawAnalyzerKStandard::~AliCaloRawAnalyzerKStandard()
{
  //  delete fTf1;
}


AliCaloFitResults
AliCaloRawAnalyzerKStandard::Evaluate( const vector<AliCaloBunchInfo>  &bunchlist, const UInt_t altrocfg1,  const UInt_t altrocfg2 )
{
  //Evaluation Amplitude and TOF
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

	
//____________________________________________________________________________ 
void
 AliCaloRawAnalyzerKStandard::FitRaw(const Int_t firstTimeBin, const Int_t lastTimeBin, Float_t & amp, Float_t & time, Float_t & chi2, Bool_t & fitDone) const 
{ 
  // Fits the raw signal time distribution
  int nsamples = lastTimeBin - firstTimeBin + 1;
  fitDone = kFALSE;
  if (nsamples < 3) { return; } 
  
  TGraph *gSig =  new TGraph( nsamples); 
 
  for (int i=0; i<nsamples; i++) 
    {
      Int_t timebin = firstTimeBin + i;    
      gSig->SetPoint(i, timebin, GetReversed(timebin)); 
    }
  
  TF1 * signalF = new TF1("signal", AliEMCALRawResponse::RawResponseFunction, 0, TIMEBINS , 5);
  
  signalF->SetParameters(10.,5., TAU  ,ORDER,0.); //set all defaults once, just to be safe
  signalF->SetParNames("amp","t0","tau","N","ped");
  signalF->FixParameter(2,TAU); 
  signalF->FixParameter(3,ORDER); 
  signalF->FixParameter(4, 0); 
  signalF->SetParameter(1, time);
  signalF->SetParameter(0, amp);
  signalF->SetParLimits(0, 0.5*amp, 2*amp );
  signalF->SetParLimits(1, time - 4, time + 4); 
      
  try {			
    gSig->Fit(signalF, "QROW"); // Note option 'W': equal errors on all points
    amp  = signalF->GetParameter(0); 
    time = signalF->GetParameter(1);
    chi2 = signalF->GetChisquare();
    fitDone = kTRUE;
  }
  catch (const std::exception & e) {
    AliError( Form("TGraph Fit exception %s", e.what()) ); 
    // stay with default amp and time in case of exception, i.e. no special action required
    fitDone = kFALSE;
  }

  delete signalF;
  delete gSig; // delete TGraph
  return;
}


