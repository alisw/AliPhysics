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

// Extraction of amplitude and peak position
// FRom CALO raw data using
// least square fit for the
// Moment assuming identical and 
// independent errors (equivalent with chi square)
// 

#include "AliCaloRawAnalyzerLMSOffline.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliLog.h"
#include "TMath.h"
#include <stdexcept>
#include <iostream>
#include "TF1.h"
#include "TGraph.h"
#include "AliEMCALRawUtils.h"
#include <TRandom.h>

//#include "AliCaloRawAnalyzerLMS.h"

#include "AliCaloRawAnalyzerFactory.h"

using namespace std;


//#define BAD 4  //CRAP PTH
 
ClassImp( AliCaloRawAnalyzerLMSOffline )


AliCaloRawAnalyzerLMSOffline::AliCaloRawAnalyzerLMSOffline() : AliCaloRawAnalyzer("Chi Square Fit OFFLINE", "LMSOffline"),
  fNoiseThreshold(0),
							       fRawAnalyzer(0), fSmearFactor()
{
  fRawAnalyzer = AliCaloRawAnalyzerFactory::CreateAnalyzer(kLMS);
  fAlgo = Algo::kLMSOffline;
}


AliCaloRawAnalyzerLMSOffline::~AliCaloRawAnalyzerLMSOffline()
{
  //  delete fTf1;
}


AliCaloFitResults
AliCaloRawAnalyzerLMSOffline::Evaluate( const vector<AliCaloBunchInfo>  &bunchlist, const UInt_t altrocfg1,  const UInt_t altrocfg2 )
{
  // fRawAnalyzer setup
  fRawAnalyzer->SetNsampleCut(5); // requirement for fits to be done, for the new methods
  fRawAnalyzer->SetOverflowCut(OVERFLOWCUT);
  fRawAnalyzer->SetAmpCut(fNoiseThreshold);
  fRawAnalyzer->SetFitArrayCut(fNoiseThreshold);
  
  //  fRawAnalyzer->SetIsZeroSuppressed(true); // TMP - should use stream->IsZeroSuppressed(), or altro cfg registers later

  fRawAnalyzer->SetIsZeroSuppressed(fIsZerosupressed); // TMP - should use stream->IsZeroSuppressed(), or altro cfg registers later


  AliCaloFitResults fitResults;
  Float_t time = 0; 
  Float_t amp  = 0; 
  short timeEstimate  = 0;
  Float_t ampEstimate = 0;
  Bool_t fitDone = kFALSE;
  fitResults = fRawAnalyzer->Evaluate( bunchlist, altrocfg1, altrocfg2 ); 
  amp          = fitResults.GetAmp();
  time         = fitResults.GetTime();
  timeEstimate = fitResults.GetMaxTimebin();
  ampEstimate  = fitResults.GetMaxSig();
  
  if (fitResults.GetStatus() == Ret::kFitPar) 
    {
      fitDone = kTRUE;
    }
  
  if ( fitDone ) 
    { // brief sanity check of fit results	    
      Float_t ampAsymm = (amp - ampEstimate)/(amp + ampEstimate);
      Float_t timeDiff = time - timeEstimate;
      
      if ( (TMath::Abs(ampAsymm) > 0.1) || (TMath::Abs(timeDiff) > 2) ) 
	{
	  amp     = ampEstimate;
	  time    = timeEstimate; 
	  fitDone = kFALSE;
	} 
    } // fitDone
  
  if (amp >= fNoiseThreshold) 
    { // something to be stored
      if ( ! fitDone) 
	{ // smear ADC with +- 0.5 uniform (avoid discrete effects)
	  fSmearFactor =  (0.5 - gRandom->Rndm() );
	  amp += fSmearFactor; // Rndm generates a number in ]0,1]
	}
      // go from time-bin units to physical time fgtimetrigger
      time = time *TIMEBINWITH; // skip subtraction of fgTimeTrigger?
      // subtract RCU L1 phase (L1Phase is in seconds) w.r.t. L0:
      //   time -= in.GetL1Phase();
      time -= fL1Phase;
      time = time/  TIMEBINWITH; 
    }
  fitResults.SetTime(time);
  fitResults.SetAmp(amp);
  return fitResults ; 
}

