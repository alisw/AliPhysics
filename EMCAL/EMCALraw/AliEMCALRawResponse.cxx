// -*- mode: c++ -*-

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "TF1.h"
#include "TMath.h"
#include <TRandom.h>

#include "AliEMCALRawResponse.h"
#include "AliCaloConstants.h"

using namespace CALO;
using namespace EMCAL;
using namespace ALTRO;

Double_t AliEMCALRawResponse::fgTimeTrigger  = 600E-9 ;   // the time of the trigger as approximately seen in the data
Int_t    AliEMCALRawResponse::fgThreshold         = 1;
Int_t    AliEMCALRawResponse::fgPedestalValue     = 0;  // pedestal value for digits2raw, default generate ZS data
Double_t AliEMCALRawResponse::fgFEENoise          = 3.; // 3 ADC channels of noise (sampled)



ClassImp(AliEMCALRawResponse)


AliEMCALRawResponse::AliEMCALRawResponse()
{
  //comment
}

AliEMCALRawResponse::~AliEMCALRawResponse()
{

}


Double_t 
AliEMCALRawResponse::RawResponseFunction(Double_t *x, Double_t *par)
{
  // step response of the emcal electronincs
  Double_t signal = 0.;
  Double_t tau    = par[2];
  Double_t n      = par[3];
  Double_t ped    = par[4];
  Double_t xx     = ( x[0] - par[1] + tau ) / tau ;
  if (xx <= 0) 
    signal = ped ;  
  else 
    {  
      signal = ped + par[0] * TMath::Power(xx , n) * TMath::Exp(n * (1 - xx )) ; 
    }
  return signal ;  
}


Bool_t  
AliEMCALRawResponse::RawSampledResponse(const Double_t dtime, const Double_t damp, Int_t * adcH, 
					    Int_t * adcL, const Int_t keyErr)
{
  // step response of the emcal electronincs
  Bool_t lowGain = kFALSE ; 
  TF1 signalF("signal", RawResponseFunction, 0, TIMEBINS, 5);
  signalF.SetParameter(0, damp) ; 
  signalF.SetParameter(1, (dtime + fgTimeTrigger)/ TIMEBINWITH) ; 
  signalF.SetParameter(2, TAU) ; 
  signalF.SetParameter(3, ORDER);
  signalF.SetParameter(4, fgPedestalValue);
	
  Double_t signal=0.0, noise=0.0;
  for (Int_t iTime = 0; iTime <  TIMEBINS; iTime++) {
    signal = signalF.Eval(iTime) ;  
    
    if(keyErr>0) {
      noise = gRandom->Gaus(0.,fgFEENoise);
      signal += noise; 
    }
	  
    adcH[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    if ( adcH[iTime] > MAXBINVALUE ){  // larger than 10 bits 
      adcH[iTime] = MAXBINVALUE ;
      lowGain = kTRUE ; 
    }
    signal /= HGLGFACTOR;
    adcL[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    if ( adcL[iTime] > MAXBINVALUE )  // larger than 10 bits 
      adcL[iTime] = MAXBINVALUE ;
  }
  return lowGain ; 
}

