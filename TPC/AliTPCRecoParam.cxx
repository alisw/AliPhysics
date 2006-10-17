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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with TPC reconstruction parameters                                  //
//                                                                           //  
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTPCRecoParam.h"

ClassImp(AliTPCRecoParam)




//_____________________________________________________________________________
AliTPCRecoParam::AliTPCRecoParam():
  fCtgRange(1.05),       
  fMaxSnpTracker(0.95),
  fMaxSnpTrack(0.999),
  fBYMirror(kTRUE),
  fFirstBin(0),
  fLastBin(-1),
  fBCalcPedestal(kFALSE),
  fBDoUnfold(kTRUE),
  fDumpAmplitudeMin(100),
  fMaxNoise(3.),
  fMaxC(0.3),
  fBSpecialSeeding(kFALSE),
  fBKinkFinder(kTRUE)
{
  //
  // constructor
  //
}

//_____________________________________________________________________________
AliTPCRecoParam::~AliTPCRecoParam() 
{
  //
  // destructor
  //  
}




AliTPCRecoParam *AliTPCRecoParam::GetLowFluxParam(){
  //
  // make default reconstruction  parameters for low  flux env.
  //
  AliTPCRecoParam *param = new AliTPCRecoParam;
  param->fCtgRange = 10;
  param->fFirstBin = 0;
  param->fLastBin  = 1000;
  return param;
}

AliTPCRecoParam *AliTPCRecoParam::GetHighFluxParam(){
  //
  // make reco parameters for high flux env.
  //
  AliTPCRecoParam *param = new AliTPCRecoParam;
  param->fCtgRange = 1.05;
  param->fFirstBin = 0;
  param->fLastBin  = 1000;
  return param;
}

AliTPCRecoParam *AliTPCRecoParam::GetLaserTestParam(Bool_t bPedestal){
  //
  // special setting for laser
  //
  AliTPCRecoParam *param = new AliTPCRecoParam;
  param->fCtgRange = 10.05;
  param->fFirstBin = 0;
  param->fLastBin  = 1000;
  param->fBCalcPedestal = bPedestal;
  param->fBDoUnfold     = kFALSE;
  param->fBKinkFinder   = kFALSE;
  param->fMaxSnpTracker = 0.98;
  param->fMaxC          = 0.02;
  param->fBSpecialSeeding = kTRUE;
  param->fBYMirror      = kFALSE;
  return param;
}

AliTPCRecoParam *AliTPCRecoParam::GetCosmicTestParam(Bool_t bPedestal){
  //
  // special setting for cosmic 
  // 
  AliTPCRecoParam *param = new AliTPCRecoParam;
  param->fCtgRange = 10.05;    // full TPC
  param->fFirstBin = 60;
  param->fLastBin  = 1000;
  param->fBCalcPedestal = bPedestal;
  param->fBDoUnfold     = kFALSE;
  param->fBSpecialSeeding = kTRUE;
  param->fMaxC          = 0.07;
  param->fBKinkFinder   = kFALSE;
  param->fBYMirror      = kFALSE;
  return param;
}



