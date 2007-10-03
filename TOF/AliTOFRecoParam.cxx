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
// Class with TOF reconstruction parameters                                  //
//                                                                           //  
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTOFRecoParam.h"

ClassImp(AliTOFRecoParam)

//_____________________________________________________________________________
AliTOFRecoParam::AliTOFRecoParam():
  TObject(),       
  fTimeZero(kFALSE),       
  fTimeZerofromT0(kFALSE),       
  fTimeZerofromTOF(kFALSE),       
  fTimeWalkCorr(kFALSE),       
  fApplyPbPbCuts(kFALSE),       
  fWindowSizeMaxY(50.),
  fWindowSizeMaxZ(35.),
  fWindowScaleFact(5.),
  fDistanceCut(3.),
  fSensRadius(378.5),
  fStepSize(0.1),
  fMaxChi2(10.),
  fMaxChi2TRD(150.),
  fTimeResolution(80.),
  fTimeNSigma(5.)
{
  //
  // constructor
  //
}
//_____________________________________________________________________________
AliTOFRecoParam::AliTOFRecoParam(const AliTOFRecoParam &p):
  TObject(),       
  fTimeZero(kFALSE),       
  fTimeZerofromT0(kFALSE),       
  fTimeZerofromTOF(kFALSE),       
  fTimeWalkCorr(kFALSE),       
  fApplyPbPbCuts(kFALSE),       
  fWindowSizeMaxY(50.),
  fWindowSizeMaxZ(35.),
  fWindowScaleFact(5.),
  fDistanceCut(3.),
  fSensRadius(378.5),
  fStepSize(0.1),
  fMaxChi2(10.),
  fMaxChi2TRD(150.),
  fTimeResolution(80.),
  fTimeNSigma(5.)
 { 
  //copy Ctor

   fTimeZero=p.fTimeZero;       
   fTimeZerofromT0=p.fTimeZerofromT0;
   fTimeZerofromTOF=p.fTimeZerofromTOF;       
   fTimeWalkCorr=p.fTimeWalkCorr;       
   fApplyPbPbCuts=p.fApplyPbPbCuts;       
   fWindowSizeMaxY=p.fWindowSizeMaxY;
   fWindowSizeMaxZ=p.fWindowSizeMaxZ;
   fWindowScaleFact=p.fWindowScaleFact;
   fDistanceCut=p.fDistanceCut;
   fSensRadius=p.fSensRadius;
   fStepSize=p.fStepSize;
   fMaxChi2=p.fMaxChi2;
   fMaxChi2TRD=p.fMaxChi2TRD;
   fTimeResolution=p.fTimeResolution;
   fTimeNSigma=p.fTimeNSigma;   
}
//_____________________________________________________________________________
AliTOFRecoParam& AliTOFRecoParam::operator=(const AliTOFRecoParam &p)
{
  //
  // assign. operator
  //
   this->fTimeZero=p.fTimeZero;       
   this->fTimeZerofromT0=p.fTimeZerofromT0;
   this->fTimeZerofromTOF=p.fTimeZerofromTOF;       
   this->fTimeWalkCorr=p.fTimeWalkCorr;       
   this->fApplyPbPbCuts=p.fApplyPbPbCuts;       
   this->fWindowSizeMaxY=p.fWindowSizeMaxY;
   this->fWindowSizeMaxZ=p.fWindowSizeMaxZ;
   this->fDistanceCut=p.fDistanceCut;
   this->fWindowScaleFact=p.fWindowScaleFact;
   this->fStepSize=p.fStepSize;
   this->fSensRadius=p.fSensRadius;
   this->fMaxChi2=p.fMaxChi2;
   this->fMaxChi2TRD=p.fMaxChi2TRD;
   this->fTimeResolution=p.fTimeResolution;
   this->fTimeNSigma=p.fTimeNSigma;   
   return *this;
}
//_____________________________________________________________________________
AliTOFRecoParam::~AliTOFRecoParam() 
{
  //
  // destructor
  //  
}

//_____________________________________________________________________________
AliTOFRecoParam *AliTOFRecoParam::GetPbPbparam(){
  //
  // set default reconstruction parameters for PbPb.
  //
  AliTOFRecoParam *param = new AliTOFRecoParam();
  param->fApplyPbPbCuts = kTRUE;
  param->fWindowScaleFact = 5.;
  param->fDistanceCut = 3.;
  return param;
}

//_____________________________________________________________________________
AliTOFRecoParam *AliTOFRecoParam::GetPPparam(){
  //
  // set default reconstruction parameters for PP.
  //
  AliTOFRecoParam *param = new AliTOFRecoParam();
  param->fApplyPbPbCuts = kFALSE;
  param->fWindowScaleFact = 5.;
  param->fDistanceCut = 10.;
  return param;
}
