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
  fTimeZero(kFALSE),       
  fTimeZerofromT0(kFALSE),       
  fTimeZerofromTOF(kFALSE),       
  fTimeWalkCorr(kFALSE),       
  fApplyPbPbCuts(kFALSE),       
  fWindowSizeMaxY(50.),
  fWindowSizeMaxZ(35.),
  fWindowScaleFact(3.),
  fDistanceCut(3.),
  fSensRadius(378.),
  fStepSize(0.1),
  fMaxChi2(150.),
  fTimeResolution(80.),
  fTimeNSigma(150.)
{
  //
  // constructor
  //
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
  param->fWindowScaleFact = 3.;
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
