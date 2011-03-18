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


#include "AliLog.h"

#include "AliTOFRecoParam.h"

ClassImp(AliTOFRecoParam)

//_____________________________________________________________________________
AliTOFRecoParam::AliTOFRecoParam():
  AliDetectorRecoParam(),       
  fTimeZero(kFALSE),       
  fTimeZerofromT0(kFALSE),       
  fTimeZerofromTOF(kFALSE),       
  fTimeWalkCorr(kFALSE),       
  fApplyPbPbCuts(kFALSE),       
  fWindowSizeMaxY(50.),
  fWindowSizeMaxZ(35.),
  fWindowScaleFact(5.),
  fDistanceCut(3.),
  fSensRadius(379.5),
  fStepSize(0.1),
  fMaxChi2(10.),
  fMaxChi2TRD(150.),
  fTimeResolution(80.),
  fTimeNSigma(5.),
  fMaxDeltaTime(2),
  fStartTimeType(1)// 1 -> AliESDpid::kTOF_T0
{
  //
  // constructor
  //
  SetNameTitle("TOF","TOF");
}
//_____________________________________________________________________________
AliTOFRecoParam::AliTOFRecoParam(const AliTOFRecoParam &p):
  AliDetectorRecoParam(p),       
  fTimeZero(p.fTimeZero),       
  fTimeZerofromT0(p.fTimeZerofromT0),       
  fTimeZerofromTOF(p.fTimeZerofromTOF),       
  fTimeWalkCorr(p.fTimeWalkCorr),       
  fApplyPbPbCuts(p.fApplyPbPbCuts),       
  fWindowSizeMaxY(p.fWindowSizeMaxY),
  fWindowSizeMaxZ(p.fWindowSizeMaxZ),
  fWindowScaleFact(p.fWindowScaleFact),
  fDistanceCut(p.fDistanceCut),
  fSensRadius(p.fSensRadius),
  fStepSize(p.fStepSize),
  fMaxChi2(p.fMaxChi2),
  fMaxChi2TRD(p.fMaxChi2TRD),
  fTimeResolution(p.fTimeResolution),
  fTimeNSigma(p.fTimeNSigma),
  fMaxDeltaTime(p.fMaxDeltaTime),
  fStartTimeType(p.fStartTimeType)
 { 
  //copy Ctor


}
//_____________________________________________________________________________
AliTOFRecoParam& AliTOFRecoParam::operator=(const AliTOFRecoParam &p)
{
  //
  // assign. operator
  //

  if (this == &p)
    return *this;
  
  AliDetectorRecoParam::operator=(p);
  fTimeZero=p.fTimeZero;       
  fTimeZerofromT0=p.fTimeZerofromT0;
  fTimeZerofromTOF=p.fTimeZerofromTOF;       
  fTimeWalkCorr=p.fTimeWalkCorr;       
  fApplyPbPbCuts=p.fApplyPbPbCuts;       
  fWindowSizeMaxY=p.fWindowSizeMaxY;
  fWindowSizeMaxZ=p.fWindowSizeMaxZ;
  fDistanceCut=p.fDistanceCut;
  fWindowScaleFact=p.fWindowScaleFact;
  fStepSize=p.fStepSize;
  fSensRadius=p.fSensRadius;
  fMaxChi2=p.fMaxChi2;
  fMaxChi2TRD=p.fMaxChi2TRD;
  fTimeResolution=p.fTimeResolution;
  fTimeNSigma=p.fTimeNSigma;   
  fMaxDeltaTime=p.fMaxDeltaTime;
  fStartTimeType=p.fStartTimeType;
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
  param->fMaxDeltaTime = 2;
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
  param->fMaxDeltaTime = 2;
  return param;
}

//_____________________________________________________________________________
void AliTOFRecoParam::PrintParameters() const
{
  //
  // Printing of the used TOF reconstruction parameters
  //

  AliInfo(Form(" Use Time Zero info in Pid: %i", fTimeZero));
  AliInfo(Form(" Use Time Zero as determined by T0: %i", fTimeZerofromT0));
  AliInfo(Form(" Use Time Zero as determined from TOF: %i",
	       fTimeZerofromTOF));
  AliInfo(Form(" Correct for signal time walk in z: %i", fTimeWalkCorr));
  AliInfo(Form(" Apply high flux cuts: %i", fApplyPbPbCuts));

  AliInfo(Form(" Cluster search window - size, Y: %f cm", fWindowSizeMaxY));
  AliInfo(Form(" Cluster search window - size, Z: %f cm", fWindowSizeMaxZ));
  AliInfo(Form(" Cluster search window - scale factor: %f",
	       fWindowScaleFact));

  AliInfo(Form(" Cut on the closest approach distance: %f", fDistanceCut));
  AliInfo(Form(" Average radius of sensitive volumes: %f cm", fSensRadius));
  AliInfo(Form(" Propagation step size: %f cm", fStepSize));
  AliInfo(Form(" Maximum X2 track-tof clusters: %f", fMaxChi2));

  AliInfo(Form(" Maximum X2 track-tof clusters (TRD): %f", fMaxChi2TRD));
  AliInfo(Form(" Time resolution for responce function in PID: %f ps",
	       fTimeResolution));
  AliInfo(Form(" N-Sigma Range used for responce function in PID: %f",
	       fTimeNSigma));
  AliInfo(Form(" Max time difference between two (or more) time measurements"
	       " coming from two (or more) neighbouring pads"
	       " to define if they are coming from the same ionization event"
	       " or not: %2d bins -> %f ps",
	       fMaxDeltaTime, fMaxDeltaTime*AliTOFGeometry::TdcBinWidth()))

}
