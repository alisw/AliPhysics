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
//
/*
  The reconstruction parameters are used in the AliTPCclustererMI and AliTPCtrackerMI
  
  They are retrieved:
  0. User speciefied it in reconstruction macro
  1. if (not 0) from OCDB  - AliTPCcalibDB::GetRecoParam(eventtype)
  2. if (not 0 or 1) default parameter - High flux enevironment used  

  FIXME:
  In the future  reconstruction parameters should be changed on event basis
  But for the moment, event types are still not defined 


  // Setting for systematic errors addition
  [0] - systematic RMSY
  [1] - systematic RMSZ
  [2] - systematic RMSSNP
  [3] - systematic RMSTheta
  [4] - systematic RMSCuravture -  systematic error in 1/cm not in 1/pt
  //
  //  How to add it example - 3 mm systematic error y, 3 cm systematic error z (drift)
  Double_t sysError[5]={0.3,3, 0.3/150., 3./150.,0.3/(150*150.)}
  param->SetSystematicError(sysError);

*/
                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTPCRecoParam.h"

ClassImp(AliTPCRecoParam)


Bool_t AliTPCRecoParam::fgUseTimeCalibration=kTRUE; // flag usage the time dependent calibration
                                      // to be switched off for pass 0 reconstruction
                                      // Use static function, other option will be to use 
                                      // additional specific storage ?

//_____________________________________________________________________________
AliTPCRecoParam::AliTPCRecoParam():
  AliDetectorRecoParam(),
  fUseHLTClusters(4),  // use HLTorRAW data
  fBClusterSharing(kTRUE),
  fCtgRange(1.05),       
  fMaxSnpTracker(0.95),
  fMaxSnpTrack(0.999),
  fUseOuterDetectors(kFALSE),
  fDumpSignal(kFALSE),
  fFirstBin(0),
  fLastBin(-1),
  fBCalcPedestal(kFALSE),
  fBDoUnfold(kTRUE),
  fDumpAmplitudeMin(100),
  fMaxNoise(2.),
  //
  fUseOnePadCluster(kTRUE),
  fUseHLTOnePadCluster(kFALSE),
  fMinMaxCutAbs(4.),
  fMinLeftRightCutAbs(6.),
  fMinUpDownCutAbs(6.),
  //
  fMinMaxCutSigma(4.),
  fMinLeftRightCutSigma(7.),
  fMinUpDownCutSigma(8.),
  fMaxC(0.3),
  fBSpecialSeeding(kFALSE),
  fBKinkFinder(kTRUE),
  fLastSeedRowSec(120),
  fSeedGapPrim(6),
  fSeedGapSec(6),
  fUseFieldCorrection(2),      // use field correction
  fUseComposedCorrection(kFALSE),      // use field correction
  fUseRPHICorrection(0),      // use rphi correction
  fUseRadialCorrection(0),    // use radial correction
  fUseQuadrantAlignment(0),   // use quadrant alignment
  fUseSectorAlignment(0),     // use sector alignment
  fUseDriftCorrectionTime(1), // use drift correction time
  fUseDriftCorrectionGY(1),   // use drif correction global y
  fUseGainCorrectionTime(0),  // use gain correction time
  fUseExBCorrection(1),  // use ExB correction
  fUseMultiplicityCorrectionDedx(kTRUE), // use Dedx multiplicity correction
  fUseAlignmentTime(kTRUE),              // use time dependent alignment correction
  //
  fUseTotCharge(kTRUE),          // switch use total or max charge
  fMinFraction(0.01),           // truncated mean - lower threshold
  fMaxFaction(0.7),            // truncated mean - upper threshold
  fNeighborRowsDedx(2),           // neighbour rows for below threshold dEdx calculation
  fUseTOFCorrection(kTRUE),
  fUseSystematicCorrelation(kTRUE)
{
  //
  // constructor
  //
  SetName("TPC");
  SetTitle("TPC");
  for (Int_t i=0;i<5;i++) fSystematicErrors[i]=0;
  fCutSharedClusters[0]=0.5; // maximal allowed fraction of shared clusters - shorter track
  fCutSharedClusters[1]=0.25; // maximal allowed fraction of shared clusters - longer  track
  fClusterMaxRange[0]=1;     // y - pad      range
  fClusterMaxRange[1]=1;     // z - time bin range
  fKinkAngleCutChi2[0]=9;    // angular cut for kink finder - to create a kink
                             // ~ about 5 % rate  for high pt kink finder
  fKinkAngleCutChi2[1]=12;    // angular cut for kink finder - to use the partial track                             // form kink 
                             // ~ about 2 % rate  for high pt kink finder
}

//_____________________________________________________________________________
AliTPCRecoParam::~AliTPCRecoParam() 
{
  //
  // destructor
  //  
}

void AliTPCRecoParam::Print(const Option_t* /*option*/) const{
  //
  //
  //
  AliTPCRecoParam::Dump();
  printf("Systematic errors:\n");
  const char * cherrs[5]={"sy=","sz=","ssnp=","stheta=","s1pt="};
  for (Int_t i=0; i<5; i++){
    printf("%s%f\n",cherrs[i],fSystematicErrors[i]);
  }
}


AliTPCRecoParam *AliTPCRecoParam::GetLowFluxParam(){
  //
  // make default reconstruction  parameters for low  flux env.
  //
  AliTPCRecoParam *param = new AliTPCRecoParam;
  param->fCtgRange = 10;
  param->fFirstBin = 0;
  param->fLastBin  = 1000;
  param->SetName("Low Flux");
  param->SetTitle("Low Flux");
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
  param->fUseTotCharge=kFALSE;
  param->SetName("High Flux");
  param->SetTitle("High Flux");
  return param;
}

AliTPCRecoParam *AliTPCRecoParam::GetHLTParam(){
  //
  // make reco parameters for high flux env.
  //
  AliTPCRecoParam *param = new AliTPCRecoParam;
  param->fCtgRange = 1.05;
  param->fFirstBin = 80;
  param->fLastBin  = 1000;  
  param->fMaxSnpTracker = 0.9; 
  param->fMaxC          = 0.06; 
  //
  param->SetName("Hlt Param");
  param->SetTitle("Hlt Param"); 
  param->fBKinkFinder   = kFALSE;
  return param;
}

AliTPCRecoParam *AliTPCRecoParam::GetLaserTestParam(Bool_t bPedestal){
  //
  // special setting for laser
  //
  AliTPCRecoParam *param = new AliTPCRecoParam;
  param->fDumpSignal=kTRUE;
  param->fCtgRange = 10.05;
  param->fFirstBin = 0;
  param->fLastBin  = 1000;
  param->fBCalcPedestal = bPedestal;
  param->fBDoUnfold     = kFALSE;
  param->fDumpAmplitudeMin = 150;
  param->fBKinkFinder   = kFALSE;
  param->fMaxSnpTracker = 0.98;
  param->fMaxC          = 0.02;
  param->fBSpecialSeeding = kTRUE;
  param->fUseTOFCorrection=kFALSE;
  param->fUseHLTClusters=1; // always RAW data
  //
  //
  param->SetName("Laser Flux");
  param->SetTitle("Laser Flux");
  return param;
}

AliTPCRecoParam *AliTPCRecoParam::GetCosmicTestParam(Bool_t bPedestal){
  //
  // special setting for cosmic 
  // 
  AliTPCRecoParam *param = new AliTPCRecoParam;
  param->fDumpSignal=kTRUE;
  param->fCtgRange = 10.05;    // full TPC
  param->fFirstBin = 60;
  param->fLastBin  = 1000;
  param->fBCalcPedestal = bPedestal;
  param->fBDoUnfold     = kFALSE;
  param->fBSpecialSeeding = kTRUE;
  param->fMaxC          = 0.07;
  param->fBKinkFinder   = kFALSE;
  param->fUseTOFCorrection =kFALSE;
  param->SetName("Cosmic Flux");
  param->SetTitle("Cosmic Flux");

  return param;
}


Bool_t  AliTPCRecoParam::GetUseTimeCalibration(){ 
  //
  // get
  //
  return fgUseTimeCalibration;
}
void    AliTPCRecoParam::SetUseTimeCalibration(Bool_t useTimeCalibration) {
  //
  // set 
  //
  fgUseTimeCalibration = useTimeCalibration;
}

