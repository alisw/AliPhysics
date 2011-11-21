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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Parameter class for the TRD reconstruction                               //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"

#include "AliTRDrecoParam.h"

ClassImp(AliTRDrecoParam)


//______________________________________________________________
AliTRDrecoParam::AliTRDrecoParam()
  :AliDetectorRecoParam()
  ,fkdNchdy(12.)
  ,fkMaxTheta(1.0)	
  ,fkMaxPhi(2.0) 
  ,fkRoad0y(6.0)
  ,fkRoad0z(8.5) 
  ,fkRoad1y(2.0)
  ,fkRoad1z(20.0)	
  ,fkRoad2y(3.0)
  ,fkRoad2z(20.0)
  ,fkPtThreshold(2.0) 
  ,fkPlaneQualityThreshold(5.0)// 4.2? under Investigation
  ,fkRoadzMultiplicator(1.5)
  ,fkFindable(.333)
  ,fkChi2Z(30./*14.*//*12.5*/)
  ,fkChi2Y(.25)
  ,fkChi2YSlope(7.73)
  ,fkChi2ZSlope(0.069)
  ,fkChi2YCut(0.5)
  ,fkPhiSlope(10.6)
  ,fkNMeanClusters(20.)
  ,fkNSigmaClusters(2.)
  ,fkNClusterNoise(0.)
  ,fkNMeanTracklets(5.5)
  ,fkTrackLikelihood(-15.)
  ,fNumberOfConfigs(3)
  ,fFlags(0)
  ,fRawStreamVersion("DEFAULT")
  ,fMinMaxCutSigma(4.)
  ,fMinLeftRightCutSigma(8.)
  ,fClusMaxThresh(4.5)
  ,fClusSigThresh(3.5)
  ,fTCnexp(1)
  ,fRecEveryNTB(1)
  ,fNumberOfPresamples(0)
  ,fNumberOfPostsamples(0)
{
  //
  // Default constructor
  //
  fSysCovMatrix[0] = 0.; // y direction (1 cm)
  fSysCovMatrix[1] = 0.; // z direction (1 cm)
  fSysCovMatrix[2] = 0.; // snp
  fSysCovMatrix[3] = 0.; // tgl
  fSysCovMatrix[4] = 0.; // 1/pt

  // Xe tail cancellation parameters
  fTCParams[0] = 1.156; // r1
  fTCParams[1] = 0.130; // r2
  fTCParams[2] = 0.114; // c1
  fTCParams[3] = 0.624; // c2
  // Ar tail cancellation parameters
  fTCParams[4] = 6.;    // r1
  fTCParams[5] = 0.62;  // r2
  fTCParams[6] = 0.0087;// c1
  fTCParams[7] = 0.07;  // c2

  memset(fPIDThreshold, 0, AliTRDCalPID::kNMom*sizeof(Double_t));
  memset(fStreamLevel, 0, kTRDreconstructionTasks * sizeof(Int_t));

  SetPIDmethod(AliTRDPIDResponse::kLQ1D);
  SetEightSlices();
  SetImproveTracklets();
  SetLUT();
  SetTailCancelation();
}

//______________________________________________________________
AliTRDrecoParam::AliTRDrecoParam(const AliTRDrecoParam &ref)
  :AliDetectorRecoParam(ref)
  ,fkdNchdy(ref.fkdNchdy)
  ,fkMaxTheta(ref.fkMaxTheta)
  ,fkMaxPhi(ref.fkMaxPhi)
  ,fkRoad0y(ref.fkRoad0y)
  ,fkRoad0z(ref.fkRoad0z)
  ,fkRoad1y(ref.fkRoad1y)
  ,fkRoad1z(ref.fkRoad1z)
  ,fkRoad2y(ref.fkRoad2y)
  ,fkRoad2z(ref.fkRoad2z)
  ,fkPtThreshold(ref.fkPtThreshold)
  ,fkPlaneQualityThreshold(ref.fkPlaneQualityThreshold)
  ,fkRoadzMultiplicator(ref.fkRoadzMultiplicator)
  ,fkFindable(ref.fkFindable)
  ,fkChi2Z(ref.fkChi2Z)
  ,fkChi2Y(ref.fkChi2Y)
  ,fkChi2YSlope(ref.fkChi2YSlope)
  ,fkChi2ZSlope(ref.fkChi2ZSlope)
  ,fkChi2YCut(ref.fkChi2YCut)
  ,fkPhiSlope(ref.fkPhiSlope)
  ,fkNMeanClusters(ref.fkNMeanClusters)
  ,fkNSigmaClusters(ref.fkNSigmaClusters)
  ,fkNClusterNoise(ref.fkNClusterNoise)
  ,fkNMeanTracklets(ref.fkNMeanTracklets)
  ,fkTrackLikelihood(ref.fkTrackLikelihood)
  ,fNumberOfConfigs(ref.fNumberOfConfigs)
  ,fFlags(ref.fFlags)
  ,fRawStreamVersion(ref.fRawStreamVersion)
  ,fMinMaxCutSigma(ref.fMinMaxCutSigma)
  ,fMinLeftRightCutSigma(ref.fMinLeftRightCutSigma)
  ,fClusMaxThresh(ref.fClusMaxThresh)
  ,fClusSigThresh(ref.fClusSigThresh)
  ,fTCnexp(ref.fTCnexp)
  ,fRecEveryNTB(ref.fRecEveryNTB)
  ,fNumberOfPresamples(ref.fNumberOfPresamples)
  ,fNumberOfPostsamples(ref.fNumberOfPostsamples)
{
  //
  // Copy constructor
  //
  memcpy(fSysCovMatrix, ref.fSysCovMatrix, 5*sizeof(Double_t));
  memcpy(fTCParams, ref.fTCParams, 8*sizeof(Double_t));
  memcpy(fPIDThreshold, ref.fPIDThreshold, AliTRDCalPID::kNMom*sizeof(Double_t));
  memcpy(fStreamLevel, ref.fStreamLevel, kTRDreconstructionTasks * sizeof(Int_t));
}

//______________________________________________________________
const AliTRDrecoParam& AliTRDrecoParam::operator=(const AliTRDrecoParam &ref)
{
  //
  // assignment operator
  //

  if(this == &ref) return *this;
  AliDetectorRecoParam::operator=(ref);
  fkdNchdy              = ref.fkdNchdy;
  fkMaxTheta            = ref.fkMaxTheta;
  fkMaxPhi              = ref.fkMaxPhi;
  fkRoad0y              = ref.fkRoad0y;
  fkRoad0z              = ref.fkRoad0z;
  fkRoad1y              = ref.fkRoad1y;
  fkRoad1z              = ref.fkRoad1z;
  fkRoad2y              = ref.fkRoad2y;
  fkRoad2z              = ref.fkRoad2z;
  fkPtThreshold         = ref.fkPtThreshold;
  fkPlaneQualityThreshold= ref.fkPlaneQualityThreshold;
  fkRoadzMultiplicator  = ref.fkRoadzMultiplicator;
  fkFindable            = ref.fkFindable;
  fkChi2Z               = ref.fkChi2Z;
  fkChi2Y               = ref.fkChi2Y;
  fkChi2YSlope          = ref.fkChi2YSlope;
  fkChi2ZSlope          = ref.fkChi2ZSlope;
  fkChi2YCut            = ref.fkChi2YCut;
  fkPhiSlope            = ref.fkPhiSlope;
  fkNMeanClusters       = ref.fkNMeanClusters;
  fkNSigmaClusters      = ref.fkNSigmaClusters;
  fkNClusterNoise       = ref.fkNClusterNoise;
  fkNMeanTracklets      = ref.fkNMeanTracklets;
  fkTrackLikelihood     = ref.fkTrackLikelihood;
  fNumberOfConfigs      = ref.fNumberOfConfigs;
  fFlags                = ref.fFlags;
  fRawStreamVersion     = ref.fRawStreamVersion;
  fMinMaxCutSigma       = ref.fMinMaxCutSigma;
  fMinLeftRightCutSigma = ref.fMinLeftRightCutSigma;
  fClusMaxThresh        = ref.fClusMaxThresh;
  fClusSigThresh        = ref.fClusSigThresh;
  fTCnexp               = ref.fTCnexp;
  fRecEveryNTB          = ref.fRecEveryNTB;
  fNumberOfPresamples   = ref.fNumberOfPresamples;
  fNumberOfPostsamples  = ref.fNumberOfPostsamples;

  memcpy(fSysCovMatrix, ref.fSysCovMatrix, 5*sizeof(Double_t));
  memcpy(fTCParams, ref.fTCParams, 8*sizeof(Double_t));
  memcpy(fPIDThreshold, ref.fPIDThreshold, AliTRDCalPID::kNMom*sizeof(Double_t));
  memcpy(fStreamLevel, ref.fStreamLevel, kTRDreconstructionTasks * sizeof(Int_t));
  return *this;
}

//______________________________________________________________
AliTRDrecoParam *AliTRDrecoParam::GetLowFluxParam()
{
  //
  // Parameters for the low flux environment
  //

  AliTRDrecoParam *rec = new AliTRDrecoParam();
  rec->fkdNchdy = 12.; // pp in TRD
  rec->SetVertexConstrained();
  rec->SetCheckTimeConsistency();
  return rec;

}

//______________________________________________________________
AliTRDrecoParam *AliTRDrecoParam::GetLowFluxHLTParam()
{
  //
  // Parameters for the high flux environment in HLT
  //

  AliTRDrecoParam *rec = GetLowFluxParam();
  rec->fNumberOfConfigs = 2;
  return rec;

}

//______________________________________________________________
AliTRDrecoParam *AliTRDrecoParam::GetHighFluxParam()
{
  //
  // Parameters for the high flux environment
  //

  AliTRDrecoParam *rec = new AliTRDrecoParam();
  rec->fkdNchdy = 4000.; // PbPb in TRD
  rec->SetVertexConstrained();
  rec->SetCheckTimeConsistency();
  return rec;

}

//______________________________________________________________
AliTRDrecoParam *AliTRDrecoParam::GetHighFluxHLTParam()
{
  //
  // Parameters for the high flux environment in HLT
  //

  AliTRDrecoParam *rec = GetHighFluxParam();
  rec->fNumberOfConfigs = 1;
  return rec;

}

//______________________________________________________________
AliTRDrecoParam *AliTRDrecoParam::GetCosmicTestParam()
{
  //
  // Parameters for the cosmics data
  //

  AliTRDrecoParam *par = new AliTRDrecoParam();
  par->fSysCovMatrix[0] = 2.; // y direction (1 cm)
  par->fSysCovMatrix[1] = 2.; // z direction (1 cm)
  par->fkChi2YSlope     = 0.11853;
  par->fkChi2ZSlope     = 0.04527;
  par->fkChi2YCut       = 1.;
  par->fkPhiSlope       = 10.; //3.17954;
  par->fkMaxTheta       = 2.1445;
  par->fkMaxPhi         = 2.7475;
  par->fkNMeanClusters  = 12.89;
  par->fkNSigmaClusters = 2.095;
  par->fkRoadzMultiplicator = 3.;
  par->fStreamLevel[kTracker] = 1;
  par->SetCheckTimeConsistency();
  return par;

}


//______________________________________________________________
Float_t AliTRDrecoParam::GetNClusters() const
{
  // Estimate the number of clusters in the TRD detector
  
  Float_t nclusters = (fkNMeanClusters + 2*fkNSigmaClusters)*fkNMeanTracklets*fkdNchdy;
  nclusters *= 1.+fkNClusterNoise;
  return nclusters;
}

//______________________________________________________________
void AliTRDrecoParam::SetPIDLQslices(Int_t s)
{
// Setting number of slices used by the PID LQ method s={1, 2}
// If PID NN is set this function will change to PID LQ.
 
  if(IsPIDNeuralNetwork()){
    AliWarning("PID set to NN. Changing to LQ.");
    SetPIDNeuralNetwork(kFALSE);
  } 

  switch(s){
  case 1: 
    if(TESTBIT(fFlags, kLQ2D)) CLRBIT(fFlags, kLQ2D);
    break;
  case 2:
    SETBIT(fFlags, kLQ2D);
    break;
  default:
    AliWarning(Form("N[%d] PID LQ slices not implemented. Using default 2.", s));
    SETBIT(fFlags, kLQ2D);
    break;
  }
}

