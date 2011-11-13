/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliITSRecoParam.h"
#include "AliESDV0Params.h"
#include "AliLog.h"

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with ITS reconstruction parameters                                  //
// Origin: andrea.dainese@lnl.infn.it                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



ClassImp(AliITSRecoParam)

const Int_t AliITSRecoParam::fgkLayersNotToSkip[AliITSgeomTGeo::kNLayers]={0,0,0,0,0,0};
const Double_t AliITSRecoParam::fgkrInsideShield[2]={7.5,25.0};
const Double_t AliITSRecoParam::fgkrOutsideShield[2]={10.5,30.0};
const Double_t AliITSRecoParam::fgkdshield[2]={0.0097,0.0034};
const Double_t AliITSRecoParam::fgkX0shield[2]={38.6,42.0};

//_____________________________________________________________________________
AliITSRecoParam::AliITSRecoParam() : AliDetectorRecoParam(),
fTracker(0),
fITSonly(kFALSE),
fVertexer(0),
fClusterFinder(0),
fPID(0),
fVtxr3DZCutWide(0.),
fVtxr3DRCutWide(0.),
fVtxr3DZCutNarrow(0.),
fVtxr3DRCutNarrow(0.),
fVtxr3DPhiCutLoose(0.),
fVtxr3DPhiCutTight(0.),
fVtxr3DDCACut(0.),
fVtxr3DPileupAlgo(1),
fVtxr3DHighMultAlgo(1),
fMaxSnp(1.),
fNSigmaYLayerForRoadY(0),
fNSigmaRoadY(0),
fNSigmaZLayerForRoadZ(0),
fNSigmaRoadZ(0),
fNSigma2RoadZC(0),
fNSigma2RoadYC(0),
fNSigma2RoadZNonC(0),
fNSigma2RoadYNonC(0),
fRoadMisal(0),
fMaxNormChi2NonCForHypothesis(0),
fMaxChi2(0),
fMaxRoad(0),
fMaxChi2In(0),
fChi2PerCluster(0),
fSearchForExtras(kTRUE),			     
fXV(0), 
fYV(0),
fZV(0),
fSigmaXV(0),
fSigmaYV(0),
fSigmaZV(0),
fVertexCut(0),
fMaxDZforPrimTrk(0),
fMaxDZToUseConstraint(0), 
fMaxDforV0dghtrForProlongation(0),
fMaxDForProlongation(0),
fMaxDZForProlongation(0),
fMinPtForProlongation(0),
fAddVirtualClustersInDeadZone(kFALSE),
fZWindowDeadZone(0),
fSigmaXDeadZoneHit2(0),
fSigmaZDeadZoneHit2(0),
fXPassDeadZoneHits(0),
fSkipSubdetsNotInTriggerCluster(kTRUE),
fUseTGeoInTracker(3),
fStepSizeTGeo(5.),
fAllowSharedClusters(kTRUE),
fClusterErrorsParam(1),
fComputePlaneEff(kFALSE),
fHistoPlaneEff(kFALSE),
fUseTrackletsPlaneEff(kFALSE),
fMCTrackletsPlaneEff(kFALSE),
fBkgTrackletsPlaneEff(kFALSE),
fTrackleterPhiWindowL1(0.10),
fTrackleterPhiWindowL2(0.07),
fTrackleterZetaWindowL1(0.6),
fTrackleterZetaWindowL2(0.4),
fUpdateOncePerEventPlaneEff(kTRUE),
fMinContVtxPlaneEff(3),
fIPlanePlaneEff(0),
fReadPlaneEffFromOCDB(kFALSE),
fMinPtPlaneEff(0),
fMaxMissingClustersPlaneEff(5),
fMaxMissingClustersOutPlaneEff(5),
fRequireClusterInOuterLayerPlaneEff(kFALSE),
fRequireClusterInInnerLayerPlaneEff(kFALSE),
fOnlyConstraintPlaneEff(kFALSE),
fNSigXFromBoundaryPlaneEff(1.),
fNSigZFromBoundaryPlaneEff(1.),
fImproveWithVertex(kFALSE),
fExtendedEtaAcceptance(kFALSE),
fUseBadZonesFromOCDB(kTRUE),
fUseSingleBadChannelsFromOCDB(kFALSE),
fMinFractionOfBadInRoad(0),
fAllowProlongationWithEmptyRoad(kFALSE),
fInwardFlagSA(kFALSE),
fOuterStartLayerSA(2),
fInnerStartLayerSA(3),
fMinNPointsSA(3),
fFactorSAWindowSizes(1.),
fNLoopsSA(32),
fMinPhiSA(0.002),
fMaxPhiSA(0.0145),
fMinLambdaSA(0.003),
fMaxLambdaSA(0.008),
fMinClusterChargeSA(0.),
fSAOnePointTracks(kFALSE),
fSAUseAllClusters(kFALSE),
fMaxSPDcontrForSAToUseAllClusters(1000000),
fSAUsedEdxInfo(kFALSE),
fSelectBestMIP03(kFALSE),
fFlagFakes(kFALSE),
fUseImproveKalman(kFALSE),
fFindV0s(kTRUE),
fStoreLikeSignV0s(kFALSE),
fUseUnfoldingInClusterFinderSPD(kFALSE),
fUseUnfoldingInClusterFinderSDD(kTRUE),
fUseUnfoldingInClusterFinderSSD(kFALSE),
fUseBadChannelsInClusterFinderSSD(kFALSE),
fUseSDDCorrectionMaps(kTRUE),
fUseSDDClusterSizeSelection(kFALSE),
fMinClusterChargeSDD(0.),
fUseChargeMatchingInClusterFinderSSD(kTRUE),
fTrackleterPhiWindow(0.08),
fTrackleterThetaWindow(0.025),
fTrackleterPhiShift(0.0045),
fTrackleterRemoveClustersFromOverlaps(kFALSE),
fTrackleterPhiOverlapCut(0.005),
fTrackleterZetaOverlapCut(0.05),
fTrackleterPhiRotationAngle(0.0),
fTrackleterNStdDev(1.5),
fScaleDTBySin2T(kFALSE),
fUseCosmicRunShiftsSSD(kFALSE),
fSPDRemoveNoisyFlag(kTRUE),
fSPDRemoveDeadFlag(kTRUE),
fVertexerFastSmearX(0.005),
fVertexerFastSmearY(0.005),
fVertexerFastSmearZ(0.01),
fAlignFilterCosmics(kFALSE),
fAlignFilterCosmicMergeTracks(kTRUE),
fAlignFilterMinITSPoints(4),
fAlignFilterMinITSPointsMerged(4),
fAlignFilterOnlyITSSATracks(kTRUE),
fAlignFilterOnlyITSTPCTracks(kFALSE),
fAlignFilterSkipExtra(kFALSE),
fAlignFilterMaxMatchingAngle(0.085),
fAlignFilterMinAngleWrtModulePlanes(0.52),
fAlignFilterMinPt(0.),
fAlignFilterMaxPt(1.e10),
fAlignFilterFillQANtuples(kTRUE),
//
fMultCutPxDrSPDin(0.1),
fMultCutPxDrSPDout(0.15),
fMultCutPxDz(0.2),
fMultCutDCArz(0.5),
fMultCutMinElectronProbTPC(0.5),
fMultCutMinElectronProbESD(0.1),
fMultCutMinP(0.05),
fMultCutMinRGamma(2.),
fMultCutMinRK0(1.),
fMultCutMinPointAngle(0.98),
fMultCutMaxDCADauther(0.5),
fMultCutMassGamma(0.03),
fMultCutMassGammaNSigma(5.),
fMultCutMassK0(0.03),
fMultCutMassK0NSigma(5.),
fMultCutChi2cGamma(2.),
fMultCutChi2cK0(2.),
fMultCutGammaSFromDecay(-10.),
fMultCutK0SFromDecay(-10.),
fMultCutMaxDCA(1.),
//
fCorrectLorentzAngleSPD(kTRUE),
fTanLorentzAngleHolesSPD(0.017455), // tan(1 degree)
fCorrectLorentzAngleSSD(kTRUE),
fTanLorentzAngleHolesSSD(0.016),  // tan(0.94 degrees)
fTanLorentzAngleElectronsSSD(0.068), // tan(3.98 degrees)
//
fOptReco("All"),
fESDV0Params(NULL)
{
  //
  // constructor
  //
  SetName("ITS");
  SetTitle("ITS");

  SetLayersParameters();
  for(Int_t i=0; i<AliITSgeomTGeo::kNLayers; i++) fLayersToSkip[i]=0;
  SetUseTGeoInTracker(3);
  SetStepSizeTGeo(5.);
  SetAllowSharedClusters(kTRUE);
  SetFindV0s(kTRUE);
  SetAddVirtualClustersInDeadZone(kFALSE);
  SetUseAmplitudeInfo(kTRUE);
  SetClusterErrorsParam(1);
  SetClusterMisalError(0.);
  SetClusterMisalErrorBOn(0.);
  SetVertexer3DDefaults();

  for(Int_t i=0; i<AliITSgeomTGeo::kNLayers; i++) fAlignFilterUseLayer[i]=kTRUE;
  fESDV0Params = new AliESDV0Params();
  
  fESDV0Params->SetMaxDist0(0.1);
  fESDV0Params->SetMaxDist1(0.1);
  fESDV0Params->SetMaxDist(1.);
  fESDV0Params->SetMinPointAngle(0.85);
  fESDV0Params->SetMinPointAngle2(0.99);
  fESDV0Params->SetMinR(0.5);
  fESDV0Params->SetMaxR(220.);
  fESDV0Params->SetMinPABestConst(0.9999);
  fESDV0Params->SetMaxRBestConst(10.);
  fESDV0Params->SetCausality0Cut(0.19);
  fESDV0Params->SetLikelihood01Cut(0.45);
  fESDV0Params->SetLikelihood1Cut(0.5);  
  fESDV0Params->SetCombinedCut(0.55);
  fESDV0Params->SetMinClFullTrk(5.0);
  fESDV0Params->SetMinTgl0(1.05);
  fESDV0Params->SetMinRTgl0(40.0);

  fESDV0Params->SetMinNormDistForbTgl0(3.0);
  fESDV0Params->SetMinClForb0(4.5);
  fESDV0Params->SetMinNormDistForb1(3.0);
  fESDV0Params->SetMinNormDistForb2(2.0);
  fESDV0Params->SetMinNormDistForb3(1.0);
  fESDV0Params->SetMinNormDistForb4(4.0);
  fESDV0Params->SetMinNormDistForb5(5.0);
  fESDV0Params->SetMinNormDistForbProt(2.0);
  fESDV0Params->SetMaxPidProbPionForb(0.5);

  fESDV0Params->SetMinRTPCdensity(40.);
  fESDV0Params->SetMaxRTPCdensity0(110.);
  fESDV0Params->SetMaxRTPCdensity10(120.);
  fESDV0Params->SetMaxRTPCdensity20(130.);
  fESDV0Params->SetMaxRTPCdensity30(140.);

  fESDV0Params->SetMinTPCdensity(0.6);
  fESDV0Params->SetMinTgl1(1.1);
  fESDV0Params->SetMinTgl2(1.0);
  fESDV0Params->SetMinchi2before0(16.);
  fESDV0Params->SetMinchi2before1(16.);
  fESDV0Params->SetMinchi2after0(16.);
  fESDV0Params->SetMinchi2after1(16.);
  fESDV0Params->SetAddchi2SharedCl(18.);
  fESDV0Params->SetAddchi2NegCl0(25.);
  fESDV0Params->SetAddchi2NegCl1(30.);
  fESDV0Params->SetSigp0Par0(0.0001);
  fESDV0Params->SetSigp0Par1(0.001);
  fESDV0Params->SetSigp0Par2(0.1);
  fESDV0Params->SetSigpPar0(0.5);
  fESDV0Params->SetSigpPar1(0.6);
  fESDV0Params->SetSigpPar2(0.4);
  fESDV0Params->SetMaxDcaLh0(0.5);
  fESDV0Params->SetStreamLevel(0);
  fESDV0Params->SetChi2KF(100);
  fESDV0Params->SetRobustChi2KF(100);
  
}
//_____________________________________________________________________________
AliITSRecoParam::~AliITSRecoParam() 
{
  //
  // destructor
  //  
  if(fESDV0Params){
    delete fESDV0Params;
    fESDV0Params=NULL;
  }
}
//_____________________________________________________________________________
AliITSRecoParam *AliITSRecoParam::GetHighFluxParam() 
{
  //
  // make default reconstruction  parameters for hig  flux env.
  //
  AliITSRecoParam *param = new AliITSRecoParam();
  param->SetVertexer3DDefaults();
  param->SetSPDVertexerPileupAlgoOff();
  // use of bads from OCDB
  param->SetUseBadZonesFromOCDB(kTRUE);
  param->SetUseSingleBadChannelsFromOCDB(kFALSE);
  // use pointing to vertex during prolongation
  param->SetImproveWithVertex(kTRUE);
  // extended eta acceptance
  param->SetExtendedEtaAcceptance(kFALSE);
  // allow to skip layer if no cluster and no bad
  param->SetAllowProlongationWithEmptyRoad(kFALSE);
  // set event specie
  param->SetEventSpecie(AliRecoParam::kHighMult);

  param->fMaxSnp = 0.95;

  param->fNSigmaYLayerForRoadY = 4.;
  param->fNSigmaRoadY = 7.5;
  param->fNSigmaZLayerForRoadZ = 4.;
  param->fNSigmaRoadZ = 7.5;

  param->fNSigma2RoadZC = 60.; //7.75^2
  param->fNSigma2RoadYC = 60.; //7.75^2
  param->fNSigma2RoadZNonC = 50.; //7.07^2
  param->fNSigma2RoadYNonC = 50.; //7.07^2

  param->fMaxChi2PerCluster[0] = 11.; //7
  param->fMaxChi2PerCluster[1] = 12.; //5
  param->fMaxChi2PerCluster[2] = 12.; //8
  param->fMaxChi2PerCluster[3] = 5.;  //8
  param->fMaxChi2PerCluster[4] = 12.; //6.5

  param->fMaxNormChi2NonC[0] = 7.;
  param->fMaxNormChi2NonC[1] = 8.;
  param->fMaxNormChi2NonC[2] = 8.;
  param->fMaxNormChi2NonC[3] = 11.;
  param->fMaxNormChi2NonC[4] = 14.;
  param->fMaxNormChi2NonC[5] = 25.;

  param->fMaxNormChi2C[0] = 11.;
  param->fMaxNormChi2C[1] = 13.;
  param->fMaxNormChi2C[2] = 15.;
  param->fMaxNormChi2C[3] = 18.;
  param->fMaxNormChi2C[4] = 30.;
  param->fMaxNormChi2C[5] = 35.;

  param->fMaxNormChi2NonCForHypothesis = 7.;
  
  param->fMaxChi2 = 35.;

  param->fMaxChi2s[0] = 25.; //40   
  param->fMaxChi2s[1] = 25.; //40   
  param->fMaxChi2s[2] = 25.; //40   
  param->fMaxChi2s[3] = 25.; //40   
  param->fMaxChi2s[4] = 40.; //40   
  param->fMaxChi2s[5] = 50.; //40

  param->fMaxRoad = 6.;

  // not used
  param->fMaxChi2In = 16.;
   
  param->fMaxChi2sR[0] = 10.;   
  param->fMaxChi2sR[1] = 10.;   
  param->fMaxChi2sR[2] = 10.;   
  param->fMaxChi2sR[3] = 10.;   
  param->fMaxChi2sR[4] = 30.;   
  param->fMaxChi2sR[5] = 40.;   


  param->fChi2PerCluster = 9.;
  // not used
  param->fSearchForExtras = kFALSE;

  param->fXV = 0.;
  param->fYV = 0.;
  param->fZV = 0.;
  param->fSigmaXV = 0.0050;
  param->fSigmaYV = 0.0050;
  param->fSigmaZV = 0.0100;

  param->fVertexCut = 25.;

  param->fMaxDZforPrimTrk = 0.4;
  param->fMaxDZToUseConstraint = 3.;

  param->fMaxDforV0dghtrForProlongation = 30.;
  param->fMaxDForProlongation = 40.;//10.;
  param->fMaxDZForProlongation = 60.;//20.;
  param->fMinPtForProlongation = 0.120;

  param->fZWindowDeadZone = 2.0;
  param->fSigmaXDeadZoneHit2 = 0.004/12.;
  param->fSigmaZDeadZoneHit2 = 0.001/12.;
  param->fXPassDeadZoneHits = 0.018;
  //
  // mult reco
  param->fMultCutPxDrSPDin = 0.1;
  param->fMultCutPxDrSPDout = 0.15;
  param->fMultCutPxDz = 0.2;
  param->fMultCutDCArz = 0.5;
  param->fMultCutMinElectronProbTPC = 0.5;
  param->fMultCutMinElectronProbESD = 0.1;
  param->fMultCutMinP = 0.05;
  param->fMultCutMinRGamma = 2.;
  param->fMultCutMinRK0 = 1.;
  param->fMultCutMinPointAngle = 0.98;
  param->fMultCutMaxDCADauther = 0.5;
  param->fMultCutMassGamma = 0.03;
  param->fMultCutMassGammaNSigma = 5.;
  param->fMultCutMassK0 = 0.03;
  param->fMultCutMassK0NSigma = 5.;
  param->fMultCutChi2cGamma = 2.;
  param->fMultCutChi2cK0 = 2.;
  param->fMultCutGammaSFromDecay = -10.;
  param->fMultCutK0SFromDecay = -10.;
  param->fMultCutMaxDCA = 1.;  
  //
  // trackleter
  param->fTrackleterPhiWindow = 0.06;
  param->fScaleDTBySin2T = kTRUE;
  //
  param->fSelectBestMIP03 = kFALSE;//kTRUE;
  param->fFlagFakes       = kTRUE;
  param->fUseImproveKalman= kFALSE;
  //
  return param;
}
//_____________________________________________________________________________
AliITSRecoParam *AliITSRecoParam::GetLowFluxParam() 
{
  //
  // make default reconstruction  parameters for low  flux env.
  //
  AliITSRecoParam *param = new AliITSRecoParam();
  param->SetVertexer3DDefaults();

  // full use of bads from OCDB
  param->SetUseBadZonesFromOCDB(kTRUE);
  param->SetUseSingleBadChannelsFromOCDB(kTRUE);
  // extended eta acceptance
  param->SetExtendedEtaAcceptance(kTRUE);
  // allow to skip layer if no cluster and no bad
  param->SetAllowProlongationWithEmptyRoad(kTRUE);
  // set event specie
  param->SetEventSpecie(AliRecoParam::kLowMult);

  param->fMaxSnp = 0.95;

  param->fNSigmaYLayerForRoadY = 4.;
  param->fNSigmaRoadY = 7.5;
  param->fNSigmaZLayerForRoadZ = 4.;
  param->fNSigmaRoadZ = 7.5;

  param->fNSigma2RoadZC = 60.; //7.75^2
  param->fNSigma2RoadYC = 60.; //7.75^2
  param->fNSigma2RoadZNonC = 50.; //7.07^2
  param->fNSigma2RoadYNonC = 50.; //7.07^2

  param->fMaxChi2PerCluster[0] = 11.; //7
  param->fMaxChi2PerCluster[1] = 12.; //5
  param->fMaxChi2PerCluster[2] = 12.; //8
  param->fMaxChi2PerCluster[3] = 5.;  //8
  param->fMaxChi2PerCluster[4] = 12.; //6.5

  param->fMaxNormChi2NonC[0] = 7.;
  param->fMaxNormChi2NonC[1] = 8.;
  param->fMaxNormChi2NonC[2] = 8.;
  param->fMaxNormChi2NonC[3] = 11.;
  param->fMaxNormChi2NonC[4] = 14.;
  param->fMaxNormChi2NonC[5] = 25.;

  param->fMaxNormChi2C[0] = 11.;
  param->fMaxNormChi2C[1] = 13.;
  param->fMaxNormChi2C[2] = 15.;
  param->fMaxNormChi2C[3] = 18.;
  param->fMaxNormChi2C[4] = 30.;
  param->fMaxNormChi2C[5] = 35.;

  param->fMaxNormChi2NonCForHypothesis = 7.;
  
  param->fMaxChi2 = 35.;

  param->fMaxChi2s[0] = 25.; //40   
  param->fMaxChi2s[1] = 25.; //40   
  param->fMaxChi2s[2] = 25.; //40   
  param->fMaxChi2s[3] = 25.; //40   
  param->fMaxChi2s[4] = 40.; //40   
  param->fMaxChi2s[5] = 50.; //40

  param->fMaxRoad = 6.;

  // not used
  param->fMaxChi2In = 16.;
   
  param->fMaxChi2sR[0] = 10.;   
  param->fMaxChi2sR[1] = 10.;   
  param->fMaxChi2sR[2] = 10.;   
  param->fMaxChi2sR[3] = 10.;   
  param->fMaxChi2sR[4] = 30.;   
  param->fMaxChi2sR[5] = 40.;   

  param->fChi2PerCluster = 9.;
  // not used

  param->fXV = 0.;
  param->fYV = 0.;
  param->fZV = 0.;
  param->fSigmaXV = 0.0050;
  param->fSigmaYV = 0.0050;
  param->fSigmaZV = 0.0100;

  param->fVertexCut = 25.;

  param->fMaxDZforPrimTrk = 0.4;
  param->fMaxDZToUseConstraint = 3.;

  param->fMaxDforV0dghtrForProlongation = 30.;
  param->fMaxDForProlongation = 40.;//10.;
  param->fMaxDZForProlongation = 60.;//20.;
  param->fMinPtForProlongation = 0.120;

  param->fZWindowDeadZone = 2.0;
  param->fSigmaXDeadZoneHit2 = 0.004/12.;
  param->fSigmaZDeadZoneHit2 = 0.001/12.;
  param->fXPassDeadZoneHits = 0.018;
  param->SetNLoopsSA(20);
  param->fMaxPhiSA = 0.07;
  param->fMaxLambdaSA = 0.04;

 
  param->GetESDV0Params()->SetMinRTgl0(0.5);
  param->GetESDV0Params()->SetMaxDist(1.5);
  param->GetESDV0Params()->SetMaxDcaLh0(1.5);
  param->GetESDV0Params()->SetMaxRBestConst(80);
  param->GetESDV0Params()->SetMinPABestConst(0.99);
  param->GetESDV0Params()->SetMinNormDistForbTgl0(1.);
  param->GetESDV0Params()->SetMinNormDistForb1(2.);
  param->GetESDV0Params()->SetMinNormDistForbProt(1.);
  param->GetESDV0Params()->SetMaxPidProbPionForb(0.7);
  param->GetESDV0Params()->SetLikelihood01Cut(0.3);
  param->GetESDV0Params()->SetLikelihood1Cut(0.35);
  param->GetESDV0Params()->SetCombinedCut(0.4);

  // mult reco
  param->fMultCutPxDrSPDin = 0.1;
  param->fMultCutPxDrSPDout = 0.15;
  param->fMultCutPxDz = 0.2;
  param->fMultCutDCArz = 0.5;
  param->fMultCutMinElectronProbTPC = 0.5;
  param->fMultCutMinElectronProbESD = 0.1;
  param->fMultCutMinP = 0.05;
  param->fMultCutMinRGamma = 2.;
  param->fMultCutMinRK0 = 1.;
  param->fMultCutMinPointAngle = 0.98;
  param->fMultCutMaxDCADauther = 0.5;
  param->fMultCutMassGamma = 0.03;
  param->fMultCutMassGammaNSigma = 5.;
  param->fMultCutMassK0 = 0.03;
  param->fMultCutMassK0NSigma = 5.;
  param->fMultCutChi2cGamma = 2.;
  param->fMultCutChi2cK0 = 2.;
  param->fMultCutGammaSFromDecay = -10.;
  param->fMultCutK0SFromDecay = -10.;
  param->fMultCutMaxDCA = 1.;  
  //

  return param;
}
//_____________________________________________________________________________
AliITSRecoParam *AliITSRecoParam::GetCosmicTestParam() 
{
  //
  // make default reconstruction  parameters for cosmics
  //
  AliITSRecoParam *param = new AliITSRecoParam();

  // vertexer for cosmics
  param->SetVertexer(2);

  param->SetClusterErrorsParam(2);
  param->SetFindV0s(kFALSE);
  param->SetAddVirtualClustersInDeadZone(kFALSE);
  param->SetUseAmplitudeInfo(kFALSE);

  // set event specie
  param->SetEventSpecie(AliRecoParam::kCosmic);

  // full use of bads from OCDB
  param->SetUseBadZonesFromOCDB(kTRUE);
  param->SetUseSingleBadChannelsFromOCDB(kTRUE);

  // find independently ITS SA tracks 
  param->SetSAUseAllClusters();
  param->SetOuterStartLayerSA(AliITSgeomTGeo::GetNLayers()-2);

  //****** COSMICS 2009 (same as COSMICS 2008) *********************

  // to maximize efficiency
  param->SetAllowProlongationWithEmptyRoad();
  param->SetMinNPointsSA(2);

  // larger seach windows for SA (in case of large misalignments)
  param->SetNLoopsSA(32);
  param->SetFactorSAWindowSizes(20);

  // additional error due to misal (B off)
  param->SetClusterMisalErrorY(1.0,1.0,1.0,1.0,1.0,1.0); // [cm]
  param->SetClusterMisalErrorZ(1.0,1.0,1.0,1.0,1.0,1.0); // [cm]
  // additional error due to misal (B on)
  param->SetClusterMisalErrorYBOn(0.0,0.0,0.1,0.1,0.1,0.1); // [cm]
  param->SetClusterMisalErrorZBOn(0.1,0.1,0.1,0.1,0.1,0.1); // [cm]


  // SDD configuration 
  param->fUseSDDCorrectionMaps = kFALSE;
  param->fUseSDDClusterSizeSelection=kTRUE;
  param->fMinClusterChargeSDD=30.;
  

  // alignment data filter
  param->SetAlignFilterCosmics(kTRUE);
  param->SetAlignFilterCosmicMergeTracks(kTRUE); 
  param->SetAlignFilterMinITSPoints(1);
  param->SetAlignFilterMinITSPointsMerged(3);
  param->SetAlignFilterOnlyITSSATracks(kTRUE);
  param->SetAlignFilterOnlyITSTPCTracks(kFALSE);
  param->SetAlignFilterSkipExtra(kFALSE);
  param->SetAlignFilterMaxMatchingAngle(0.085/*5deg*/);
  param->SetAlignFilterMinPt(0.2);          
  param->SetAlignFilterMaxPt(1.e10);          
  param->SetAlignFilterFillQANtuples(kTRUE);    

  //******************************************************************

  param->fMaxSnp = 0.95;

  param->fNSigmaYLayerForRoadY = 4.;
  param->fNSigmaRoadY = 7.5;
  param->fNSigmaZLayerForRoadZ = 4.;
  param->fNSigmaRoadZ = 7.5;

  param->fNSigma2RoadZC = 60.; //7.75^2
  param->fNSigma2RoadYC = 60.; //7.75^2
  param->fNSigma2RoadZNonC = 50.; //7.07^2
  param->fNSigma2RoadYNonC = 50.; //7.07^2

  param->fMaxChi2PerCluster[0] = 11.; //7
  param->fMaxChi2PerCluster[1] = 12.; //5
  param->fMaxChi2PerCluster[2] = 12.; //8
  param->fMaxChi2PerCluster[3] = 5.;  //8
  param->fMaxChi2PerCluster[4] = 12.; //6.5

  param->fMaxNormChi2NonC[0] = 7.;
  param->fMaxNormChi2NonC[1] = 8.;
  param->fMaxNormChi2NonC[2] = 8.;
  param->fMaxNormChi2NonC[3] = 11.;
  param->fMaxNormChi2NonC[4] = 14.;
  param->fMaxNormChi2NonC[5] = 25.;

  param->fMaxNormChi2C[0] = 11.;
  param->fMaxNormChi2C[1] = 13.;
  param->fMaxNormChi2C[2] = 15.;
  param->fMaxNormChi2C[3] = 18.;
  param->fMaxNormChi2C[4] = 30.;
  param->fMaxNormChi2C[5] = 35.;

  param->fMaxNormChi2NonCForHypothesis = 7.;
  
  param->fMaxChi2 = 35.;

  param->fMaxChi2s[0] = 25.; //40   
  param->fMaxChi2s[1] = 25.; //40   
  param->fMaxChi2s[2] = 25.; //40   
  param->fMaxChi2s[3] = 25.; //40   
  param->fMaxChi2s[4] = 40.; //40   
  param->fMaxChi2s[5] = 50.; //40

  param->fMaxRoad = 6.;

  // not used
  param->fMaxChi2In = 16.;
   
  param->fMaxChi2sR[0] = 10.;   
  param->fMaxChi2sR[1] = 10.;   
  param->fMaxChi2sR[2] = 10.;   
  param->fMaxChi2sR[3] = 10.;   
  param->fMaxChi2sR[4] = 30.;   
  param->fMaxChi2sR[5] = 40.;   

  param->fChi2PerCluster = 9.;
  // not used

  param->fXV = 0.;
  param->fYV = 0.;
  param->fZV = 0.;
  param->fSigmaXV = 0.0050;
  param->fSigmaYV = 0.0050;
  param->fSigmaZV = 0.0100;

  param->fVertexCut = 25.;

  param->fMaxDZforPrimTrk = 0.4;
  param->fMaxDZToUseConstraint = 3.;

  param->fMaxDforV0dghtrForProlongation = 30.;
  param->fMaxDForProlongation = 10.;
  param->fMaxDZForProlongation = 20.;
  param->fMinPtForProlongation = 0.120;

  param->fZWindowDeadZone = 2.0;
  param->fSigmaXDeadZoneHit2 = 0.004/12.;
  param->fSigmaZDeadZoneHit2 = 0.001/12.;
  param->fXPassDeadZoneHits = 0.018;

  return param;
}
//_____________________________________________________________________________
AliITSRecoParam *AliITSRecoParam::GetPlaneEffParam(Int_t i)
{
  //
  // make special reconstruction parameters for Plane Efficiency study on layer i (0,5)
  // 
  // if i=-1, then the evaluation for both pixel layers is tried with the tracklet method
  //
 if (i<-1 || i>=AliITSgeomTGeo::kNLayers) {
    printf("AliITSRecoParam::GetPlaneEffParam: index of ITS Plane nor in the range [0,5] neither =-1\n");
    printf("returning null pointer");
    return NULL;
  }
  if(i>=0) {  // Method using tracks (remove given plane from tracking)
    AliITSRecoParam *param;
    param = GetLowFluxParam();
    param->SetClusterErrorsParam(2);
    // find independently ITS SA tracks 
    param->SetSAUseAllClusters();
    param->SetOuterStartLayerSA(2);
    param->SetAllowProlongationWithEmptyRoad(kTRUE);
    // larger seach windows for SA (in case of large misalignments)
    param->SetFactorSAWindowSizes(2);

    // Misalignment syst errors decided at ITS meeting 25.03.2010
    // additional error due to misal (B off)
    param->SetClusterMisalErrorY(0.0010,0.0010,0.0300,0.0300,0.0020,0.0020); // [cm]
    param->SetClusterMisalErrorZ(0.0100,0.0100,0.0100,0.0100,0.0500,0.0500); // [cm]
    // additional error due to misal (B on)
    param->SetClusterMisalErrorYBOn(0.0010,0.0030,0.0500,0.0500,0.0020,0.0020); // [cm]
    param->SetClusterMisalErrorZBOn(0.0100,0.0100,0.0100,0.0100,0.0500,0.0500); // [cm]
    //----

    // SDD configuration 
    param->SetUseSDDCorrectionMaps(kTRUE); // changed 30.04.2010
    param->SetUseSDDClusterSizeSelection(kTRUE);
    param->SetMinClusterChargeSDD(30.);
    param->SetUseUnfoldingInClusterFinderSDD(kFALSE);
    //
    param->SetComputePlaneEff();
    param->SetLayerToSkip(i);
    param->SetIPlanePlaneEff(i);
    param->fNSigXFromBoundaryPlaneEff= 1.;
    param->fNSigZFromBoundaryPlaneEff= 1.;
    // optimized setting for SPD0 (i==0)
    if (i==0) {
      param->fMinPtPlaneEff = 0.200; // high pt particles
      param->fMaxMissingClustersPlaneEff = 2; // at most 2 layers out of 5 without cluster
      param->fMaxMissingClustersOutPlaneEff = 2; // at most 2 layers out of 5 external ones without cluster
      param->fRequireClusterInOuterLayerPlaneEff = kTRUE; // cluster on SPD1 must be
      //param->fOnlyConstraintPlaneEff = kTRUE;
    }
    if (i==1 ) {
      param->fMinPtPlaneEff = 0.200; // high pt particles
      param->fMaxMissingClustersPlaneEff = 2; // at most 2 layer out of 5 without cluster
      param->fMaxMissingClustersOutPlaneEff = 2; // at most 2 layer out of 4 external ones without cluster
      //param->fRequireClusterInOuterLayerPlaneEff = kTRUE; // cluster on SSD1 must be
    }
    if (i==2) {
      param->fMinPtPlaneEff = 0.200; // high pt particles
      param->fMaxMissingClustersPlaneEff = 2; // at most 2 layer out of 5 without cluster
      param->fMaxMissingClustersOutPlaneEff = 2; // at most 2 layer out of 3 external ones without cluster
      //param->fRequireClusterInOuterLayerPlaneEff = kTRUE;
      //param->fOnlyConstraintPlaneEff = kTRUE;
    }
    if (i==3) {
      param->fMinPtPlaneEff = 0.200; // high pt particles
      param->fMaxMissingClustersPlaneEff = 2; // at most 2 layer out of 5 without cluster
      param->fMaxMissingClustersOutPlaneEff = 1; // at most 1 layer out of 2 external ones without cluster
    }
    if (i==4) {
      param->fMinPtPlaneEff = 0.200; // high pt particles
      param->fMaxMissingClustersPlaneEff = 2; // at most 2 layer out of 5 without cluster
      param->fMaxMissingClustersOutPlaneEff = 1; // at most 1 layer out of 1 external ones without cluster
      // param->fRequireClusterInOuterLayerPlaneEff = kTRUE;
      //param->fOnlyConstraintPlaneEff = kTRUE;
    }
    if (i==5) {
      param->fMinPtPlaneEff = 0.200; // high pt particles
      param->fMaxMissingClustersPlaneEff = 2; // at most 2 layer out of 5 without cluster
    }
    //
    return param;
  }
  else if (i==-1) { // Method using tracklets
    AliITSRecoParam *param;
    param = GetLowFluxParam();
    param->SetIPlanePlaneEff(i);
    param->SetComputePlaneEff(kTRUE,kFALSE);
    param->SetUseTrackletsPlaneEff(kTRUE);
    param->SetTrackleterPhiWindowL2(0.07);
    param->SetTrackleterZetaWindowL2(0.4);
    param->SetTrackleterPhiWindowL1(0.10);
    param->SetTrackleterZetaWindowL1(0.6);
    param->SetUpdateOncePerEventPlaneEff(kTRUE);
    param->SetMinContVtxPlaneEff(3);
    return param;
  }
  else {
    AliErrorGeneral("AliITSRecoParam",Form("Unrecognised value of i %d\n",i));
    return 0;
  }
}

//_____________________________________________________________________________
void AliITSRecoParam::SetLayersParameters() 
{
  //
  // number of layers and layers spatial resolutions
  //

  // spatial resolutions of the detectors
  // y: 12 12 38 38 20 20 micron
  fSigmaY2[0]=1.44e-6;
  fSigmaY2[1]=1.44e-6;
  fSigmaY2[2]=1.444e-5;
  fSigmaY2[3]=1.444e-5;
  fSigmaY2[4]=4.0e-6;
  fSigmaY2[5]=4.0e-6;
  // z: 120 120 28 28 830 830 micron
  fSigmaZ2[0]=1.44e-4;
  fSigmaZ2[1]=1.44e-4;
  fSigmaZ2[2]=7.84e-6;
  fSigmaZ2[3]=7.84e-6;
  fSigmaZ2[4]=6.889e-3;
  fSigmaZ2[5]=6.889e-3;

  return;
}
//_____________________________________________________________________________
void AliITSRecoParam::PrintParameters() const 
{
  //
  // print parameters
  //

  printf("=============================  AliITSRecoParam::PrintParameters ");
  printf("============================= \n\n");
  for(Int_t i=0; i<AliITSgeomTGeo::kNLayers; i++) {
    if(!fLayersToSkip[i]) {
      printf("ITS Traking: using layer %d\n",i);
    } else {
      printf("ITS Traking: skipping layer %d\n",i);
    }
  }
  for(Int_t i=0; i<AliITSgeomTGeo::kNLayers; i++) {
    if(fUseAmplitudeInfo[i]) {
      printf("ITS Traking: use amplitude info for layer %d\n",i);
    } else {
      printf("ITS Traking: don't use amplitude info for layer %d\n",i);
    }
  }
  for(Int_t i=0; i<AliITSgeomTGeo::kNLayers; i++)
    printf("Layer %d:\n  sigmaY2 %f, sigmaZ2 %f\n  sigmaMisalY %f, sigmaMisalZ %f\n  max norm chi2 for non constrained tracks %f\n  max norm chi2 for constrained tracks %f\n  max predicted chi2 (cluster & track prol.) %f\n",i,fSigmaY2[i],fSigmaZ2[i],fClusterMisalErrorY[i],fClusterMisalErrorZ[i],fMaxNormChi2NonC[i],fMaxNormChi2C[i],fMaxChi2s[i]);


  Dump();

  return;
}


//_____________________________________________________________________________
Bool_t AliITSRecoParam::SetOptReco(TString r){
  // Set option for local reconstruction. 
  // The string must contain at least one of the following
  // substrings: "All", "SPD", "SDD", "SSD"
  Bool_t isFine = kFALSE;
  if(r.Contains("All") || r.Contains("SPD") || r.Contains("SDD") 
     || r.Contains("SSD")){
      isFine = kTRUE;
      fOptReco=r;
  }
  return isFine;
} 

