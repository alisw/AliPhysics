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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with ITS reconstruction parameters                                  //
// Origin: andrea.dainese@lnl.infn.it                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



ClassImp(AliITSRecoParam)

const Int_t AliITSRecoParam::fgkLayersNotToSkip[AliITSgeomTGeo::kNLayers]={0,0,0,0,0,0};
const Int_t AliITSRecoParam::fgkLastLayerToTrackTo=0;
const Int_t AliITSRecoParam::fgkMaxDetectorPerLayer=1000;
const Double_t AliITSRecoParam::fgkriw=80.0;
const Double_t AliITSRecoParam::fgkdiw=0.0053;
const Double_t AliITSRecoParam::fgkX0iw=30.0;
const Double_t AliITSRecoParam::fgkrcd=61.0;
const Double_t AliITSRecoParam::fgkdcd=0.0053;
const Double_t AliITSRecoParam::fgkX0cd=30.0;
const Double_t AliITSRecoParam::fgkyr=12.8;
const Double_t AliITSRecoParam::fgkdr=0.03;
const Double_t AliITSRecoParam::fgkzm=0.2;
const Double_t AliITSRecoParam::fgkdm=0.40;
const Double_t AliITSRecoParam::fgkrs=50.0;
const Double_t AliITSRecoParam::fgkds=0.001;
const Double_t AliITSRecoParam::fgkrInsideITSscreen=49.0;
const Double_t AliITSRecoParam::fgkrInsideSPD1=3.7;
const Double_t AliITSRecoParam::fgkrPipe=3.;
const Double_t AliITSRecoParam::fgkrInsidePipe=2.7;
const Double_t AliITSRecoParam::fgkrOutsidePipe=3.3;
const Double_t AliITSRecoParam::fgkdPipe=0.0028;
const Double_t AliITSRecoParam::fgkrInsideShield[2]={7.5,25.0};
const Double_t AliITSRecoParam::fgkrOutsideShield[2]={10.5,30.0};
const Double_t AliITSRecoParam::fgkdshield[2]={0.0097,0.0034};
const Double_t AliITSRecoParam::fgkX0shield[2]={38.6,42.0};
const Double_t AliITSRecoParam::fgkX0Air=21.82;
const Double_t AliITSRecoParam::fgkX0Be=65.19;
const Double_t AliITSRecoParam::fgkBoundaryWidth=0.2;
const Double_t AliITSRecoParam::fgkDeltaXNeighbDets=0.5;
const Double_t AliITSRecoParam::fgkSPDdetzlength=6.960; // 7.072-2*0.056
const Double_t AliITSRecoParam::fgkSPDdetxlength=1.298; // 1.410-2*0.056

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
fUseTGeoInTracker(3),
fStepSizeTGeo(0.1),
fAllowSharedClusters(kTRUE),
fClusterErrorsParam(1),
fComputePlaneEff(kFALSE),
fHistoPlaneEff(kFALSE),
fIPlanePlaneEff(0),
fReadPlaneEffFromOCDB(kFALSE),
fMinPtPlaneEff(0),
fMaxMissingClustersPlaneEff(0),
fRequireClusterInOuterLayerPlaneEff(kFALSE),
fRequireClusterInInnerLayerPlaneEff(kFALSE),
fOnlyConstraintPlaneEff(kFALSE),
fExtendedEtaAcceptance(kFALSE),
fUseBadZonesFromOCDB(kTRUE),
fUseSingleBadChannelsFromOCDB(kFALSE),
fMinFractionOfBadInRoad(0),
fAllowProlongationWithEmptyRoad(kFALSE),
fOuterStartLayerSA(0),
fFactorSAWindowSizes(1.),
fNLoopsSA(33),
fMinPhiSA(0.002),
fMaxPhiSA(0.0145),
fMinLambdaSA(0.003),
fMaxLambdaSA(0.008),
fSAOnePointTracks(kFALSE),
fSAUseAllClusters(kFALSE),
fFindV0s(kTRUE),
fUseUnfoldingInClusterFinderSPD(kFALSE),
fUseUnfoldingInClusterFinderSDD(kTRUE),
fUseUnfoldingInClusterFinderSSD(kTRUE),
fUseChargeMatchingInClusterFinderSSD(kTRUE),
fTrackleterOnlyOneTrackletPerC2(kTRUE),
fTrackleterPhiWindow(0.08),
fTrackleterZetaWindow(1.00),
fTrackleterRemoveClustersFromOverlaps(kFALSE),
fTrackleterPhiOverlapCut(0.005),
fTrackleterZetaOverlapCut(0.05),
fUseCosmicRunShiftsSSD(kFALSE)
{
  //
  // constructor
  //
  SetName("ITS");
  SetTitle("ITS");

  SetLayersParameters();
  for(Int_t i=0; i<AliITSgeomTGeo::kNLayers; i++) fLayersToSkip[i]=0;
  SetUseTGeoInTracker(3);
  SetStepSizeTGeo(0.1);
  SetAllowSharedClusters(kTRUE);
  SetFindV0s(kTRUE);
  SetAddVirtualClustersInDeadZone(kFALSE);
  SetUseAmplitudeInfo(kTRUE);
  SetClusterErrorsParam(1);
  SetClusterMisalError(0.);
  SetVertexer3DDefaults();
}
//_____________________________________________________________________________
AliITSRecoParam::~AliITSRecoParam() 
{
  //
  // destructor
  //  
}
//_____________________________________________________________________________
AliITSRecoParam *AliITSRecoParam::GetHighFluxParam() 
{
  //
  // make default reconstruction  parameters for hig  flux env.
  //
  AliITSRecoParam *param = new AliITSRecoParam();
  param->SetVertexerZ();

  // use of bads from OCDB
  param->SetUseBadZonesFromOCDB(kTRUE);
  param->SetUseSingleBadChannelsFromOCDB(kFALSE);
  // extended eta acceptance
  param->SetExtendedEtaAcceptance(kFALSE);
  // allow to skip layer if no cluster and no bad
  param->SetAllowProlongationWithEmptyRoad(kFALSE);


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
  param->SetNLoopsSA(10);

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

  // find independently ITS SA tracks 
  param->SetSAUseAllClusters();
  param->SetOuterStartLayerSA(AliITSgeomTGeo::GetNLayers()-2);

  // larger seach windows for SA (in case of large misalignments)
  param->SetFactorSAWindowSizes(3.); 

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
  // make special reconstruction parameters for Plane Efficiency study on layer i
  //
  if (i<0 || i>=AliITSgeomTGeo::kNLayers) {
    printf("AliITSRecoParam::GetPlaneEffParam: index of ITS Plane not in the range [0,5]\n");
    printf("returning null pointer");
    return NULL;
  }
  AliITSRecoParam *param;
  param = GetHighFluxParam();
  param->SetComputePlaneEff();
  param->SetLayerToSkip(i);
  param->SetIPlanePlaneEff(i);
  // optimized setting for SPD0 (i==0)
  if (i==0 || i==1) {
    param->fMinPtPlaneEff = 0.200; // high pt particles
    param->fMaxMissingClustersPlaneEff = 1; // at most 1 layer out of 5 without cluster
    param->fRequireClusterInOuterLayerPlaneEff = kTRUE; // cluster on SPD1 must be
    //param->fOnlyConstraintPlaneEff = kTRUE;
  }
  if (i==2 || i==3) {
    param->fMinPtPlaneEff = 0.200; // high pt particles
    param->fMaxMissingClustersPlaneEff = 1; // at most 1 layer out of 5 without cluster
    param->fRequireClusterInOuterLayerPlaneEff = kTRUE;
    //param->fOnlyConstraintPlaneEff = kTRUE;
  }
  if (i==4) {
    param->fMinPtPlaneEff = 0.200; // high pt particles
    param->fMaxMissingClustersPlaneEff = 0; // at most 1 layer out of 5 without cluster
    param->fRequireClusterInOuterLayerPlaneEff = kTRUE;
    //param->fOnlyConstraintPlaneEff = kTRUE;
  }
  if (i==5) {
    param->fMinPtPlaneEff = 0.200; // high pt particles
  }
  //
  return param;
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
