#ifndef ALIITSRECOPARAM_H
#define ALIITSRECOPARAM_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with ITS reconstruction parameters                                  //
// Origin: andrea.dainese@lnl.infn.it                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TObject.h"

//--------------- move from AliITSrecoV2.h ---------------------------    
const Int_t kMaxLayer = 6;

const Int_t kLayersNotToSkip[6]={0,0,0,0,0,0};
const Int_t kLastLayerToTrackTo=0;

const Int_t kMaxClusterPerLayer=7000*10;
const Int_t kMaxClusterPerLayer5=7000*10*2/5;
const Int_t kMaxClusterPerLayer10=7000*10*2/10;
const Int_t kMaxClusterPerLayer20=7000*10*2/20;
const Int_t kMaxDetectorPerLayer=1000;
//------------- end of move from AliITSrecoV2.h --------------------

const Double_t kriw=80.0,kdiw=0.0053,kX0iw=30.0; // TPC inner wall
const Double_t krcd=61.0,kdcd=0.0053,kX0cd=30.0; // TPC "central drum"
const Double_t kyr=12.8,kdr=0.03; // rods
const Double_t kzm=0.2,kdm=0.40;  // membrane
const Double_t krs=50.0,kds=0.001; // ITS screen
const Double_t krInsideITSscreen=49.0; // inside ITS screen

const Double_t krInsideSPD1=3.7; // inside SPD
const Double_t krPipe=3.; // beam pipe radius
const Double_t krInsidePipe=2.7; // inside beam pipe
const Double_t krOutsidePipe=3.3; // outside beam pipe
const Double_t kdPipe=0.0023; // beam pipe thickness

const Double_t kX0Air=21.82;
const Double_t kX0Be=65.19;
const Double_t kX0shieldSDD=38.6;
const Double_t kX0shieldSPD=42.0;

const Double_t kdshieldSDD=0.0034;
const Double_t krshieldSPD=7.5,kdshieldSPD=0.0097;


const Double_t kBoundaryWidth=0.2; // to define track at detector boundary 
const Double_t kDeltaXNeighbDets=0.5; // max difference in radius between 
                                      // neighbouring detectors

// Size of the SPD sensitive volumes (ladders), for dead zones treatment
const Double_t kSPDdetzlength=6.960; // 7.072-2*0.056
const Double_t kSPDdetxlength=1.298; // 1.410-2*0.056

class AliITSRecoParam : public TObject
{
 public: 
  AliITSRecoParam();
  virtual ~AliITSRecoParam();

  static AliITSRecoParam *GetLowFluxParam();// make reco parameters for low flux env.
  static AliITSRecoParam *GetHighFluxParam();// make reco parameters for high flux env. 
  static AliITSRecoParam *GetCosmicTestParam();// special setting for cosmic  

  
  Double_t GetSigmaY2(Int_t i) const { return fSigmaY2[i]; }
  Double_t GetSigmaZ2(Int_t i) const { return fSigmaZ2[i]; }

  Double_t GetMaxSnp() const { return fMaxSnp; }

  Double_t GetNSigmaYLayerForRoadY() const { return fNSigmaYLayerForRoadY; }
  Double_t GetNSigmaRoadY() const { return fNSigmaRoadY; }
  Double_t GetNSigmaZLayerForRoadZ() const { return fNSigmaZLayerForRoadZ; }
  Double_t GetNSigmaRoadZ() const { return fNSigmaRoadZ; }
  Double_t GetNSigma2RoadYC() const { return fNSigma2RoadYC; }
  Double_t GetNSigma2RoadZC() const { return fNSigma2RoadZC; }
  Double_t GetNSigma2RoadYNonC() const { return fNSigma2RoadYNonC; }
  Double_t GetNSigma2RoadZNonC() const { return fNSigma2RoadZNonC; }

  Double_t GetChi2PerCluster() const { return fChi2PerCluster; }
  Double_t GetMaxChi2PerCluster(Int_t i) const { return fMaxChi2PerCluster[i]; }
  Double_t GetMaxNormChi2NonC(Int_t i) const { return fMaxNormChi2NonC[i]; }
  Double_t GetMaxNormChi2C(Int_t i) const { return fMaxNormChi2C[i]; }
  Double_t GetMaxNormChi2NonCForHypothesis() const { return fMaxNormChi2NonCForHypothesis; }
  Double_t GetMaxChi2() const { return fMaxChi2; }
  Double_t GetMaxChi2s(Int_t i) const { return fMaxChi2s[i]; }
  Double_t GetMaxChi2sR(Int_t i) const { return fMaxChi2sR[i]; }
  Double_t GetMaxChi2In() const { return fMaxChi2In; }
  Double_t GetMaxRoad() const { return fMaxRoad; }
  Double_t GetMaxNormChi2ForGolden(Int_t i) const { return 3.+0.5*i; }

  Double_t GetXVdef() const { return fXV; }
  Double_t GetYVdef() const { return fYV; }
  Double_t GetZVdef() const { return fZV; }
  Double_t GetSigmaXVdef() const { return fSigmaXV; }
  Double_t GetSigmaYVdef() const { return fSigmaYV; }
  Double_t GetSigmaZVdef() const { return fSigmaZV; }

  Double_t GetVertexCut() const { return fVertexCut; }
  Double_t GetMaxDZforPrimTrk() const { return fMaxDZforPrimTrk; }
  Double_t GetMaxDZToUseConstraint() const { return fMaxDZToUseConstraint; }
  Double_t GetMaxDforV0dghtrForProlongation() const { return fMaxDforV0dghtrForProlongation; }
  Double_t GetMaxDForProlongation() const { return fMaxDForProlongation; }
  Double_t GetMaxDZForProlongation() const { return fMaxDZForProlongation; }
  Double_t GetMinPtForProlongation() const { return fMinPtForProlongation; }

  void   SetAddVirtualClustersInDeadZone(Bool_t add=kTRUE) { fAddVirtualClustersInDeadZone=add; return; }  
  Bool_t GetAddVirtualClustersInDeadZone() const { return fAddVirtualClustersInDeadZone; }  
  Double_t GetZWindowDeadZone() const { return fZWindowDeadZone; }
  Double_t GetSigmaXDeadZoneHit2() const { return fSigmaXDeadZoneHit2; }
  Double_t GetSigmaZDeadZoneHit2() const { return fSigmaZDeadZoneHit2; }
  Double_t GetXPassDeadZoneHits() const { return fXPassDeadZoneHits; }



  void   SetUseTGeoInTracker(Bool_t use=kTRUE) { fUseTGeoInTracker=use; return; }
  Bool_t GetUseTGeoInTracker() const { return fUseTGeoInTracker; }
  
  void   SetAllowSharedClusters(Bool_t allow=kTRUE) { fAllowSharedClusters=allow; return; }
  Bool_t GetAllowSharedClusters() const { return fAllowSharedClusters; }

  void   SetUseNominalClusterErrors(Bool_t nominal=kTRUE) { fUseNominalClusterErrors=nominal; return; }
  Bool_t GetUseNominalClusterErrors() const { return fUseNominalClusterErrors; }
  void   SetUseAmplitudeInfo(Bool_t use=kTRUE) { for(Int_t i=0;i<6;i++) fUseAmplitudeInfo[i]=use; return; }
  void   SetUseAmplitudeInfo(Int_t ilay,Bool_t use) { fUseAmplitudeInfo[ilay]=use; return; }
  Bool_t GetUseAmplitudeInfo(Int_t ilay) const { return fUseAmplitudeInfo[ilay]; }


  void   SetFindV0s(Bool_t find=kTRUE) { fFindV0s=find; return; }
  Bool_t GetFindV0s() const { return fFindV0s; }

  void SetLayersParameters();
  //
 protected:
  //
  // spatial resolutions of the detectors
  Double_t fSigmaY2[kMaxLayer]; // y
  Double_t fSigmaZ2[kMaxLayer]; // z
  //
  Double_t fMaxSnp; // maximum of sin(phi)  (MI)
  //
  // search road (MI)
  Double_t fNSigmaYLayerForRoadY; // y
  Double_t fNSigmaRoadY;  // y
  Double_t fNSigmaZLayerForRoadZ; // z
  Double_t fNSigmaRoadZ; // z
  Double_t fNSigma2RoadZC; // z
  Double_t fNSigma2RoadYC; // y
  Double_t fNSigma2RoadZNonC; // z
  Double_t fNSigma2RoadYNonC; // y
  //
  // chi2 cuts
  Double_t fMaxChi2PerCluster[kMaxLayer-1]; // max chi2 for MIP (MI)
  Double_t fMaxNormChi2NonC[kMaxLayer]; //max norm chi2 for non constrained tracks (MI)
  Double_t fMaxNormChi2C[kMaxLayer];  //max norm chi2 for constrained tracks (MI)
  Double_t fMaxNormChi2NonCForHypothesis; //max norm chi2 (on layers 0,1,2) for hypotheis to be kept (MI)
  Double_t fMaxChi2; // used to initialize variables needed to find minimum chi2 (MI,V2)
  Double_t fMaxChi2s[kMaxLayer];   // max predicted chi2 (cluster & track prol.) (MI)
  //
  Double_t fMaxRoad;   // (V2)
  //
  Double_t fMaxChi2In; // (NOT USED)
  Double_t fMaxChi2sR[kMaxLayer];  // (NOT USED) 
  Double_t fChi2PerCluster; // (NOT USED)
  //
  // default primary vertex (MI,V2)
  Double_t fXV;  // x
  Double_t fYV;  // y
  Double_t fZV;  // z
  Double_t fSigmaXV; // x
  Double_t fSigmaYV; // y
  Double_t fSigmaZV; // z
  Double_t fVertexCut; // (V2)
  Double_t fMaxDZforPrimTrk; // maximum (imp. par.)/(1+layer) to define 
                             // a primary and apply vertex constraint (MI)
  Double_t fMaxDZToUseConstraint; // maximum (imp. par.) for tracks to be 
                                  // prolonged with constraint
  // cuts to decide if trying to prolong a TPC track (MI)
  Double_t fMaxDforV0dghtrForProlongation; // max. rphi imp. par. cut for V0 daughter
  //
  Double_t fMaxDForProlongation; // max. rphi imp. par. cut
  Double_t fMaxDZForProlongation; // max. 3D imp. par. cut
  Double_t fMinPtForProlongation; // min. pt cut

  // parameters to create "virtual" clusters in SPD dead zone (MI)
  Bool_t   fAddVirtualClustersInDeadZone; // add if kTRUE
  Double_t fZWindowDeadZone; // window size
  Double_t fSigmaXDeadZoneHit2; // x error virtual cls
  Double_t fSigmaZDeadZoneHit2; // z error virtual cls
  Double_t fXPassDeadZoneHits;  // x distance between clusters


  Bool_t fUseTGeoInTracker; // use TGeo to get material budget in tracker MI
  Bool_t fAllowSharedClusters; // if kFALSE don't set to kITSin tracks with shared clusters (MI)
  Bool_t fUseNominalClusterErrors; // if kFALSE don't modify errors using AliITSClusterParam (MI)
  Bool_t fUseAmplitudeInfo[6]; // use cluster charge in cluster-track matching (SDD,SSD) (MI)

  Bool_t fFindV0s;  // flag to enable V0 finder (MI)

  ClassDef(AliITSRecoParam,1) // ITS reco parameters
};

#endif
