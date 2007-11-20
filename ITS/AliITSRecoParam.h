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
#include "AliITSgeomTGeo.h"

class AliITSRecoParam : public TObject
{
 public: 
  AliITSRecoParam();
  virtual ~AliITSRecoParam();

  static AliITSRecoParam *GetLowFluxParam();// make reco parameters for low flux env.
  static AliITSRecoParam *GetHighFluxParam();// make reco parameters for high flux env. 
  static AliITSRecoParam *GetCosmicTestParam();// special setting for cosmic  

  static Int_t GetLayersNotToSkip(Int_t i) { return fgkLayersNotToSkip[i]; }
  static Int_t GetLastLayerToTrackTo() { return fgkLastLayerToTrackTo; }
  static Int_t GetMaxClusterPerLayer() { return fgkMaxClusterPerLayer; }
  static Int_t GetMaxClusterPerLayer5() { return fgkMaxClusterPerLayer5; }
  static Int_t GetMaxClusterPerLayer10() { return fgkMaxClusterPerLayer10; }
  static Int_t GetMaxClusterPerLayer20() { return fgkMaxClusterPerLayer20; }
  static Int_t GetMaxDetectorPerLayer() { return fgkMaxDetectorPerLayer; }
  static Double_t Getriw() { return fgkriw; }
  static Double_t Getdiw() { return fgkdiw; }
  static Double_t GetX0iw() { return fgkX0iw; }
  static Double_t Getrcd() { return fgkrcd; }
  static Double_t Getdcd() { return fgkdcd; }
  static Double_t GetX0cd() { return fgkX0cd; }
  static Double_t Getyr() { return fgkyr; }
  static Double_t Getdr() { return fgkdr; }
  static Double_t Getzm() { return fgkzm; }
  static Double_t Getdm() { return fgkdm; }
  static Double_t Getrs() { return fgkrs; }
  static Double_t Getds() { return fgkds; }
  static Double_t GetrInsideITSscreen() { return fgkrInsideITSscreen; }
  static Double_t GetrInsideSPD1() { return fgkrInsideSPD1; }
  static Double_t GetrPipe() { return fgkrPipe; }
  static Double_t GetrInsidePipe() { return fgkrInsidePipe; }
  static Double_t GetrOutsidePipe() { return fgkrOutsidePipe; }
  static Double_t GetdPipe() { return fgkdPipe; }
  static Double_t GetrInsideShield(Int_t i) { return fgkrInsideShield[i]; }
  static Double_t GetrOutsideShield(Int_t i) { return fgkrOutsideShield[i]; }
  static Double_t Getdshield(Int_t i) { return fgkdshield[i]; }
  static Double_t GetX0shield(Int_t i) { return fgkX0shield[i]; }
  static Double_t GetX0Air() { return fgkX0Air; }
  static Double_t GetX0Be() { return fgkX0Be; }
  static Double_t GetBoundaryWidth() { return fgkBoundaryWidth; }
  static Double_t GetDeltaXNeighbDets() { return fgkDeltaXNeighbDets; }
  static Double_t GetSPDdetzlength() { return fgkSPDdetzlength; }
  static Double_t GetSPDdetxlength() { return fgkSPDdetxlength; }

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



  void   SetUseTGeoInTracker(Int_t use=1) { fUseTGeoInTracker=use; return; }
  Int_t  GetUseTGeoInTracker() const { return fUseTGeoInTracker; }
  
  void   SetAllowSharedClusters(Bool_t allow=kTRUE) { fAllowSharedClusters=allow; return; }
  Bool_t GetAllowSharedClusters() const { return fAllowSharedClusters; }

  void   SetClusterErrorsParam(Int_t param=1) { fClusterErrorsParam=param; return; }
  Int_t  GetClusterErrorsParam() const { return fClusterErrorsParam; }
  void   SetUseAmplitudeInfo(Bool_t use=kTRUE) { for(Int_t i=0;i<AliITSgeomTGeo::kNLayers;i++) fUseAmplitudeInfo[i]=use; return; }
  void   SetUseAmplitudeInfo(Int_t ilay,Bool_t use) { fUseAmplitudeInfo[ilay]=use; return; }
  Bool_t GetUseAmplitudeInfo(Int_t ilay) const { return fUseAmplitudeInfo[ilay]; }
  void   SetExtendedEtaAcceptance(Bool_t ext=kTRUE) { fExtendedEtaAcceptance=ext; return; }
  Bool_t GetExtendedEtaAcceptance() const { return fExtendedEtaAcceptance; }

  void   SetFactorSAWindowSizes(Double_t fact=1.) { fFactorSAWindowSizes=fact; return; }
  Double_t GetFactorSAWindowSizes() const { return fFactorSAWindowSizes; }

  void   SetFindV0s(Bool_t find=kTRUE) { fFindV0s=find; return; }
  Bool_t GetFindV0s() const { return fFindV0s; }

  void   SetLayersParameters();

  void   SetLayerToSkip(Int_t i) { fLayersToSkip[i]=1; return; }
  Int_t  GetLayersToSkip(Int_t i) const { return fLayersToSkip[i]; }

  //

  enum {fgkMaxClusterPerLayer=70000}; //7000*10;   // max clusters per layer
  enum {fgkMaxClusterPerLayer5=28000};//7000*10*2/5;  // max clusters per layer
  enum {fgkMaxClusterPerLayer10=14000};//7000*10*2/10; // max clusters per layer
  enum {fgkMaxClusterPerLayer20=7000};//7000*10*2/20; // max clusters per layer

 protected:
  //
  static const Int_t fgkLayersNotToSkip[AliITSgeomTGeo::kNLayers]; // array with layers not to skip
  static const Int_t fgkLastLayerToTrackTo;  // innermost layer
  static const Int_t fgkMaxDetectorPerLayer; // max clusters per layer
  static const Double_t fgkriw;              // TPC inner wall radius
  static const Double_t fgkdiw;              // TPC inner wall x/X0
  static const Double_t fgkX0iw;             // TPC inner wall X0 
  static const Double_t fgkrcd;              // TPC central drum radius
  static const Double_t fgkdcd;              // TPC central drum x/X0
  static const Double_t fgkX0cd;             // TPC central drum X0
  static const Double_t fgkyr;               // TPC rods y (tracking c.s.)
  static const Double_t fgkdr;               // TPC rods x/X0
  static const Double_t fgkzm;               // TPC membrane z
  static const Double_t fgkdm;               // TPC membrane x/X0
  static const Double_t fgkrs;               // ITS screen radius
  static const Double_t fgkds;               // ITS screed x/X0
  static const Double_t fgkrInsideITSscreen; // inside ITS screen radius
  static const Double_t fgkrInsideSPD1;      // inside SPD1 radius
  static const Double_t fgkrPipe;            // pipe radius
  static const Double_t fgkrInsidePipe;      // inside pipe radius
  static const Double_t fgkrOutsidePipe;     // outside pipe radius
  static const Double_t fgkdPipe;            // pipe x/X0
  static const Double_t fgkrInsideShield[2]; // inside SPD (0) SDD (1) shield radius
  static const Double_t fgkrOutsideShield[2]; // outside SPD (0) SDD (1) shield radius
  static const Double_t fgkdshield[2];        // SPD (0) SDD (1) shield x/X0
  static const Double_t fgkX0shield[2];       // SPD (0) SDD (1) shield X0
  static const Double_t fgkX0Air;             // air X0
  static const Double_t fgkX0Be;              // Berillium X0
  static const Double_t fgkBoundaryWidth;     // to define track at detector boundary
  static const Double_t fgkDeltaXNeighbDets;  // max difference in radius between neighbouring detectors 
  static const Double_t fgkSPDdetzlength;     // SPD ladder length in z
  static const Double_t fgkSPDdetxlength;     // SPD ladder length in x

  Int_t fLayersToSkip[AliITSgeomTGeo::kNLayers]; // array with layers to skip (MI,SA)

  // spatial resolutions of the detectors
  Double_t fSigmaY2[AliITSgeomTGeo::kNLayers]; // y
  Double_t fSigmaZ2[AliITSgeomTGeo::kNLayers]; // z
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
  Double_t fMaxChi2PerCluster[AliITSgeomTGeo::kNLayers-1]; // max chi2 for MIP (MI)
  Double_t fMaxNormChi2NonC[AliITSgeomTGeo::kNLayers]; //max norm chi2 for non constrained tracks (MI)
  Double_t fMaxNormChi2C[AliITSgeomTGeo::kNLayers];  //max norm chi2 for constrained tracks (MI)
  Double_t fMaxNormChi2NonCForHypothesis; //max norm chi2 (on layers 0,1,2) for hypotheis to be kept (MI)
  Double_t fMaxChi2; // used to initialize variables needed to find minimum chi2 (MI,V2)
  Double_t fMaxChi2s[AliITSgeomTGeo::kNLayers];   // max predicted chi2 (cluster & track prol.) (MI)
  //
  Double_t fMaxRoad;   // (V2)
  //
  Double_t fMaxChi2In; // (NOT USED)
  Double_t fMaxChi2sR[AliITSgeomTGeo::kNLayers];  // (NOT USED) 
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


  Int_t fUseTGeoInTracker; // use TGeo to get material budget in tracker MI
  Bool_t fAllowSharedClusters; // if kFALSE don't set to kITSin tracks with shared clusters (MI)
  Int_t fClusterErrorsParam; // parametrization for cluster errors (MI), see AliITSRecoParam::GetError()
  Bool_t fUseAmplitudeInfo[AliITSgeomTGeo::kNLayers]; // use cluster charge in cluster-track matching (SDD,SSD) (MI)
  Bool_t fExtendedEtaAcceptance;  // enable jumping from TPC to SPD at large eta (MI)
  Double_t fFactorSAWindowSizes; // larger window sizes in SA

  Bool_t fFindV0s;  // flag to enable V0 finder (MI)

  ClassDef(AliITSRecoParam,1) // ITS reco parameters
};

#endif
