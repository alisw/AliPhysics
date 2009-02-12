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


#include "AliDetectorRecoParam.h"
#include "AliITSgeomTGeo.h"

class AliITSRecoParam : public AliDetectorRecoParam
{
 public: 
  AliITSRecoParam();
  virtual ~AliITSRecoParam();

  static AliITSRecoParam *GetLowFluxParam();// make reco parameters for low flux env.
  static AliITSRecoParam *GetHighFluxParam();// make reco parameters for high flux env. 
  static AliITSRecoParam *GetCosmicTestParam();// special setting for cosmic  
  static AliITSRecoParam *GetPlaneEffParam(Int_t i);// special setting for Plane Efficiency studies

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

  void PrintParameters() const; 

  void     SetTracker(Int_t tracker=0) { fTracker=tracker; }
  void     SetTrackerDefault() { SetTracker(0); } // = MI and SA
  void     SetTrackerMI() { SetTracker(1); }
  void     SetTrackerV2() { SetTracker(2); }
  Int_t    GetTracker() const { return fTracker; }
  void     SetTrackerSAOnly(Bool_t flag=kTRUE) { fITSonly=flag; }
  Bool_t   GetTrackerSAOnly() const { return fITSonly; }
  void     SetVertexer(Int_t vertexer=0) { fVertexer=vertexer; }
  void     SetVertexer3D() { SetVertexer(0); }
  void     SetVertexerZ() { SetVertexer(1); }
  void     SetVertexerCosmics() { SetVertexer(2); }
  void     SetVertexerIons() { SetVertexer(3); }
  void     SetVertexerSmearMC() { SetVertexer(4); }
  void     SetVertexerFixedOnTDI() {SetVertexer(5);} // for injection tests
  void     SetVertexerFixedOnTED() {SetVertexer(6);} // for injection tests
  Int_t    GetVertexer() const { return fVertexer; }
  void     SetClusterFinder(Int_t cf=0) { fClusterFinder=cf; }
  void     SetClusterFinderV2() { SetClusterFinder(0); }
  void     SetClusterFinderOrig() { SetClusterFinder(1); }
  Int_t    GetClusterFinder() const { return fClusterFinder; }
  void     SetPID(Int_t pid=0) {fPID=pid;}
  void     SetDefaultPID() {SetPID(0);}
  void     SetLandauFitPID() {SetPID(1);}
  Int_t    GetPID() const {return fPID;}

  void     SetVertexer3DFiducialRegions(Float_t dzwid=20.0, Float_t drwid=2.5, Float_t dznar=0.5, Float_t drnar=0.5){
    SetVertexer3DWideFiducialRegion(dzwid,drwid);
    SetVertexer3DNarrowFiducialRegion(dznar,drnar);
  }
  void     SetVertexer3DWideFiducialRegion(Float_t dz=20.0, Float_t dr=2.5){
    fVtxr3DZCutWide=dz; fVtxr3DRCutWide=dr;
  }
  void     SetVertexer3DNarrowFiducialRegion(Float_t dz=0.5, Float_t dr=0.5){
    fVtxr3DZCutNarrow=dz; fVtxr3DRCutNarrow=dr;
  }
  void     SetVertexer3DDeltaPhiCuts(Float_t dphiloose=0.5, Float_t dphitight=0.01){
    fVtxr3DPhiCutLoose=dphiloose;
    fVtxr3DPhiCutTight=dphitight;
  }
  void     SetVertexer3DDCACut(Float_t dca=0.1){
    fVtxr3DDCACut=dca;
  }
  void SetVertexer3DDefaults(){
    SetVertexer3DFiducialRegions();
    SetVertexer3DDeltaPhiCuts();
    SetVertexer3DDCACut();    
  }

  Float_t  GetVertexer3DWideFiducialRegionZ() const {return fVtxr3DZCutWide;}
  Float_t  GetVertexer3DWideFiducialRegionR() const {return fVtxr3DRCutWide;}
  Float_t  GetVertexer3DNarrowFiducialRegionZ() const {return fVtxr3DZCutNarrow;}
  Float_t  GetVertexer3DNarrowFiducialRegionR() const {return fVtxr3DRCutNarrow;}
  Float_t  GetVertexer3DLooseDeltaPhiCut() const {return fVtxr3DPhiCutLoose;}
  Float_t  GetVertexer3DTightDeltaPhiCut() const {return fVtxr3DPhiCutTight;}
  Float_t  GetVertexer3DDCACut() const {return fVtxr3DDCACut;}
  

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
  Double_t GetRoadMisal() const { return fRoadMisal; }
  void     SetRoadMisal(Double_t road=0) { fRoadMisal=road; }

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
  void   SetStepSizeTGeo(Double_t size=0.1) { fStepSizeTGeo=size; return; }
  Double_t GetStepSizeTGeo() const { return fStepSizeTGeo; }
  
  void   SetAllowSharedClusters(Bool_t allow=kTRUE) { fAllowSharedClusters=allow; return; }
  Bool_t GetAllowSharedClusters() const { return fAllowSharedClusters; }

  void   SetClusterErrorsParam(Int_t param=1) { fClusterErrorsParam=param; return; }
  Int_t  GetClusterErrorsParam() const { return fClusterErrorsParam; }
  void   SetClusterMisalErrorY(Float_t e0,Float_t e1,Float_t e2,Float_t e3,Float_t e4,Float_t e5) { fClusterMisalErrorY[0]=e0; fClusterMisalErrorY[1]=e1; fClusterMisalErrorY[2]=e2; fClusterMisalErrorY[3]=e3; fClusterMisalErrorY[4]=e4; fClusterMisalErrorY[5]=e5; return; }
  void   SetClusterMisalErrorZ(Float_t e0,Float_t e1,Float_t e2,Float_t e3,Float_t e4,Float_t e5) { fClusterMisalErrorZ[0]=e0; fClusterMisalErrorZ[1]=e1; fClusterMisalErrorZ[2]=e2; fClusterMisalErrorZ[3]=e3; fClusterMisalErrorZ[4]=e4; fClusterMisalErrorZ[5]=e5; return; }
  void   SetClusterMisalError(Float_t err=0.) { SetClusterMisalErrorY(err,err,err,err,err,err); SetClusterMisalErrorZ(err,err,err,err,err,err); }
  Float_t GetClusterMisalErrorY(Int_t i) const { return fClusterMisalErrorY[i]; }
  Float_t GetClusterMisalErrorZ(Int_t i) const { return fClusterMisalErrorZ[i]; }

  void   SetUseAmplitudeInfo(Bool_t use=kTRUE) { for(Int_t i=0;i<AliITSgeomTGeo::kNLayers;i++) fUseAmplitudeInfo[i]=use; return; }
  void   SetUseAmplitudeInfo(Int_t ilay,Bool_t use) { fUseAmplitudeInfo[ilay]=use; return; }
  Bool_t GetUseAmplitudeInfo(Int_t ilay) const { return fUseAmplitudeInfo[ilay]; }
// Option for Plane Efficiency evaluation
  void   SetComputePlaneEff(Bool_t eff=kTRUE, Bool_t his=kTRUE)
      { fComputePlaneEff=eff; fHistoPlaneEff=his; return; }
  Bool_t GetComputePlaneEff() const { return fComputePlaneEff; }
  Bool_t GetHistoPlaneEff() const { return fHistoPlaneEff; }
  void   SetIPlanePlaneEff(Int_t i=0) {if(i<0 || i>=AliITSgeomTGeo::kNLayers) return; fIPlanePlaneEff=i; }
  Int_t  GetIPlanePlaneEff() const {return fIPlanePlaneEff;}
  void   SetReadPlaneEffFrom0CDB(Bool_t read=kTRUE) { fReadPlaneEffFromOCDB=read; }
  Bool_t GetReadPlaneEffFromOCDB() const { return fReadPlaneEffFromOCDB; }
  void   SetMinPtPlaneEff(Bool_t ptmin=0.) { fMinPtPlaneEff=ptmin; }
  Double_t GetMinPtPlaneEff() const { return fMinPtPlaneEff; }
  void   SetMaxMissingClustersPlaneEff(Int_t max=0) { fMaxMissingClustersPlaneEff=max;}
  Int_t  GetMaxMissingClustersPlaneEff() const {return fMaxMissingClustersPlaneEff;}
  void   SetRequireClusterInOuterLayerPlaneEff(Bool_t out=kTRUE) { fRequireClusterInOuterLayerPlaneEff=out;}
  Bool_t GetRequireClusterInOuterLayerPlaneEff() const {return fRequireClusterInOuterLayerPlaneEff;}
  void   SetRequireClusterInInnerLayerPlaneEff(Bool_t in=kTRUE) { fRequireClusterInInnerLayerPlaneEff=in;}
  Bool_t GetRequireClusterInInnerLayerPlaneEff() const {return fRequireClusterInInnerLayerPlaneEff;}
  void   SetOnlyConstraintPlaneEff(Bool_t con=kFALSE) { fOnlyConstraintPlaneEff=con; }
  Bool_t GetOnlyConstraintPlaneEff() const { return fOnlyConstraintPlaneEff; }
  //
  void   SetExtendedEtaAcceptance(Bool_t ext=kTRUE) { fExtendedEtaAcceptance=ext; return; }
  Bool_t GetExtendedEtaAcceptance() const { return fExtendedEtaAcceptance; }
  void   SetAllowProlongationWithEmptyRoad(Bool_t allow=kTRUE) { fAllowProlongationWithEmptyRoad=allow; return; }  
  Bool_t GetAllowProlongationWithEmptyRoad() const { return fAllowProlongationWithEmptyRoad; }

  void   SetUseBadZonesFromOCDB(Bool_t use=kTRUE) { fUseBadZonesFromOCDB=use; return; }
  Bool_t GetUseBadZonesFromOCDB() const { return fUseBadZonesFromOCDB; }

  void   SetUseSingleBadChannelsFromOCDB(Bool_t use=kTRUE) { fUseSingleBadChannelsFromOCDB=use; return; }
  Bool_t GetUseSingleBadChannelsFromOCDB() const { return fUseSingleBadChannelsFromOCDB; }

  void   SetMinFractionOfBadInRoad(Float_t frac=0) { fMinFractionOfBadInRoad=frac; return; }
  Float_t GetMinFractionOfBadInRoad() const { return fMinFractionOfBadInRoad; }

  void   SetOuterStartLayerSA(Int_t lay) { fOuterStartLayerSA=lay; return; }
  Int_t  GetOuterStartLayerSA() const { return fOuterStartLayerSA; }
  void   SetFactorSAWindowSizes(Double_t fact=1.) { fFactorSAWindowSizes=fact; return; }
  Double_t GetFactorSAWindowSizes() const { return fFactorSAWindowSizes; }

  void SetNLoopsSA(Int_t nl=10) {fNLoopsSA=nl;}
  Int_t GetNLoopsSA() const { return fNLoopsSA;}
  void SetPhiLimitsSA(Double_t phimin,Double_t phimax){
    fMinPhiSA=phimin; fMaxPhiSA=phimax;
  }
  Double_t GetMinPhiSA() const {return fMinPhiSA;}
  Double_t GetMaxPhiSA() const {return fMaxPhiSA;}
  void SetLambdaLimitsSA(Double_t lambmin,Double_t lambmax){
    fMinLambdaSA=lambmin; fMaxLambdaSA=lambmax;
  }
  Double_t GetMinLambdaSA() const {return fMinLambdaSA;}
  Double_t GetMaxLambdaSA() const {return fMaxLambdaSA;}
  
  void   SetSAMinClusterCharge(Float_t minq=0.) {fMinClusterChargeSA=minq;}
  Float_t GetSAMinClusterCharge() const {return fMinClusterChargeSA;}

  void   SetSAOnePointTracks() { fSAOnePointTracks=kTRUE; return; }
  Bool_t GetSAOnePointTracks() const { return fSAOnePointTracks; }

  void   SetSAUseAllClusters() { fSAUseAllClusters=kTRUE; return; }
  Bool_t GetSAUseAllClusters() const { return fSAUseAllClusters; }

  void   SetFindV0s(Bool_t find=kTRUE) { fFindV0s=find; return; }
  Bool_t GetFindV0s() const { return fFindV0s; }

  void   SetLayersParameters();

  void   SetLayerToSkip(Int_t i) { fLayersToSkip[i]=1; return; }
  Int_t  GetLayersToSkip(Int_t i) const { return fLayersToSkip[i]; }

  void   SetUseUnfoldingInClusterFinderSPD(Bool_t use=kTRUE) { fUseUnfoldingInClusterFinderSPD=use; return; }
  Bool_t GetUseUnfoldingInClusterFinderSPD() const { return fUseUnfoldingInClusterFinderSPD; }
  void   SetUseUnfoldingInClusterFinderSDD(Bool_t use=kTRUE) { fUseUnfoldingInClusterFinderSDD=use; return; }
  Bool_t GetUseUnfoldingInClusterFinderSDD() const { return fUseUnfoldingInClusterFinderSDD; }
  void   SetUseUnfoldingInClusterFinderSSD(Bool_t use=kTRUE) { fUseUnfoldingInClusterFinderSSD=use; return; }
  Bool_t GetUseUnfoldingInClusterFinderSSD() const { return fUseUnfoldingInClusterFinderSSD; }

  void   SetUseSDDClusterSizeSelection(Bool_t use=kTRUE) {fUseSDDClusterSizeSelection=use;}
  Bool_t GetUseSDDClusterSizeSelection() const {return fUseSDDClusterSizeSelection;}
  void   SetMinClusterChargeSDD(Float_t qcut=0.){fMinClusterChargeSDD=qcut;}
  Float_t GetMinClusterChargeSDD() const {return fMinClusterChargeSDD;}

  void   SetUseChargeMatchingInClusterFinderSSD(Bool_t use=kTRUE) { fUseChargeMatchingInClusterFinderSSD=use; return; }
  Bool_t GetUseChargeMatchingInClusterFinderSSD() const { return fUseChargeMatchingInClusterFinderSSD; }

  void   SetUseCosmicRunShiftsSSD(Bool_t use=kFALSE) { fUseCosmicRunShiftsSSD=use; return; }
  Bool_t GetUseCosmicRunShiftsSSD() const { return fUseCosmicRunShiftsSSD; }

  // SPD Tracklets (D. Elia)
  void    SetTrackleterOnlyOneTrackletPerC2(Bool_t use= kTRUE) {fTrackleterOnlyOneTrackletPerC2=use; return; }
  Bool_t  GetTrackleterOnlyOneTrackletPerC2() const { return fTrackleterOnlyOneTrackletPerC2; }
  void    SetTrackleterPhiWindow(Float_t w=0.08) {fTrackleterPhiWindow=w;}
  void    SetTrackleterZetaWindow(Float_t w=1.) {fTrackleterZetaWindow=w;}
  Float_t GetTrackleterPhiWindow() const {return fTrackleterPhiWindow;}
  Float_t GetTrackleterZetaWindow() const {return fTrackleterZetaWindow;}
  void    SetTrackleterRemoveClustersFromOverlaps(Bool_t use=kTRUE) { fTrackleterRemoveClustersFromOverlaps=use; return; }
  Bool_t  GetTrackleterRemoveClustersFromOverlaps() const { return fTrackleterRemoveClustersFromOverlaps; }
  void    SetTrackleterPhiOverlapCut(Float_t w=0.005) {fTrackleterPhiOverlapCut=w;}
  void    SetTrackleterZetaOverlapCut(Float_t w=0.05) {fTrackleterZetaOverlapCut=w;}
  Float_t GetTrackleterPhiOverlapCut() const {return fTrackleterPhiOverlapCut;}
  Float_t GetTrackleterZetaOverlapCut() const {return fTrackleterZetaOverlapCut;}

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


  Int_t  fTracker;  // ITS tracker to be used (see AliITSReconstructor)
  Bool_t fITSonly;  // tracking only in ITS (no TPC)
  Int_t  fVertexer; // ITS vertexer to be used (see AliITSReconstructor)
  Int_t  fClusterFinder; // ITS cf to be used (see AliITSReconstructor)
  Int_t  fPID;      // ITS PID method to be used (see AliITSReconstructor)


  Float_t fVtxr3DZCutWide;    // Z extension of the wide fiducial region for vertexer 3D
  Float_t fVtxr3DRCutWide;    // R extension of the wide fiducial region for vertexer 3D
  Float_t fVtxr3DZCutNarrow;  // Z extension of the narrow fiducial region for vertexer 3D
  Float_t fVtxr3DRCutNarrow;  // R extension of the narrow fiducial region for vertexer 3D
  Float_t fVtxr3DPhiCutLoose; // loose deltaPhi cut to define tracklets in vertexer 3D
  Float_t fVtxr3DPhiCutTight; // tight deltaPhi cut to define tracklets in vertexer 3D
  Float_t fVtxr3DDCACut;      // cut on tracklet-to-tracklet DCA in vertexer3D

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

  Double_t fRoadMisal; // [cm] increase of road for misalignment (MI)
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
  Double_t fStepSizeTGeo; // step size (cm)
                     // in AliITStrackerMI::CorrectFor*Material methods
  Bool_t fAllowSharedClusters; // if kFALSE don't set to kITSin tracks with shared clusters (MI)
  Int_t fClusterErrorsParam; // parametrization for cluster errors (MI), see AliITSRecoParam::GetError()
  Float_t fClusterMisalErrorY[AliITSgeomTGeo::kNLayers]; // [cm] additional error on cluster Y pos. due to misalignment (MI,SA)
  Float_t fClusterMisalErrorZ[AliITSgeomTGeo::kNLayers]; // [cm] additional error on cluster Z pos. due to misalignment (MI,SA)

  Bool_t fUseAmplitudeInfo[AliITSgeomTGeo::kNLayers]; // use cluster charge in cluster-track matching (SDD,SSD) (MI)

  // Plane Efficiency evaluation
  Bool_t fComputePlaneEff;  // flag to enable computation of PlaneEfficiency
  Bool_t fHistoPlaneEff;  // flag to enable auxiliary PlaneEff histograms (e.g. residual distributions)
  Int_t  fIPlanePlaneEff; // index of the plane (in the range [0,5])  to study the efficiency
  Bool_t fReadPlaneEffFromOCDB; // enable initial reading of Plane Eff statistics from OCDB
                               // The analized events would be used to increase the statistics
  Double_t fMinPtPlaneEff;  // minimum p_t of the track to be used for Plane Efficiency evaluation
  Int_t  fMaxMissingClustersPlaneEff;  // max n. of (other) layers without a cluster associated to the track
  Bool_t fRequireClusterInOuterLayerPlaneEff; // if kTRUE, then only tracks with an associated cluster on the closest
  Bool_t fRequireClusterInInnerLayerPlaneEff; // outer/inner layer are used. It has no effect for outermost/innermost layer
  Bool_t fOnlyConstraintPlaneEff;  // if kTRUE, use only constrained tracks at primary vertex for Plane Eff.

  Bool_t fExtendedEtaAcceptance;  // enable jumping from TPC to SPD at large eta (MI)
  Bool_t fUseBadZonesFromOCDB; // enable using OCDB info on dead modules and chips (MI)
  Bool_t fUseSingleBadChannelsFromOCDB; // enable using OCDB info on bad single SPD pixels and SDD anodes (MI)
  Float_t fMinFractionOfBadInRoad; // to decide whether to skip the layer (MI)
  Bool_t fAllowProlongationWithEmptyRoad; // allow to prolong even if road is empty (MI)
  Int_t fOuterStartLayerSA;      // outer ITS layer to start track in SA
  Double_t fFactorSAWindowSizes; // larger window sizes in SA
  Int_t fNLoopsSA;               // number of loops in tracker SA
  Double_t fMinPhiSA;               // minimum phi value for SA windows
  Double_t fMaxPhiSA;               // maximum phi value for SA windows
  Double_t fMinLambdaSA;            // minimum lambda value for SA windows
  Double_t fMaxLambdaSA;            // maximum lambda value for SA windows
  Float_t  fMinClusterChargeSA;     // minimum SDD,SSD cluster charge for SA tarcker
  Bool_t fSAOnePointTracks; // one-cluster tracks in SA (only for cosmics!)
  Bool_t fSAUseAllClusters; // do not skip clusters used by MI (same track twice in AliESDEvent!)

  Bool_t fFindV0s;  // flag to enable V0 finder (MI)

  // cluster unfolding in ITS cluster finders
  Bool_t fUseUnfoldingInClusterFinderSPD; // SPD
  Bool_t fUseUnfoldingInClusterFinderSDD; // SDD
  Bool_t fUseUnfoldingInClusterFinderSSD; // SSD

  Bool_t  fUseSDDClusterSizeSelection; // cut on SDD cluster size
  Float_t fMinClusterChargeSDD; // cut on SDD cluster charge

  Bool_t fUseChargeMatchingInClusterFinderSSD; // SSD

  // SPD Tracklets (D. Elia)
  Bool_t  fTrackleterOnlyOneTrackletPerC2;         // Allow only one tracklet per cluster in the outer layer
  Float_t fTrackleterPhiWindow;                    // Search window in phi
  Float_t fTrackleterZetaWindow;                   // Search window in eta
  Bool_t  fTrackleterRemoveClustersFromOverlaps;   // Option to skip clusters in the overlaps
  Float_t fTrackleterPhiOverlapCut;                // Fiducial window in phi for overlap cut
  Float_t fTrackleterZetaOverlapCut;               // Fiducial window in eta for overlap cut
  Bool_t fUseCosmicRunShiftsSSD; // SSD time shifts for cosmic run 2007/2008 (use for data taken up to 18 sept 2008)


  ClassDef(AliITSRecoParam,15) // ITS reco parameters
};

#endif
