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
//#include "AliESDV0Params.h"

class AliESDV0Params;

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
  static Int_t GetMaxClusterPerLayer() { return kMaxClusterPerLayer; }
  static Int_t GetMaxClusterPerLayer5() { return kMaxClusterPerLayer5; }
  static Int_t GetMaxClusterPerLayer10() { return kMaxClusterPerLayer10; }
  static Int_t GetMaxClusterPerLayer20() { return kMaxClusterPerLayer20; }
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
  void     SetVertexerSmearMC(Float_t smearx=0.005, Float_t smeary=0.005, Float_t smearz=0.01) { 
    fVertexerFastSmearX=smearx;  fVertexerFastSmearY=smeary; fVertexerFastSmearZ=smearz; SetVertexer(4); 
  }
  void     SetVertexerFixedOnTDI() {SetVertexer(5);} // for injection tests
  void     SetVertexerFixedOnTED() {SetVertexer(6);} // for injection tests
  Int_t    GetVertexer() const { return fVertexer; }
  Float_t  GetVertexerFastSmearX() const {return fVertexerFastSmearX;}
  Float_t  GetVertexerFastSmearY() const {return fVertexerFastSmearY;}
  Float_t  GetVertexerFastSmearZ() const {return fVertexerFastSmearZ;}

  void     SetClusterFinder(Int_t cf=0) { fClusterFinder=cf; }
  void     SetClusterFinderV2() { SetClusterFinder(0); }
  void     SetClusterFinderOrig() { SetClusterFinder(1); }
  Int_t    GetClusterFinder() const { return fClusterFinder; }
  void     SetPID(Int_t pid=0) {fPID=pid;}
  void     SetDefaultPID() {SetPID(0);}
  void     SetLandauFitPID() {SetPID(1);}
  Int_t    GetPID() const {return fPID;}

  void     SetVertexer3DFiducialRegions(Float_t dzwid=40.0, Float_t drwid=2.5, Float_t dznar=0.5, Float_t drnar=0.5){
    SetVertexer3DWideFiducialRegion(dzwid,drwid);
    SetVertexer3DNarrowFiducialRegion(dznar,drnar);
  }
  void     SetVertexer3DWideFiducialRegion(Float_t dz=40.0, Float_t dr=2.5){
    fVtxr3DZCutWide=dz; fVtxr3DRCutWide=dr;
  }
  void     SetVertexer3DNarrowFiducialRegion(Float_t dz=0.5, Float_t dr=0.5){
    fVtxr3DZCutNarrow=dz; fVtxr3DRCutNarrow=dr;
  }
  void     SetVertexer3DDeltaPhiCuts(Float_t dphiloose=0.5, Float_t dphitight=0.025){
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
  void SetSPDVertexerPileupAlgoOff(){fVtxr3DPileupAlgo=3;}
  void SetSPDVertexerPileupAlgoZ(){fVtxr3DPileupAlgo=0;}
  void SetSPDVertexerPileupAlgo3DTwoSteps(){fVtxr3DPileupAlgo=1;}
  void SetSPDVertexerPileupAlgo3DOneShot(){fVtxr3DPileupAlgo=2;}
  void SetSPDVertexerHighMultAlgoDownscale(){fVtxr3DHighMultAlgo=0;}
  void SetSPDVertexerHighMultAlgoTraces(){fVtxr3DHighMultAlgo=1;}
  //
  Bool_t   GetSelectBestMIP03()                 const {return fSelectBestMIP03;}
  Bool_t   GetFlagFakes()                       const {return fFlagFakes;}
  Bool_t   GetUseImproveKalman()                const {return fUseImproveKalman;}
  void     SetSelectBestMIP03(Bool_t v=kTRUE)         {fSelectBestMIP03 = v;}
  void     SetFlagFakes(Bool_t v=kTRUE)               {fFlagFakes = v;}
  void     SetUseImproveKalman(Bool_t v=kTRUE)        {fUseImproveKalman = v;}
  //
  Float_t  GetVertexer3DWideFiducialRegionZ() const {return fVtxr3DZCutWide;}
  Float_t  GetVertexer3DWideFiducialRegionR() const {return fVtxr3DRCutWide;}
  Float_t  GetVertexer3DNarrowFiducialRegionZ() const {return fVtxr3DZCutNarrow;}
  Float_t  GetVertexer3DNarrowFiducialRegionR() const {return fVtxr3DRCutNarrow;}
  Float_t  GetVertexer3DLooseDeltaPhiCut() const {return fVtxr3DPhiCutLoose;}
  Float_t  GetVertexer3DTightDeltaPhiCut() const {return fVtxr3DPhiCutTight;}
  Float_t  GetVertexer3DDCACut() const {return fVtxr3DDCACut;}
  Int_t    GetSPDVertexerPileupAlgo() const {return fVtxr3DPileupAlgo;}
  UChar_t  GetSPDVertexerHighMultAlgo() const {return fVtxr3DHighMultAlgo;}

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

  void     SetSearchForExtraClusters(Bool_t opt=kTRUE){ fSearchForExtras=opt; }
  Double_t GetSearchForExtraClusters() const { return fSearchForExtras; }

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

  Bool_t   GetSkipSubdetsNotInTriggerCluster() const { return fSkipSubdetsNotInTriggerCluster; }
  void     SetSkipSubdetsNotInTriggerCluster(Bool_t flag=kTRUE) { fSkipSubdetsNotInTriggerCluster=flag; }

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
  void   SetClusterMisalErrorYBOn(Float_t e0,Float_t e1,Float_t e2,Float_t e3,Float_t e4,Float_t e5) { fClusterMisalErrorYBOn[0]=e0; fClusterMisalErrorYBOn[1]=e1; fClusterMisalErrorYBOn[2]=e2; fClusterMisalErrorYBOn[3]=e3; fClusterMisalErrorYBOn[4]=e4; fClusterMisalErrorYBOn[5]=e5; return; }
  void   SetClusterMisalErrorZBOn(Float_t e0,Float_t e1,Float_t e2,Float_t e3,Float_t e4,Float_t e5) { fClusterMisalErrorZBOn[0]=e0; fClusterMisalErrorZBOn[1]=e1; fClusterMisalErrorZBOn[2]=e2; fClusterMisalErrorZBOn[3]=e3; fClusterMisalErrorZBOn[4]=e4; fClusterMisalErrorZBOn[5]=e5; return; }
  void   SetClusterMisalErrorBOn(Float_t err=0.) { SetClusterMisalErrorYBOn(err,err,err,err,err,err); SetClusterMisalErrorZBOn(err,err,err,err,err,err); }
  Float_t GetClusterMisalErrorY(Int_t i,Double_t b=0.) const { return (TMath::Abs(b)<0.0001 ? fClusterMisalErrorY[i] : fClusterMisalErrorYBOn[i]); }
  Float_t GetClusterMisalErrorZ(Int_t i,Double_t b=0.) const { return (TMath::Abs(b)<0.0001 ? fClusterMisalErrorZ[i] : fClusterMisalErrorZBOn[i]); }

  void   SetUseAmplitudeInfo(Bool_t use=kTRUE) { for(Int_t i=0;i<AliITSgeomTGeo::kNLayers;i++) fUseAmplitudeInfo[i]=use; return; }
  void   SetUseAmplitudeInfo(Int_t ilay,Bool_t use) { fUseAmplitudeInfo[ilay]=use; return; }
  Bool_t GetUseAmplitudeInfo(Int_t ilay) const { return fUseAmplitudeInfo[ilay]; }
// Option for Plane Efficiency evaluation
  void   SetComputePlaneEff(Bool_t eff=kTRUE, Bool_t his=kTRUE)
      { fComputePlaneEff=eff; fHistoPlaneEff=his; return; }
  Bool_t GetComputePlaneEff() const { return fComputePlaneEff; }
  Bool_t GetHistoPlaneEff() const { return fHistoPlaneEff; }
  void    SetUseTrackletsPlaneEff(Bool_t use=kTRUE) {fUseTrackletsPlaneEff=use; return;}
  Bool_t  GetUseTrackletsPlaneEff() const {return fUseTrackletsPlaneEff;}
  void    SetOptTrackletsPlaneEff(Bool_t mc=kFALSE,Bool_t bkg=kFALSE)
           {fMCTrackletsPlaneEff=mc;fBkgTrackletsPlaneEff=bkg; return;}
  Bool_t  GetMCTrackletsPlaneEff() const {return fMCTrackletsPlaneEff;}
  Bool_t  GetBkgTrackletsPlaneEff() const {return fBkgTrackletsPlaneEff;}
  void    SetTrackleterPhiWindowL1(Float_t w=0.10) {fTrackleterPhiWindowL1=w; return;}
  Float_t GetTrackleterPhiWindowL1() const {return fTrackleterPhiWindowL1;}
  void    SetTrackleterPhiWindowL2(Float_t w=0.07) {fTrackleterPhiWindowL2=w; return;}
  Float_t GetTrackleterPhiWindowL2() const {return fTrackleterPhiWindowL2;}
  void    SetTrackleterZetaWindowL1(Float_t w=0.6) {fTrackleterZetaWindowL1=w; return;}
  Float_t GetTrackleterZetaWindowL1() const {return fTrackleterZetaWindowL1;}
  void    SetTrackleterZetaWindowL2(Float_t w=0.40) {fTrackleterZetaWindowL2=w; return;}
  Float_t GetTrackleterZetaWindowL2() const {return fTrackleterZetaWindowL2;}
  void    SetUpdateOncePerEventPlaneEff(Bool_t use=kTRUE) {fUpdateOncePerEventPlaneEff=use; return;}
  Bool_t  GetUpdateOncePerEventPlaneEff() const {return fUpdateOncePerEventPlaneEff;}
  void    SetMinContVtxPlaneEff(Int_t n=3) {fMinContVtxPlaneEff=n; return;}
  Int_t   GetMinContVtxPlaneEff() const {return fMinContVtxPlaneEff;}
  void   SetIPlanePlaneEff(Int_t i=0) {if(i<-1 || i>=AliITSgeomTGeo::kNLayers) return; fIPlanePlaneEff=i; }
  Int_t  GetIPlanePlaneEff() const {return fIPlanePlaneEff;}
  void   SetReadPlaneEffFrom0CDB(Bool_t read=kTRUE) { fReadPlaneEffFromOCDB=read; }
  Bool_t GetReadPlaneEffFromOCDB() const { return fReadPlaneEffFromOCDB; }
  void   SetMinPtPlaneEff(Bool_t ptmin=0.) { fMinPtPlaneEff=ptmin; }
  Double_t GetMinPtPlaneEff() const { return fMinPtPlaneEff; }
  void   SetMaxMissingClustersPlaneEff(Int_t max=0) { fMaxMissingClustersPlaneEff=max;}
  Int_t  GetMaxMissingClustersPlaneEff() const {return fMaxMissingClustersPlaneEff;}
  void   SetMaxMissingClustersOutPlaneEff(Int_t max=0) { fMaxMissingClustersOutPlaneEff=max;}
  Int_t  GetMaxMissingClustersOutPlaneEff() const {return fMaxMissingClustersOutPlaneEff;}
  void   SetRequireClusterInOuterLayerPlaneEff(Bool_t out=kTRUE) { fRequireClusterInOuterLayerPlaneEff=out;}
  Bool_t GetRequireClusterInOuterLayerPlaneEff() const {return fRequireClusterInOuterLayerPlaneEff;}
  void   SetRequireClusterInInnerLayerPlaneEff(Bool_t in=kTRUE) { fRequireClusterInInnerLayerPlaneEff=in;}
  Bool_t GetRequireClusterInInnerLayerPlaneEff() const {return fRequireClusterInInnerLayerPlaneEff;}
  void   SetOnlyConstraintPlaneEff(Bool_t con=kFALSE) { fOnlyConstraintPlaneEff=con; }
  Bool_t GetOnlyConstraintPlaneEff() const { return fOnlyConstraintPlaneEff; }
  void SetNSigXFromBoundaryPlaneEff(Double_t nsigx=1.) {fNSigXFromBoundaryPlaneEff=nsigx;}
  Double_t GetNSigXFromBoundaryPlaneEff() const {return fNSigXFromBoundaryPlaneEff;}
  void SetNSigZFromBoundaryPlaneEff(Double_t nsigz=1.) {fNSigZFromBoundaryPlaneEff=nsigz;}
  Double_t GetNSigZFromBoundaryPlaneEff() const {return fNSigZFromBoundaryPlaneEff;}
  //
  void   SetImproveWithVertex(Bool_t impr=kFALSE) { fImproveWithVertex=impr; return; }
  Bool_t GetImproveWithVertex() const { return fImproveWithVertex; }
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

  void   SetOutwardFindingSA() {fInwardFlagSA=kFALSE;}
  void   SetInwardFindingSA() {fInwardFlagSA=kTRUE;}
  Bool_t GetInwardFindingSA() const {return fInwardFlagSA;}
  void   SetOuterStartLayerSA(Int_t lay) { fOuterStartLayerSA=lay; return; }
  Int_t  GetOuterStartLayerSA() const { return fOuterStartLayerSA; }
  void   SetInnerStartLayerSA(Int_t lay) { fInnerStartLayerSA=lay; return; }
  Int_t  GetInnerStartLayerSA() const { return fInnerStartLayerSA; }
  void   SetMinNPointsSA(Int_t np) { fMinNPointsSA=np; return; }
  Int_t  GetMinNPointsSA() const { return fMinNPointsSA;}
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

  void   SetSAUseAllClusters(Bool_t opt=kTRUE) { fSAUseAllClusters=opt; return; }
  Bool_t GetSAUseAllClusters() const { return fSAUseAllClusters; }

  void SetMaxSPDcontrForSAToUseAllClusters(Int_t contr=50) { fMaxSPDcontrForSAToUseAllClusters=contr; return; }
  Int_t GetMaxSPDcontrForSAToUseAllClusters() const { return fMaxSPDcontrForSAToUseAllClusters; }

  void   SetSAUsedEdxInfo(Bool_t opt=kTRUE) { fSAUsedEdxInfo=opt; return; }
  Bool_t GetSAUsedEdxInfo() const { return fSAUsedEdxInfo; }

  void   SetFindV0s(Bool_t find=kTRUE) { fFindV0s=find; return; }
  Bool_t GetFindV0s() const { return fFindV0s; }

  void SetStoreLikeSignV0s(Bool_t like=kFALSE) { fStoreLikeSignV0s=like; return; } 
  Bool_t GetStoreLikeSignV0s() const { return fStoreLikeSignV0s; } 

  void   SetLayersParameters();

  void   SetLayerToSkip(Int_t i) { fLayersToSkip[i]=1; return; }
  Int_t  GetLayersToSkip(Int_t i) const { return fLayersToSkip[i]; }

  void   SetUseUnfoldingInClusterFinderSPD(Bool_t use=kTRUE) { fUseUnfoldingInClusterFinderSPD=use; return; }
  Bool_t GetUseUnfoldingInClusterFinderSPD() const { return fUseUnfoldingInClusterFinderSPD; }
  void   SetUseUnfoldingInClusterFinderSDD(Bool_t use=kTRUE) { fUseUnfoldingInClusterFinderSDD=use; return; }
  Bool_t GetUseUnfoldingInClusterFinderSDD() const { return fUseUnfoldingInClusterFinderSDD; }
  void   SetUseUnfoldingInClusterFinderSSD(Bool_t use=kTRUE) { fUseUnfoldingInClusterFinderSSD=use; return; }
  Bool_t GetUseUnfoldingInClusterFinderSSD() const { return fUseUnfoldingInClusterFinderSSD; }

  void   SetUseBadChannelsInClusterFinderSSD(Bool_t use=kFALSE) { fUseBadChannelsInClusterFinderSSD=use; return; }
  Bool_t GetUseBadChannelsInClusterFinderSSD() const  { return fUseBadChannelsInClusterFinderSSD;  }   

  void   SetUseSDDCorrectionMaps(Bool_t use=kTRUE) {fUseSDDCorrectionMaps=use;}
  Bool_t GetUseSDDCorrectionMaps() const {return fUseSDDCorrectionMaps;}
  void   SetUseSDDClusterSizeSelection(Bool_t use=kTRUE) {fUseSDDClusterSizeSelection=use;}
  Bool_t GetUseSDDClusterSizeSelection() const {return fUseSDDClusterSizeSelection;}
  void   SetMinClusterChargeSDD(Float_t qcut=0.){fMinClusterChargeSDD=qcut;}
  Float_t GetMinClusterChargeSDD() const {return fMinClusterChargeSDD;}

  void   SetUseChargeMatchingInClusterFinderSSD(Bool_t use=kTRUE) { fUseChargeMatchingInClusterFinderSSD=use; return; }
  Bool_t GetUseChargeMatchingInClusterFinderSSD() const { return fUseChargeMatchingInClusterFinderSSD; }

  void   SetUseCosmicRunShiftsSSD(Bool_t use=kFALSE) { fUseCosmicRunShiftsSSD=use; return; }
  Bool_t GetUseCosmicRunShiftsSSD() const { return fUseCosmicRunShiftsSSD; }

  // SPD Tracklets (D. Elia)
  void    SetTrackleterPhiWindow(Float_t w=0.08) {fTrackleterPhiWindow=w;}
  void    SetTrackleterThetaWindow(Float_t w=0.025) {fTrackleterThetaWindow=w;}
  void    SetTrackleterPhiShift(Float_t w=0.0045) {fTrackleterPhiShift=w;}
  Float_t GetTrackleterPhiWindow() const {return fTrackleterPhiWindow;}
  Float_t GetTrackleterThetaWindow() const {return fTrackleterThetaWindow;}
  Float_t GetTrackleterPhiShift() const {return fTrackleterPhiShift;}
  void    SetTrackleterRemoveClustersFromOverlaps(Bool_t use=kTRUE) { fTrackleterRemoveClustersFromOverlaps=use; return; }
  Bool_t  GetTrackleterRemoveClustersFromOverlaps() const { return fTrackleterRemoveClustersFromOverlaps; }
  void    SetTrackleterPhiOverlapCut(Float_t w=0.005) {fTrackleterPhiOverlapCut=w;}
  void    SetTrackleterZetaOverlapCut(Float_t w=0.05) {fTrackleterZetaOverlapCut=w;}
  Float_t GetTrackleterPhiOverlapCut() const {return fTrackleterPhiOverlapCut;}
  Float_t GetTrackleterZetaOverlapCut() const {return fTrackleterZetaOverlapCut;}
  void    SetTrackleterPhiRotationAngle(Float_t w=0.0) {fTrackleterPhiRotationAngle=w;}
  Float_t GetTrackleterPhiRotationAngle() const {return fTrackleterPhiRotationAngle;}
  //
  void    SetTrackleterNStdDevCut(Float_t f=1.)          {fTrackleterNStdDev = f<0.01 ? 0.01 : f;}
  Float_t GetTrackleterNStdDevCut()               const  {return fTrackleterNStdDev;}
  void    SetTrackleterScaleDThetaBySin2T(Bool_t v=kFALSE)  {fScaleDTBySin2T = v;}
  Bool_t  GetTrackleterScaleDThetaBySin2T()       const  {return fScaleDTBySin2T;}
  //
  void   SetSPDRemoveNoisyFlag(Bool_t value) {fSPDRemoveNoisyFlag = value;}
  Bool_t GetSPDRemoveNoisyFlag() const {return fSPDRemoveNoisyFlag;}
  void   SetSPDRemoveDeadFlag(Bool_t value) {fSPDRemoveDeadFlag = value;}
  Bool_t GetSPDRemoveDeadFlag() const {return fSPDRemoveDeadFlag;}
  
  //
  void    SetAlignFilterCosmics(Bool_t b=kTRUE) {fAlignFilterCosmics=b;}
  void    SetAlignFilterCosmicMergeTracks(Bool_t b=kTRUE) {fAlignFilterCosmicMergeTracks=b;} 
  void    SetAlignFilterMinITSPoints(Int_t n=4) {fAlignFilterMinITSPoints=n;}
  void    SetAlignFilterMinITSPointsMerged(Int_t n=4) {fAlignFilterMinITSPointsMerged=n;}
  void    SetAlignFilterOnlyITSSATracks(Bool_t b=kTRUE) {fAlignFilterOnlyITSSATracks=b;}
  void    SetAlignFilterOnlyITSTPCTracks(Bool_t b=kFALSE) {fAlignFilterOnlyITSTPCTracks=b;}
  void    SetAlignFilterUseLayer(Int_t ilay,Bool_t use) {fAlignFilterUseLayer[ilay]=use;}
  void    SetAlignFilterSkipExtra(Bool_t b=kFALSE) {fAlignFilterSkipExtra=b;}
  void    SetAlignFilterMaxMatchingAngle(Float_t max=0.085/*5deg*/) {fAlignFilterMaxMatchingAngle=max;}
  void    SetAlignFilterMinAngleWrtModulePlanes(Float_t min=0.52/*30deg*/) {fAlignFilterMinAngleWrtModulePlanes=min;}
  void    SetAlignFilterMinPt(Float_t min=0.) {fAlignFilterMinPt=min;}          
  void    SetAlignFilterMaxPt(Float_t max=1.e10) {fAlignFilterMaxPt=max;}          
  void    SetAlignFilterFillQANtuples(Bool_t b=kTRUE) {fAlignFilterFillQANtuples=b;}     
  Bool_t  GetAlignFilterCosmics() const {return fAlignFilterCosmics;}
  Bool_t  GetAlignFilterCosmicMergeTracks() const {return fAlignFilterCosmicMergeTracks;} 
  Int_t   GetAlignFilterMinITSPoints() const {return fAlignFilterMinITSPoints;}
  Int_t   GetAlignFilterMinITSPointsMerged() const {return fAlignFilterMinITSPointsMerged;}
  Bool_t  GetAlignFilterOnlyITSSATracks() const {return fAlignFilterOnlyITSSATracks;}
  Bool_t  GetAlignFilterOnlyITSTPCTracks() const {return fAlignFilterOnlyITSTPCTracks;}
  Bool_t  GetAlignFilterUseLayer(Int_t i) const {return fAlignFilterUseLayer[i];}
  Bool_t  GetAlignFilterSkipExtra() const {return fAlignFilterSkipExtra;}
  Float_t GetAlignFilterMaxMatchingAngle() const {return fAlignFilterMaxMatchingAngle;}
  Float_t GetAlignFilterMinAngleWrtModulePlanes() const {return fAlignFilterMinAngleWrtModulePlanes;}
  Float_t GetAlignFilterMinPt() const {return fAlignFilterMinPt;}          
  Float_t GetAlignFilterMaxPt() const {return fAlignFilterMaxPt;}          
  Bool_t  GetAlignFilterFillQANtuples() const {return fAlignFilterFillQANtuples;}     

  // Multiplicity Reconstructor
  Float_t GetMultCutPxDrSPDin()                 const {return fMultCutPxDrSPDin;}
  Float_t GetMultCutPxDrSPDout()                const {return fMultCutPxDrSPDout;}
  Float_t GetMultCutPxDz()                      const {return fMultCutPxDz;}
  Float_t GetMultCutDCArz()                     const {return fMultCutDCArz;}
  Float_t GetMultCutMinElectronProbTPC()        const {return fMultCutMinElectronProbTPC;}
  Float_t GetMultCutMinElectronProbESD()        const {return fMultCutMinElectronProbESD;}
  Float_t GetMultCutMinP()                      const {return fMultCutMinP;}
  Float_t GetMultCutMinRGamma()                 const {return fMultCutMinRGamma;}
  Float_t GetMultCutMinRK0()                    const {return fMultCutMinRK0;}
  Float_t GetMultCutMinPointAngle()             const {return fMultCutMinPointAngle;}
  Float_t GetMultCutMaxDCADauther()             const {return fMultCutMaxDCADauther;}
  Float_t GetMultCutMassGamma()                 const {return fMultCutMassGamma;}
  Float_t GetMultCutMassGammaNSigma()           const {return fMultCutMassGammaNSigma;}
  Float_t GetMultCutMassK0()                    const {return fMultCutMassK0;}
  Float_t GetMultCutMassK0NSigma()              const {return fMultCutMassK0NSigma;}
  Float_t GetMultCutChi2cGamma()                const {return fMultCutChi2cGamma;}
  Float_t GetMultCutChi2cK0()                   const {return fMultCutChi2cK0;}
  Float_t GetMultCutGammaSFromDecay()           const {return fMultCutGammaSFromDecay;}
  Float_t GetMultCutK0SFromDecay()              const {return fMultCutK0SFromDecay;}
  Float_t GetMultCutMaxDCA()                    const {return fMultCutMaxDCA;}
  //
  void    SetMultCutPxDrSPDin(Float_t v=0.1)             { fMultCutPxDrSPDin = v;}
  void    SetMultCutPxDrSPDout(Float_t v=0.15)           { fMultCutPxDrSPDout = v;}
  void    SetMultCutPxDz(Float_t v=0.2)                  { fMultCutPxDz = v;}
  void    SetMultCutDCArz(Float_t v=0.5)                 { fMultCutDCArz = v;}
  void    SetMultCutMinElectronProbTPC(Float_t v=0.5)    { fMultCutMinElectronProbTPC = v;}
  void    SetMultCutMinElectronProbESD(Float_t v=0.1)    { fMultCutMinElectronProbESD = v;}
  void    SetMultCutMinP(Float_t v=0.05)                 { fMultCutMinP = v;}
  void    SetMultCutMinRGamma(Float_t v=2.)              { fMultCutMinRGamma = v;}
  void    SetMultCutMinRK0(Float_t v=1.)                 { fMultCutMinRK0 = v;}
  void    SetMultCutMinPointAngle(Float_t v=0.98)        { fMultCutMinPointAngle = v;}
  void    SetMultCutMaxDCADauther(Float_t v=0.5)         { fMultCutMaxDCADauther = v;}
  void    SetMultCutMassGamma(Float_t v=0.03)            { fMultCutMassGamma = v;}
  void    SetMultCutMassGammaNSigma(Float_t v=5.)        { fMultCutMassGammaNSigma = v;}
  void    SetMultCutMassK0(Float_t v=0.03)               { fMultCutMassK0 = v;}
  void    SetMultCutMassK0NSigma(Float_t v=5.)           { fMultCutMassK0NSigma = v;}
  void    SetMultCutChi2cGamma(Float_t v=2.)             { fMultCutChi2cGamma = v;}
  void    SetMultCutChi2cK0(Float_t v=2.)                { fMultCutChi2cK0 = v;}
  void    SetMultCutGammaSFromDecay(Float_t v=-10.)      { fMultCutGammaSFromDecay = v;}
  void    SetMultCutK0SFromDecay(Float_t v=-10.)         { fMultCutK0SFromDecay = v;}
  void    SetMultCutMaxDCA(Float_t v=1.)                 { fMultCutMaxDCA = v;}
  //
  AliESDV0Params *GetESDV0Params() const {return fESDV0Params;}
  //
  // Lorentz angle
  Bool_t  GetCorrectLorentzAngleSPD() const {return fCorrectLorentzAngleSPD;}
  Float_t GetTanLorentzAngleHolesSPD() const {return fTanLorentzAngleHolesSPD;}
  Bool_t  GetCorrectLorentzAngleSSD() const {return fCorrectLorentzAngleSSD;}
  Float_t GetTanLorentzAngleHolesSSD() const {return fTanLorentzAngleHolesSSD;}
  Float_t GetTanLorentzAngleElectronsSSD() const {return fTanLorentzAngleElectronsSSD;}

  void SetCorrectLorentzAngleSPD(Bool_t flag) {fCorrectLorentzAngleSPD=flag;}
  void SetTanLorentzAngleHolesSPD(Float_t la) {fTanLorentzAngleHolesSPD=la;}
  void SetCorrectLorentzAngleSSD(Bool_t flag) {fCorrectLorentzAngleSSD=flag;}
  void SetTanLorentzAngleHolesSSD(Float_t la) {fTanLorentzAngleHolesSSD=la;}
  void SetTanLorentzAngleElectronsSSD(Float_t la) {fTanLorentzAngleElectronsSSD=la;}

  //
  enum {kMaxClusterPerLayer=70000}; //7000*10;   // max clusters per layer
  enum {kMaxClusterPerLayer5=28000};//7000*10*2/5;  // max clusters per layer
  enum {kMaxClusterPerLayer10=14000};//7000*10*2/10; // max clusters per layer
  enum {kMaxClusterPerLayer20=7000};//7000*10*2/20; // max clusters per layer

 protected:
  //
  static const Int_t fgkLayersNotToSkip[AliITSgeomTGeo::kNLayers]; // array with layers not to skip
  static const Int_t fgkLastLayerToTrackTo=0;     // innermost layer
  static const Int_t fgkMaxDetectorPerLayer=1000; // max clusters per layer
  static const Double_t fgkriw=80.0;              // TPC inner wall radius
  static const Double_t fgkdiw=0.0053;            // TPC inner wall x/X0
  static const Double_t fgkX0iw=30.0;             // TPC inner wall X0 
  static const Double_t fgkrcd=61.0;              // TPC central drum radius
  static const Double_t fgkdcd=0.0053;            // TPC central drum x/X0
  static const Double_t fgkX0cd=30.0;             // TPC central drum X0
  static const Double_t fgkyr=12.8;               // TPC rods y (tracking c.s.)
  static const Double_t fgkdr=0.03;               // TPC rods x/X0
  static const Double_t fgkzm=0.2;                // TPC membrane z
  static const Double_t fgkdm=0.40;               // TPC membrane x/X0
  static const Double_t fgkrs=50.0;               // ITS screen radius
  static const Double_t fgkds=0.001;              // ITS screed x/X0
  static const Double_t fgkrInsideITSscreen=49.0; // inside ITS screen radius
  static const Double_t fgkrInsideSPD1=3.5;       // inside SPD1 radius
  static const Double_t fgkrPipe=3.;              // pipe radius
  static const Double_t fgkrInsidePipe=2.7;       // inside pipe radius
  static const Double_t fgkrOutsidePipe=3.3;      // outside pipe radius
  static const Double_t fgkdPipe=0.0028;          // pipe x/X0
  static const Double_t fgkrInsideShield[2]; // inside SPD (0) SDD (1) shield radius
  static const Double_t fgkrOutsideShield[2]; // outside SPD (0) SDD (1) shield radius
  static const Double_t fgkdshield[2];        // SPD (0) SDD (1) shield x/X0
  static const Double_t fgkX0shield[2];       // SPD (0) SDD (1) shield X0
  static const Double_t fgkX0Air=21.82;       // air X0
  static const Double_t fgkX0Be=65.19;        // Berillium X0
  static const Double_t fgkBoundaryWidth=0.2; // to define track at detector boundary
  static const Double_t fgkDeltaXNeighbDets=0.5; // max difference in radius between neighbouring detectors 
  static const Double_t fgkSPDdetzlength=6.960;     // SPD ladder length in z (=7.072-2*0.056)
  static const Double_t fgkSPDdetxlength=1.298;     // SPD ladder length in x (=1.410-2*0.056)


  Int_t  fTracker;  // ITS tracker to be used (see AliITSReconstructor)
  Bool_t fITSonly;  // tracking only in ITS (no TPC)
  Int_t  fVertexer; // ITS vertexer to be used (see AliITSReconstructor)
  Int_t  fClusterFinder; // ITS cf to be used (see AliITSReconstructor)
  Int_t  fPID;      // ITS PID method to be used (see AliITSReconstructor)


  // SPD 3D Vertexer configuration
  Float_t fVtxr3DZCutWide;    // Z extension of the wide fiducial region for vertexer 3D
  Float_t fVtxr3DRCutWide;    // R extension of the wide fiducial region for vertexer 3D
  Float_t fVtxr3DZCutNarrow;  // Z extension of the narrow fiducial region for vertexer 3D
  Float_t fVtxr3DRCutNarrow;  // R extension of the narrow fiducial region for vertexer 3D
  Float_t fVtxr3DPhiCutLoose; // loose deltaPhi cut to define tracklets in vertexer 3D
  Float_t fVtxr3DPhiCutTight; // tight deltaPhi cut to define tracklets in vertexer 3D
  Float_t fVtxr3DDCACut;      // cut on tracklet-to-tracklet DCA in vertexer3D
  Int_t   fVtxr3DPileupAlgo;  // pileup algorithm (0 = VtxZ, 1 = 3D - 2 step, 2 = 3D all in once)
  UChar_t fVtxr3DHighMultAlgo; // downscaling if 0 - traces if 1

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
  // search for extra clusters
  Bool_t   fSearchForExtras; // swicth yes/no for the search of extra-clusters in RefitInward step
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

  Bool_t fSkipSubdetsNotInTriggerCluster; // skip the subdetectors that are not in the trigger cluster

  Int_t fUseTGeoInTracker; // use TGeo to get material budget in tracker MI
  Double_t fStepSizeTGeo; // step size (cm)
                     // in AliITStrackerMI::CorrectFor*Material methods
  Bool_t fAllowSharedClusters; // if kFALSE don't set to kITSin tracks with shared clusters (MI)
  Int_t fClusterErrorsParam; // parametrization for cluster errors (MI), see AliITSRecoParam::GetError()
  Float_t fClusterMisalErrorY[AliITSgeomTGeo::kNLayers]; // [cm] additional error on cluster Y pos. due to misalignment (MI,SA)
  Float_t fClusterMisalErrorZ[AliITSgeomTGeo::kNLayers]; // [cm] additional error on cluster Z pos. due to misalignment (MI,SA)
  Float_t fClusterMisalErrorYBOn[AliITSgeomTGeo::kNLayers]; // [cm] additional error on cluster Y pos. due to misalignment (MI,SA)
  Float_t fClusterMisalErrorZBOn[AliITSgeomTGeo::kNLayers]; // [cm] additional error on cluster Z pos. due to misalignment (MI,SA)

  Bool_t fUseAmplitudeInfo[AliITSgeomTGeo::kNLayers]; // use cluster charge in cluster-track matching (SDD,SSD) (MI)

  // Plane Efficiency evaluation
  Bool_t fComputePlaneEff;  // flag to enable computation of PlaneEfficiency
  Bool_t fHistoPlaneEff;  // flag to enable auxiliary PlaneEff histograms (e.g. residual distributions)
  Bool_t fUseTrackletsPlaneEff; // flag to enable estimate of SPD PlaneEfficiency using tracklets
  Bool_t fMCTrackletsPlaneEff; // flag to enable the use of MC info for corrections (SPD PlaneEff using tracklets)
  Bool_t fBkgTrackletsPlaneEff; // flag to evaluate background instead of normal use (SPD PlaneEff using tracklets)
  Float_t fTrackleterPhiWindowL1; // Search window in phi for inner layer (1) (SPD PlaneEff using tracklets)
  Float_t fTrackleterPhiWindowL2; // Search window in phi for outer layer (2) (SPD PlaneEff using tracklets)
  Float_t fTrackleterZetaWindowL1; // Search window in zeta for inner layer (1) (SPD PlaneEff using tracklets)
  Float_t fTrackleterZetaWindowL2; // Search window in zeta for outer layer (2) (SPD PlaneEff using tracklets)
  Bool_t fUpdateOncePerEventPlaneEff; // option to update chip efficiency once/event (to avoid doubles)
  Int_t  fMinContVtxPlaneEff; // min number of contributors to ESD vtx for SPD PlaneEff using tracklets
  Int_t  fIPlanePlaneEff; // index of the plane (in the range [-1,5]) to study the efficiency (-1 ->Tracklets)
  Bool_t fReadPlaneEffFromOCDB; // enable initial reading of Plane Eff statistics from OCDB
                               // The analized events would be used to increase the statistics
  Double_t fMinPtPlaneEff;  // minimum p_t of the track to be used for Plane Efficiency evaluation
  Int_t  fMaxMissingClustersPlaneEff;  // max n. of (other) layers without a cluster associated to the track
  Int_t  fMaxMissingClustersOutPlaneEff;  // max n. of outermost layers without a cluster associated to the track
  Bool_t fRequireClusterInOuterLayerPlaneEff; // if kTRUE, then only tracks with an associated cluster on the closest
  Bool_t fRequireClusterInInnerLayerPlaneEff; // outer/inner layer are used. It has no effect for outermost/innermost layer
  Bool_t fOnlyConstraintPlaneEff;  // if kTRUE, use only constrained tracks at primary vertex for Plane Eff.
  Double_t fNSigXFromBoundaryPlaneEff;  // accept one track for PlaneEff if distance from border (in loc x or z)
  Double_t fNSigZFromBoundaryPlaneEff;  // is greater than fNSigXFromBoundaryPlaneEff * Track_precision

  Bool_t fImproveWithVertex;    // use the method AliITStrackV2::Improve() to point to the vertex during prolongation
  Bool_t fExtendedEtaAcceptance;  // enable jumping from TPC to SPD at large eta (MI)
  Bool_t fUseBadZonesFromOCDB; // enable using OCDB info on dead modules and chips (MI)
  Bool_t fUseSingleBadChannelsFromOCDB; // enable using OCDB info on bad single SPD pixels and SDD anodes (MI)
  Float_t fMinFractionOfBadInRoad; // to decide whether to skip the layer (MI)
  Bool_t fAllowProlongationWithEmptyRoad; // allow to prolong even if road is empty (MI)
  Int_t fInwardFlagSA;           // flag for inward track finding in SA
  Int_t fOuterStartLayerSA;      // outer ITS layer to start track in SA outward
  Int_t fInnerStartLayerSA;      // inner ITS layer to start track in SA inward
  Int_t fMinNPointsSA;           // min. number of ITS clusters for a SA track
  Double_t fFactorSAWindowSizes; // larger window sizes in SA
  Int_t fNLoopsSA;               // number of loops in tracker SA
  Double_t fMinPhiSA;               // minimum phi value for SA windows
  Double_t fMaxPhiSA;               // maximum phi value for SA windows
  Double_t fMinLambdaSA;            // minimum lambda value for SA windows
  Double_t fMaxLambdaSA;            // maximum lambda value for SA windows
  Float_t  fMinClusterChargeSA;     // minimum SDD,SSD cluster charge for SA tarcker
  Bool_t fSAOnePointTracks; // one-cluster tracks in SA (only for cosmics!)
  Bool_t fSAUseAllClusters; // do not skip clusters used by MI (same track twice in AliESDEvent!)
  Int_t fMaxSPDcontrForSAToUseAllClusters; // maximum nContr of SPD vertex for which trackerSA will reuse all ITS clusters
  Bool_t fSAUsedEdxInfo;   // use/not use dE/dx in ITS for assign mass hypothesis

  Bool_t fSelectBestMIP03;          // (MI) Multiply norm chi2 by interpolated one in hypthesis analysis
  Bool_t fFlagFakes;                // (MI) preform shared cluster analysis and flag candidates for fakes
  Bool_t fUseImproveKalman;         // (MI) Use ImproveKalman version of AliITSTrackV2 instead of Improve

  Bool_t fFindV0s;  // flag to enable V0 finder (MI)
  Bool_t fStoreLikeSignV0s; // flag to store like-sign V0s (MI)

  // cluster unfolding in ITS cluster finders
  Bool_t fUseUnfoldingInClusterFinderSPD; // SPD
  Bool_t fUseUnfoldingInClusterFinderSDD; // SDD
  Bool_t fUseUnfoldingInClusterFinderSSD; // SSD

  Bool_t fUseBadChannelsInClusterFinderSSD; // flag to switch on bad channels in CF SSD

  Bool_t  fUseSDDCorrectionMaps; // flag for use of SDD maps in C.F.
  Bool_t  fUseSDDClusterSizeSelection; // cut on SDD cluster size
  Float_t fMinClusterChargeSDD; // cut on SDD cluster charge

  Bool_t fUseChargeMatchingInClusterFinderSSD; // SSD

  // SPD Tracklets (D. Elia)
  Float_t fTrackleterPhiWindow;                    // Search window in phi
  Float_t fTrackleterThetaWindow;                   // Search window in theta
  Float_t fTrackleterPhiShift;                     // Phi shift reference value (at 0.5 T) 
  Bool_t  fTrackleterRemoveClustersFromOverlaps;   // Option to skip clusters in the overlaps
  Float_t fTrackleterPhiOverlapCut;                // Fiducial window in phi for overlap cut
  Float_t fTrackleterZetaOverlapCut;               // Fiducial window in eta for overlap cut
  Float_t fTrackleterPhiRotationAngle;             // Angle to rotate cluster in the SPD inner layer for combinatorial reco only
  Float_t fTrackleterNStdDev;      // cut on the number of standard deviations
  Bool_t  fScaleDTBySin2T;         // scale Dtheta by 1/sin^2(theta)

  Bool_t fUseCosmicRunShiftsSSD; // SSD time shifts for cosmic run 2007/2008 (use for data taken up to 18 sept 2008)


   // SPD flags to specify whether noisy and dead pixels 
  // should be removed at the local reconstruction step (default and safe way is true for both)
  Bool_t  fSPDRemoveNoisyFlag;  // Flag saying whether noisy pixels should be removed
  Bool_t  fSPDRemoveDeadFlag;   // Flag saying whether dead pixels should be removed
  
  // VertexerFast configuration
  Float_t fVertexerFastSmearX;  // gaussian sigma for x MC vertex smearing 
  Float_t fVertexerFastSmearY;  // gaussian sigma for y MC vertex smearing
  Float_t fVertexerFastSmearZ;  // gaussian sigma for z MC vertex smearing

  // PWG1/AliAlignmentDataFilterITS configuration
  Bool_t  fAlignFilterCosmics;            // flag for cosmics case
  Bool_t  fAlignFilterCosmicMergeTracks;  // merge cosmic tracks
  Int_t   fAlignFilterMinITSPoints;       // min points per track
  Int_t   fAlignFilterMinITSPointsMerged; // min points for merged tracks
  Bool_t  fAlignFilterOnlyITSSATracks;    // only ITS SA tracks
  Bool_t  fAlignFilterOnlyITSTPCTracks;   // only ITS+TPC tracks
  Bool_t  fAlignFilterUseLayer[AliITSgeomTGeo::kNLayers]; // layers to use 
  Bool_t  fAlignFilterSkipExtra;          // no extra cls in array
  Float_t fAlignFilterMaxMatchingAngle;   // matching for cosmics
  Float_t fAlignFilterMinAngleWrtModulePlanes; // min angle track-to-sensor
  Float_t fAlignFilterMinPt;              // min pt
  Float_t fAlignFilterMaxPt;              // max pt
  Bool_t  fAlignFilterFillQANtuples;      // fill QA ntuples  

  // Multiplicity reconstructor settings
  // cuts for flagging secondaries
  Float_t fMultCutPxDrSPDin;              // max P*DR for primaries involving at least 1 SPD
  Float_t fMultCutPxDrSPDout;             // max P*DR for primaries not involving any SPD
  Float_t fMultCutPxDz;                   // max P*DZ for primaries
  Float_t fMultCutDCArz;                  // max DR or DZ for primares
  //
  // cuts for flagging tracks in V0s
  Float_t fMultCutMinElectronProbTPC;     // min probability for e+/e- PID involving TPC
  Float_t fMultCutMinElectronProbESD;     // min probability for e+/e- PID not involving TPC
  //
  Float_t fMultCutMinP;                   // min P of V0
  Float_t fMultCutMinRGamma;              // min transv. distance from ESDVertex to V0 for gammas
  Float_t fMultCutMinRK0;                 // min transv. distance from ESDVertex to V0 for K0s
  Float_t fMultCutMinPointAngle;          // min pointing angle cosine
  Float_t fMultCutMaxDCADauther;          // max DCA of daughters at V0
  Float_t fMultCutMassGamma;              // max gamma mass
  Float_t fMultCutMassGammaNSigma;        // max standard deviations from 0 for gamma
  Float_t fMultCutMassK0;                 // max K0 mass difference from PGD value
  Float_t fMultCutMassK0NSigma;           // max standard deviations for K0 mass from PDG value
  Float_t fMultCutChi2cGamma;             // max constrained chi2 cut for gammas
  Float_t fMultCutChi2cK0;                // max constrained chi2 cut for K0s
  Float_t fMultCutGammaSFromDecay;        // min path*P for gammas
  Float_t fMultCutK0SFromDecay;           // min path*P for K0s
  Float_t fMultCutMaxDCA;                 // max DCA for V0 at ESD vertex
  //
  // Lorentz angle
  Bool_t fCorrectLorentzAngleSPD;         // flag to enable correction
  Float_t fTanLorentzAngleHolesSPD;       // angle for holes in SPD
  Bool_t fCorrectLorentzAngleSSD;         // flag to enable correction
  Float_t fTanLorentzAngleHolesSSD;       // tan(angle) for holes in SSD @ B = 0.5 T
  Float_t fTanLorentzAngleElectronsSSD;   // tan(angle) for electrons in SSD @ B = 0.5 T

 private:
  AliESDV0Params * fESDV0Params;  // declare the AliESDV0Params to be able to used in AliITSV0Finder

  AliITSRecoParam(const AliITSRecoParam & param);
  AliITSRecoParam & operator=(const AliITSRecoParam &param);

  ClassDef(AliITSRecoParam,39) // ITS reco parameters
};

#endif
