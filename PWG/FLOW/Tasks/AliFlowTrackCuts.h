/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliFlowTrackESDCuts:
// A cut class for ESD, AOD and MC particles for the flow framework
// author: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)
// mods:   Redmer A. Bertens (rbertens@cern.ch)

#ifndef ALIFLOWTRACKCUTS_H
#define ALIFLOWTRACKCUTS_H

#include <TMatrix.h>
#include <TList.h>
#include "AliFlowTrackSimpleCuts.h"
#include "AliFlowTrackSimple.h"
#include "AliESDtrackCuts.h"
#include "TMCProcess.h"
#include "AliESDtrack.h"
#include "AliMuonTrackCuts.h"  // XZhang 20120604
#include "AliPID.h"
#include "AliESDpid.h"
#include "TF2.h"


class TBrowser;
class TArrayD;
class AliVParticle;
class AliMCParticle;
class AliFlowTrack;
class AliMCEvent;
class AliInputEventHandler;
class AliVEvent;
class AliMultiplicity; 
class AliAODTracklets;  // XZhang 20120615
class AliAODTrack;
class AliESDtrack;
class AliESDPmdTrack;
class AliFlowBayesianPID;
class AliESDkink;
class AliESDv0;
class AliESDVZERO;
class AliPIDResponse;

class AliFlowTrackCuts : public AliFlowTrackSimpleCuts {

 public:
  AliFlowTrackCuts();
  AliFlowTrackCuts(const char* name);
  AliFlowTrackCuts(const AliFlowTrackCuts& someCuts);
  AliFlowTrackCuts& operator=(const AliFlowTrackCuts& someCuts);
  virtual ~AliFlowTrackCuts();

  static AliFlowTrackCuts* GetAODTrackCutsForFilterBit(UInt_t bit = 1, TString suffix = "");
  static AliFlowTrackCuts* GetStandardTPCStandaloneTrackCuts();
  static AliFlowTrackCuts* GetStandardTPCStandaloneTrackCuts2010();
  static AliFlowTrackCuts* GetStandardGlobalTrackCuts2010();
  static AliFlowTrackCuts* GetStandardITSTPCTrackCuts2009(Bool_t selPrimaries=kTRUE);
  static AliFlowTrackCuts* GetStandardVZEROOnlyTrackCuts();
  static AliFlowTrackCuts* GetStandardVZEROOnlyTrackCuts2010();
  static AliFlowTrackCuts* GetStandardVZEROOnlyTrackCuts2011();
  // beta test, don't use for results
  static AliFlowTrackCuts* GetBetaVZEROOnlyTrackCuts();
  // end of beta test
  static AliFlowTrackCuts* GetStandardMuonTrackCuts(Bool_t isMC=kFALSE, Int_t passN=2);  // XZhang 20120604

  Int_t Count(AliVEvent* event=NULL);

  enum trackParameterType { kMC, 
                            kGlobal, 
                            kTPCstandalone, 
                            kSPDtracklet,
                            kPMD,
                            kV0,    //neutral reconstructed v0 particle
                            kVZERO, //forward VZERO detector
                            kMUON,  // XZhang 20120604
                            kKink,
                            kAODFilterBit,
                            kUserA, // reserved for custom cuts
                            kUserB, // reserved for custom cuts
                            kBetaVZERO, // temporary enum for beta testing of new vzero calibration
                            kDeltaVZERO, // temporary enum for beta testing of new vzero calibration
                            kKappaVZERO, // permanent enum for beta testing of new vzero calibration, 
                            kHotfixHI   // fix for HI runs
                          };
  enum trackParameterMix  { kPure, 
                            kTrackWithMCkine, 
                            kTrackWithMCPID, 
                            kTrackWithMCpt, 
                            kTrackWithPtFromFirstMother,
                            kTrackWithTPCInnerParams,
                            kTrackWithTPCstandalone
                          };
  enum PIDsource {
                   kTPCpid,      // default TPC pid (via GetTPCpid)
                   kTOFpid,      // default TOF pid (via GetTOFpid)
                   kTOFbayesian, // TOF bayesian pid (F.Noferini)
                   kTOFbeta,     // asymmetric cuts of TOF beta signal
                   kTPCdedx,      // asymmetric cuts of TPC dedx signal
                   kTOFbetaSimple, //simple TOF only cut
                   kTPCbayesian, //bayesian cutTPC
		   kTPCNuclei,   // added by Natasha for Nuclei
                   kTPCTOFNsigma, // simple cut on combined tpc tof nsigma
                   kTPCTOFNsigmaPurity, // purity>0.8 cut on combined tpc tof nsigma
				   kTPCTPCTOFNsigma ////cut on sigma tpc below certain pt, on combined tpc tof sigma above (AOD)
                   };

  //setters (interface to AliESDtrackCuts)
  Int_t MaxSharedITSClusterCuts(AliESDtrack* track);
  Double_t MaxChi2perITSClusterCuts(AliESDtrack* track);

  void SetMaxSharedITSCluster(Int_t b){fCutITSclusterShared = kTRUE; fMaxITSclusterShared = b;}
  void SetMaxChi2perITSCluster(Double_t b){fCutITSChi2 = kTRUE; fMaxITSChi2 = b;}
  void SetCutTPCSecbound( Bool_t a, Double_t ptmin=0.2 ) {fCutTPCSecbound = a; fCutTPCSecboundMinpt=ptmin;}
  void SetCutTPCSecboundVar( Bool_t a ) {fCutTPCSecboundVar = a;}
  void SetMinNClustersTPC( Int_t a ) {fCutNClustersTPC=kTRUE; fNClustersTPCMin=a;}
  void SetMinNClustersITS( Int_t a ) {fCutNClustersITS=kTRUE; fNClustersITSMin=a;}
  void SetClusterRequirementITS( AliESDtrackCuts::Detector det,
                                 AliESDtrackCuts::ITSClusterRequirement req = AliESDtrackCuts::kOff )
                                 { InitESDcuts(); fAliESDtrackCuts->SetClusterRequirementITS(det,req); }
  void SetCutChi2PerClusterITS( Float_t a ) {fCutChi2PerClusterITS=kTRUE; fMaxChi2PerClusterITS=a;}
  void SetCutITSClusterGlobal( Bool_t a ) {fCutITSClusterGlobal=a;}
  void SetMaxChi2PerClusterTPC( Float_t a ) {fMaxChi2PerClusterTPC=a;fCutChi2PerClusterTPC=kTRUE;}
  void SetMinChi2PerClusterTPC( Float_t a ) {fMinChi2PerClusterTPC=a;fCutChi2PerClusterTPC=kTRUE;}
  void SetMaxFracSharedTPCCluster( Float_t a ) {fMaxFracSharedTPCCluster=a;fCutFracSharedTPCCluster=kTRUE;}
  void SetCutCrossedTPCRows( Int_t a, Float_t b) {fCutCrossedTPCRows=kTRUE; fMinNCrossedRows=a; fMinCrossedRowsOverFindableClusters=b;}
  void SetCutGoldenChi2( Double_t m ) {fCutGoldenChi2=kTRUE; fMaxGoldenChi2=m;}
  void SetRequireTOFSignal( Bool_t a ) {fRequireTOFSignal=a;}
  void SetMaxChi2PerClusterITS( Float_t a ) {InitESDcuts(); fAliESDtrackCuts->SetMaxChi2PerClusterITS(a);}
  void SetRequireTPCRefit( Bool_t a ) {InitESDcuts(); fAliESDtrackCuts->SetRequireTPCRefit(a);}
  void SetRequireTPCStandAlone( Bool_t a) {InitESDcuts(); fAliESDtrackCuts->SetRequireTPCStandAlone(a);}
  void SetRequireITSRefit( Bool_t a ) {InitESDcuts(); fAliESDtrackCuts->SetRequireITSRefit(a);}
  void SetRequireITSStandAlone( Bool_t a) {InitESDcuts(); fAliESDtrackCuts->SetRequireITSStandAlone(a);}
  void SetAcceptKinkDaughters( Bool_t a ) {InitESDcuts(); fAliESDtrackCuts->SetAcceptKinkDaughters(a);}
  void SetMaxDCAToVertexZ( Float_t a ) {InitESDcuts(); fAliESDtrackCuts->SetMaxDCAToVertexZ(a);fCutDCAToVertexZ=kTRUE;}
  void SetMaxDCAToVertexXY( Float_t a ) {InitESDcuts(); fAliESDtrackCuts->SetMaxDCAToVertexXY(a);fCutDCAToVertexXY=kTRUE;}
  void SetMaxDCAToVertexXYPtDep( const char* a ) {InitESDcuts(); fAliESDtrackCuts->SetMaxDCAToVertexXYPtDep(a);}
  void SetMaxDCAToVertexXYPtDepAOD( Bool_t a ) {fCutDCAToVertexXYPtDepAOD=a;}
  void SetMaxDCAToVertexXYAOD( Float_t a ) {fCutDCAToVertexXYAOD=kTRUE; fMaxDCAxyAOD=a;}
  void SetMaxDCAToVertexZAOD( Float_t a ) {fCutDCAToVertexZAOD=kTRUE; fMaxDCAzAOD=a;}
  void SetRequireSigmaToVertex(Bool_t a) {InitESDcuts(); fAliESDtrackCuts->SetRequireSigmaToVertex(a);}
  void SetMaxNsigmaToVertex(Float_t sigma=1e10) {InitESDcuts(); fAliESDtrackCuts->SetMaxNsigmaToVertex(sigma); }
  void SetDCAToVertex2D( Bool_t a ) {InitESDcuts(); fAliESDtrackCuts->SetDCAToVertex2D(a);}
  void SetEtaRange( Float_t r1, Float_t r2 ) { SetEtaMin(r1); SetEtaMax(r2); }
  void SetPtRange( Float_t r1, Float_t r2 ) { SetPtMin(r1); SetPtMax(r2); }
  void SetRequireCharge( Bool_t r ) {fRequireCharge=r;}
  void SetFakesAreOK( Bool_t b ) {fFakesAreOK=b;}
  void SetSPDtrackletDeltaPhiMax( Double_t m ) {fSPDtrackletDeltaPhiMax=m; fCutSPDtrackletDeltaPhi=kTRUE;}
  void SetSPDtrackletDeltaPhiMin( Double_t m ) {fSPDtrackletDeltaPhiMin=m; fCutSPDtrackletDeltaPhi=kTRUE;}
  void SetIgnoreTPCzRange( Double_t min, Double_t max ) 
                         { fIgnoreTPCzRange=kTRUE; fIgnoreTPCzRangeMin=min; fIgnoreTPCzRangeMax=max; }
  void SetAODfilterBit( UInt_t a ) {fAODFilterBit = a; fUseAODFilterBit = kTRUE;}  						 
  void SetMinimalTPCdedx(Double_t d=10.) {fMinimalTPCdedx=d; fCutMinimalTPCdedx=kTRUE;}
  void SetPmdDetPlane(Int_t pmdDet){fCutPmdDet=kTRUE; fPmdDet = pmdDet; }
  void SetPmdAdc(Float_t pmdAdc){fCutPmdAdc=kTRUE; fPmdAdc = pmdAdc; }
  void SetPmdNcell(Float_t pmdNcell) {fCutPmdNcell=kTRUE; fPmdNcell = pmdNcell; }						 
  void SetPriors(Float_t centr = 0); // set my favourite priors for Bayesian PID (requested if Bayesian PID is used)

  AliMuonTrackCuts *GetMuonTrackCuts() { InitMuonCuts(); return fMuonTrackCuts; }                           // XZhang 20121014
  void SetStandardMuonTrackCuts()      { InitMuonCuts(); fMuonTrackCuts->SetDefaultFilterMask(); return; }  // XZhang 20120604
  void SetIsMuonMC(Bool_t isMC)        { InitMuonCuts(); fMuonTrackCuts->SetIsMC(isMC);          return; }  // XZhang 20120604
  void SetMuonPassNumber(Int_t passN)  { InitMuonCuts(); fMuonTrackCuts->SetPassNumber(passN);   return; }  // XZhang 20121013
  void SetRunsMuon(const AliInputEventHandler* eventHandler) { if (fMuonTrackCuts) fMuonTrackCuts->SetRun(eventHandler); }  // XZhang 20120604

  void SetForceTPCstandalone(Bool_t b) {fForceTPCstandalone=b;}

  //Kinks
  void SetMinKinkAngle(Double_t a) {fMinKinkAngle=a;}
  void SetMinKinkRadius(Double_t r) {fMinKinkRadius=r;}
  void SetMaxKinkRAdius(Double_t r) {fMaxKinkRadius=r;}
  void SetMinKinkQt(Double_t m) {fMinKinkQt=m;}
  void SetMaxKinkQt(Double_t m) {fMaxKinkQt=m;}
  void SetMaxKinkInvMassKmu(Double_t m) {fMaxKinkInvMassKmu=m;}
  void SetMinKinkInvMassKmu(Double_t m) {fMinKinkInvMassKmu=m;}

  Double_t GetMinKinkAngle() const {return fMinKinkAngle;}
  Double_t GetMinKinkRadius() const {return fMinKinkRadius;}
  Double_t GetMaxKinkRadius() const {return fMaxKinkRadius;}
  Double_t GetMinKinkQt() const {return fMinKinkQt;}
  Double_t GetMaxKinkQt() const {return fMaxKinkQt;}
  Double_t GetMaxKinkInvMassKmu() const {return fMaxKinkInvMassKmu;}
  Double_t GetMinKinkInvMassKmu() const {return fMinKinkInvMassKmu;}

  Int_t GetMinNClustersTPC() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetMinNClusterTPC();}
  Int_t GetMinNClustersITS() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetMinNClustersITS();}
  AliESDtrackCuts::ITSClusterRequirement GetClusterRequirementITS( AliESDtrackCuts::Detector det ) const
                                 {if (!fAliESDtrackCuts) return AliESDtrackCuts::kOff;  return fAliESDtrackCuts->GetClusterRequirementITS(det); } 
  Float_t GetMaxChi2PerClusterTPC() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetMaxChi2PerClusterTPC();}
  Float_t GetMaxChi2PerClusterITS() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetMaxChi2PerClusterITS();}
  Bool_t GetRequireTPCRefit() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetRequireTPCRefit();}
  Bool_t GetRequireTPCStandAlone() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetRequireTPCStandAlone();}
  Bool_t GetRequireITSRefit() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetRequireITSRefit();}
  Bool_t GetRequireITSStandAlone() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetRequireITSStandAlone();}
  Bool_t GetAcceptKinkDaughters() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetAcceptKinkDaughters();}
  Float_t GetMaxDCAToVertexZ() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetMaxDCAToVertexZ();}
  Float_t GetMaxDCAToVertexXY() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetMaxDCAToVertexXY();}
  const char* GetMaxDCAToVertexXYPtDep() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetMaxDCAToVertexXYPtDep();}
  Bool_t GetRequireSigmaToVertex() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetRequireSigmaToVertex();}
  Float_t GetMaxNsigmaToVertex() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetMaxNsigmaToVertex(); }
  Bool_t GetDCAToVertex2D() const {if (!fAliESDtrackCuts) return 0; return fAliESDtrackCuts->GetDCAToVertex2D();}
  void GetEtaRange( Float_t& r1, Float_t& r2 ) const { r1=GetEtaMin(); r2=GetEtaMax(); }
  void GetPtRange( Float_t& r1, Float_t& r2 ) const { r1=GetPtMin(); r2=GetPtMax(); }
  Bool_t GetRequireCharge() const {return fRequireCharge;}
  Bool_t GetFakesAreOK() const {return fFakesAreOK;}
  Double_t GetSPDtrackletDeltaPhiMax() const {return fSPDtrackletDeltaPhiMax;}
  Double_t GetSPDtrackletDeltaPhiMin() const {return fSPDtrackletDeltaPhiMin;}
  UInt_t GetAODFilterBit() const {if (!fUseAODFilterBit) return 0; return fAODFilterBit;}
  Double_t GetMinimalTPCdedx() const {return fMinimalTPCdedx;}
  Int_t GetPmdDetPlane()const {return fPmdDet; }
  Float_t GetPmdAdc()const {return fPmdAdc;}
  Float_t GetPmdNcell() const {return fPmdNcell; }
  Float_t GetBeta(const AliVTrack* t, Bool_t QAmode = kFALSE);
  Float_t Getdedx(const AliESDtrack* t) const;
  Float_t GetBayesianProb() const {return fProbBayes;};
  AliFlowBayesianPID* GetBayesianResponse() const {return  fBayesianResponse;}

  Bool_t GetForceTPCstandalone() const {return fForceTPCstandalone;}

  void SetQA(Bool_t b=kTRUE) {if (b) DefineHistograms();}
  TList* GetQA() const {return fQA;}
  TH1* QAbefore(Int_t i) {return static_cast<TH1*>(static_cast<TList*>(fQA->At(0))->At(i));}
  TH1* QAafter(Int_t i) {return static_cast<TH1*>(static_cast<TList*>(fQA->At(1))->At(i));}  

  //MC stuff
  void SetIgnoreSignInMCPID( Bool_t b=kTRUE ) {fIgnoreSignInMCPID=b;}
  void SetCutMC( Bool_t b=kTRUE );
  void SetCutMChasTrackReferences(Bool_t b=kTRUE) {fCutMChasTrackReferences=b;}
  void SetMCprocessType( TMCProcess t ) { fMCprocessType = t; fCutMCprocessType=kTRUE; SetCutMC();}
  void SetMCisPrimary( Bool_t b=kTRUE ) { fMCisPrimary=b; fCutMCisPrimary=kTRUE; SetCutMC();}
  void SetMCPID( Int_t pid ) { fMCPID=pid; fCutMCPID=kTRUE; SetCutMC(); }
  void SetMCfirstMotherPID( Int_t pid ) { fMCfirstMotherPID=pid; fCutMCfirstMotherPID=kTRUE; SetCutMC(); }
  TMCProcess GetMCprocessType() const { return fMCprocessType; }
  Bool_t GetMCisPrimary() const {return fMCisPrimary;}
  Int_t GetMCPID() const {return fMCPID;}
  void SetRequireTransportBitForPrimaries(Bool_t b) {fRequireTransportBitForPrimaries=b; SetCutMC();}

  void SetParamType(trackParameterType paramType) {fParamType=paramType;}
  trackParameterType GetParamType() const {return fParamType;}
  static const char* GetParamTypeName(trackParameterType type);
  void SetParamMix(trackParameterMix paramMix) {fParamMix=paramMix;}
  trackParameterMix GetParamMix() const {return fParamMix;}

  virtual Bool_t IsSelected(TObject* obj, Int_t id=-666);
  virtual Bool_t IsSelectedMCtruth(TObject* obj, Int_t id=-666);
  AliVParticle* GetTrack() const {return fTrack;}
  AliMCParticle* GetMCparticle() const {return fMCparticle;}
  //AliFlowTrack* MakeFlowTrack() const;
  Bool_t FillFlowTrack(AliFlowTrack* track) const;
  //FillFlowTrackV0(TObjArray* trackCollection, Int_t trackIndex) const
  AliFlowTrack* FillFlowTrack(TObjArray* trackCollection, Int_t trackIndex) const;
  Bool_t IsPhysicalPrimary() const; 
  static Bool_t IsPhysicalPrimary(AliMCEvent* p, Int_t label, Bool_t requiretransported=kTRUE); 
  
  void SetMCevent(AliMCEvent* mcEvent) {fMCevent=mcEvent;}
  AliMCEvent* GetMCevent() const {return fMCevent;}
  void SetEvent(AliVEvent* event, AliMCEvent* mcEvent=NULL);
  AliVEvent* GetEvent() const {return fEvent;}
  Int_t GetNumberOfInputObjects() const;
  TObject* GetInputObject(Int_t i);
  void Clear(Option_t* option="");
  void ClearTrack(Option_t* option="");

  Double_t GetPmdEta(Float_t xPos, Float_t yPos, Float_t zPos);
  Double_t GetPmdPhi(Float_t xPos, Float_t yPos);  

  //PID
  void SetPID(AliPID::EParticleType pid, PIDsource s=kTOFpid, Double_t prob=0.9)
             {fParticleID=pid; fPIDsource=s; fParticleProbability=prob; fCutPID=kTRUE; InitPIDcuts();}
  AliPID::EParticleType GetParticleID() const {return fParticleID;}
  Bool_t GetCutPID() const {return fCutPID;}
  void SetTPCpidCuts(const TMatrixF* mat) {fTPCpidCuts=new TMatrixF(*mat);}
  void SetTOFpidCuts(const TMatrixF* mat) {fTOFpidCuts=new TMatrixF(*mat);}
  static const char* PIDsourceName(PIDsource s);
  AliESDpid& GetESDpid() {return fESDpid;}
  void SetAllowTOFmismatchFlag(Bool_t b=kTRUE) {fAllowTOFmismatchFlag=b;}
  Bool_t GetAllowTOFmismatchFlag() const {return fAllowTOFmismatchFlag;}
  void SetRequireStrictTOFTPCagreement(Bool_t b=kTRUE) {fRequireStrictTOFTPCagreement=b;}
  Bool_t GetRequireStrictTOFTPCagreement() const {return fRequireStrictTOFTPCagreement;}
  void SetRejectElectronsWithTPCpid(Bool_t b=kTRUE) {fCutRejectElectronsWithTPCpid=b;}
  void SetLinearizeVZEROresponse( Bool_t b=kTRUE ) {fLinearizeVZEROresponse=b;}

  //these should maybe be protected
  Bool_t PassesCuts(AliVParticle* track);
  Bool_t PassesESDcuts(AliESDtrack* track);
  Bool_t PassesAODcuts(const AliAODTrack* track, Bool_t passFid=kTRUE);
  Bool_t PassesPMDcuts(const AliESDPmdTrack* track);
  Bool_t PassesVZEROcuts(Int_t id);
  Bool_t PassesCuts(const AliFlowTrackSimple* track);
  Bool_t PassesCuts(const AliMultiplicity* track, Int_t id);
  Bool_t PassesCuts(const AliAODTracklets* track, Int_t id);  // XZhang 20120615
  Bool_t PassesCuts(const AliESDkink* kink);
  Bool_t PassesCuts(const AliESDv0* v0);
  Bool_t PassesMCcuts();
  Bool_t PassesMCcuts(AliMCEvent* mcevent, Int_t label);
  Bool_t PassesTPCdedxCut(const AliESDtrack* track);
  Bool_t PassesTPCbayesianCut(const AliESDtrack* track);
  Bool_t PassesTPCpidCut(const AliESDtrack* track) const;
  Bool_t PassesTOFbetaCut(const AliESDtrack* track);
  Bool_t PassesTOFbetaSimpleCut(const AliESDtrack* track);
  Bool_t PassesTOFpidCut(const AliESDtrack* track) const;
  Bool_t PassesESDpidCut(const AliESDtrack* track);
  Bool_t PassesAODpidCut(const AliAODTrack* track);
  Bool_t PassesMuonCuts(AliVParticle* track);  // XZhang 20120604
  Bool_t PassesTPCbayesianCut(const AliAODTrack* track);
  Bool_t PassesTOFbayesianCut(const AliAODTrack* track);
  Bool_t PassesTOFbetaCut(const AliAODTrack* track);
  Int_t  GetITStype(const AliAODTrack* track) const;
  
  void Browse(TBrowser* b);
  Long64_t Merge(TCollection* list);
    
  void SetTPCTOFNsigmaPIDPurityFunctions(Float_t purityLevel);
  void SetCentralityPercentile(Int_t centMin,Int_t centMax){fCentralityPercentileMin=centMin; fCentralityPercentileMax=centMax;}

  void SetPtTOFPIDoff(Double_t pt) {fPtTOFPIDoff = pt;} // added by B.Hohlweger
  Double_t GetPtTOFPIDoff() {return fPtTOFPIDoff;}

  //gain equalization and recentering
  void SetVZEROgainEqualisation(TH1* g) {fVZEROgainEqualization=g;}
  void SetVZEROgainEqualisationCen(TH2* g) {fVZEROgainEqualizationCen=g;}
  void SetVZEROApol(Int_t ring, Float_t f) {fVZEROApol[ring]=f;}
  void SetVZEROCpol(Int_t ring, Float_t f) {fVZEROCpol[ring]=f;}
  // set the flag for recentering (which is done in AliFlowEvent)
  void SetApplyRecentering(Bool_t r)    { fApplyRecentering = r; }
  Bool_t GetApplyRecentering() const    { return fApplyRecentering;}
  void SetVZEROgainEqualizationPerRing(Bool_t s)   {fVZEROgainEqualizationPerRing = s;}
  Bool_t GetVZEROgainEqualizationPerRing() const {return fVZEROgainEqualizationPerRing;}
  void SetDivSigma(Bool_t r)    { fDivSigma = r; }
  Bool_t GetDivSigma() const    { return fDivSigma;}
  // exclude vzero rings: 0 through 7 can be excluded by calling this setter multiple times
  // 0 corresponds to segment ID 0 through 7, etc
  // disabled vzero rings get weight 0
  void SetUseVZERORing(Int_t i, Bool_t u) {
      fUseVZERORing[i] = u;
      fVZEROgainEqualizationPerRing = kTRUE;       // must be true for this option
  }
  Bool_t GetUseVZERORing(Int_t i) const {return fUseVZERORing[i];}
  void SetChi2A(TArrayD*  Chi2A) {fChi2A = Chi2A;}  // chi vs cent for vzero A ep_2
  void SetChi3A(TArrayD*  Chi3A) {fChi3A = Chi3A;}
  void SetChi2C(TArrayD*  Chi2C) {fChi2C = Chi2C;}
  void SetChi3C(TArrayD*  Chi3C) {fChi3C = Chi3C;}

  TArrayD* GetChi2A() {return fChi2A;}  // chi vs cent for vzero A ep_2
  TArrayD* GetChi3A() {return fChi3A;} 
  TArrayD* GetChi2C() {return fChi2C;} 
  TArrayD* GetChi3C() {return fChi3C;} 

  void SetNumberOfSigmas(Float_t val) {fNsigmaCut2 = val*val;};
  Float_t GetNumberOfSigmas() const {return TMath::Sqrt(fNsigmaCut2);};
 
  void SetRun(Int_t const run) {this->fRun = run;};
  Int_t GetRun() const {return this->fRun;};

 protected:
  //AliFlowTrack* MakeFlowTrackSPDtracklet() const;
  //AliFlowTrack* MakeFlowTrackPMDtrack() const;
  //AliFlowTrack* MakeFlowTrackVZERO() const;
  //AliFlowTrack* MakeFlowTrackVParticle() const;
  Bool_t FillFlowTrackVParticle(AliFlowTrack* t) const;
  Bool_t FillFlowTrackGeneric(AliFlowTrack* t) const;
  AliFlowTrack* FillFlowTrackKink(TObjArray* trackCollection, Int_t trackIndex) const;
  AliFlowTrack* FillFlowTrackVZERO(TObjArray* trackCollection, Int_t trackIndex) const;
  AliFlowTrack* FillFlowTrackGeneric(TObjArray* trackCollection, Int_t trackIndex) const;
  AliFlowTrack* FillFlowTrackVParticle(TObjArray* trackCollection, Int_t trackIndex) const;
  void HandleESDtrack(AliESDtrack* track);
  void HandleVParticle(AliVParticle* track);
  void DefineHistograms();
  void InitPIDcuts();
  void InitESDcuts() {if (!fAliESDtrackCuts) {fAliESDtrackCuts=new AliESDtrackCuts();}}
  void InitMuonCuts() { if (!fMuonTrackCuts)  fMuonTrackCuts  =new AliMuonTrackCuts("StdMuCuts","StdMuCuts"); return; }  // XZhang 20120604
  // part added by F. Noferini
  Bool_t PassesTOFbayesianCut(const AliESDtrack* track); 
  Bool_t PassesNucleiSelection(const AliESDtrack* track);   // added by Natasha
  Bool_t PassesTPCTOFNsigmaCut(const AliAODTrack* track); 
  Bool_t PassesTPCTOFNsigmaCut(const AliESDtrack* track);
  Bool_t PassesTPCTOFNsigmaPurityCut(const AliAODTrack* track);
  Bool_t TPCTOFagree(const AliVTrack *track);
  // end part added by F. Noferini
  Bool_t PassesTPCTPCTOFNsigmaCut(const AliAODTrack* track); // added by B. Hohlweger

  //the cuts
  AliESDtrackCuts* fAliESDtrackCuts; //alianalysis cuts
  AliMuonTrackCuts* fMuonTrackCuts;  // muon selection cuts // XZhang 20120604
  TList* fQA;                        //qa histograms go here
  Bool_t fCutMC;                     //do we cut on MC?
  Bool_t fCutMChasTrackReferences;   //did we leave a trace in the detector?
  Bool_t fCutMCprocessType;          //do we cut on mc process type?
  TMCProcess fMCprocessType;         //mc process type
  Bool_t fCutMCPID;                  //cut on MC pid?
  Int_t fMCPID;                      //MC PID
  Bool_t fCutMCfirstMotherPID;       //cut on PID of first mother?
  Int_t fMCfirstMotherPID;           //PID of the first mother of track
  Bool_t fIgnoreSignInMCPID;           //when MC PID cut is set, pass also the antiparticle
  Bool_t fCutMCisPrimary;            //do we cut on primaryness?
  Bool_t fRequireTransportBitForPrimaries; //require the transport bit to be set for primaries
  Bool_t fMCisPrimary;               //is MC primary
  Bool_t fRequireCharge;          //is charged?
  Bool_t fFakesAreOK;             //are fakes (negative labels) ok?
  Bool_t fCutSPDtrackletDeltaPhi; //are we cutting on the trcklet deltaphi?
  Double_t fSPDtrackletDeltaPhiMax; //maximal deltaphi for tracklets
  Double_t fSPDtrackletDeltaPhiMin; //minimal deltaphi for tracklets
  Bool_t fIgnoreTPCzRange;   //ignore tracks going close to central membrane
  Double_t fIgnoreTPCzRangeMax; //max z to ignore
  Double_t fIgnoreTPCzRangeMin; //min z to ignore
  Bool_t fCutChi2PerClusterTPC; //cut on tpc chi2
  Float_t fMaxChi2PerClusterTPC; //max chi2 tpc/cluster
  Float_t fMinChi2PerClusterTPC; //min chi2 tpc/cluster
  Bool_t fCutFracSharedTPCCluster; //cut on fraction of shared TPC clusters
  Float_t fMaxFracSharedTPCCluster; //max fraction of shared TPC clusters
  Bool_t fCutCrossedTPCRows;     //cut on number crossed TPC rows
  Int_t fMinNCrossedRows;        //minimum number of crossed rows
  Float_t fMinCrossedRowsOverFindableClusters; //min. number of crossed rows / findable clusters
  Bool_t fCutGoldenChi2;         //cut on golden chi2 (Chi2TPCConstrainedVsGlobal)
  Float_t fMaxGoldenChi2;        //max golden chi2 (Chi2TPCConstrainedVsGlobal)
  Bool_t fRequireTOFSignal;      //require TOF signal
  Bool_t fCutNClustersTPC;       //cut on clusters?
  Int_t fNClustersTPCMax;        //max tpc ncls
  Bool_t fCutChi2PerClusterITS;  //cut on chi2 per ITS cluster
  Bool_t fCutITSClusterGlobal;   //cut on ITS clusters: either any hit on SPD or no hit on SPD and hit on first layer SDD (like global tracks)
  Float_t fMaxChi2PerClusterITS; //max chi2 per ITS cluster
  Int_t fNClustersTPCMin;        //min tpc clusters  
  Bool_t fCutNClustersITS;       //cut on clusters?
  Int_t fNClustersITSMax;        //max tpc ncls
  Int_t fNClustersITSMin;        //min tpc clusters  
  Bool_t fUseAODFilterBit;       //use AOD filter bit selection?
  UInt_t fAODFilterBit;          //AOD filter bit to select
  Bool_t fCutDCAToVertexXY;      //dca xy cut
  Bool_t fCutDCAToVertexZ;       //dca z cut
  Bool_t fCutDCAToVertexXYPtDepAOD; //dca xy cut pt dep (AOD only)
  Bool_t fCutDCAToVertexXYAOD;   //dca xy cut (AOD only)
  Float_t fMaxDCAxyAOD;          //max dca xy (AOD only)
  Bool_t fCutDCAToVertexZAOD;    //dca z cut (AOD only)
  Float_t fMaxDCAzAOD;           //max dca z (AOD only)
  Bool_t fCutMinimalTPCdedx;     //cut on minimal dedx in TPC to reject noise tracks
  Double_t fMinimalTPCdedx;       //value for minimal TPC dedx
  Bool_t fCutTPCSecbound;         // cut tracks entering TPC close to TPC sector boundaries
  Double_t fCutTPCSecboundMinpt;  // minimum pT for previous cut
  Bool_t fCutTPCSecboundVar;      // cut tracks entering TPC close to TPC sector boundaries
  TF1* fPhiCutLow; //!
  TF1* fPhiCutHigh; //!
  Bool_t fLinearizeVZEROresponse; //linearize VZERO response using AliESDUtil
 
  Int_t fCentralityPercentileMin; //centrality min
  Int_t fCentralityPercentileMax; //centrality max
  Float_t fPurityLevel; //Purity cut percentage

  Bool_t  fCutPmdDet;   //cut on PMD detector plane 
  Int_t   fPmdDet;      // value of PMD detector plane
  Bool_t  fCutPmdAdc;   //cut on cluster ADC
  Float_t fPmdAdc;      //value of cluster ADC
  Bool_t  fCutPmdNcell; //cut on cluster ncell
  Float_t fPmdNcell;    //value of cluster ncell

  Double_t fMinKinkAngle; //max kink angle
  Double_t fMinKinkRadius; //min kink radius
  Double_t fMaxKinkRadius; //max kink radius
  Double_t fMinKinkQt; //min kink qt
  Double_t fMaxKinkQt; //max kink qt
  Double_t fMinKinkInvMassKmu; //max kink inv mass
  Double_t fMaxKinkInvMassKmu; //max kink inv mass
  Bool_t fForceTPCstandalone; //use TPC parameters when applying cuts on the kink mother
  Bool_t fRequireKinkDaughters; //well, the name says it all
   
  trackParameterType fParamType;     //parameter type tu cut on
  trackParameterMix fParamMix;       //parameter mixing
  
  AliESDkink* fKink;                 //!placeholder for the current kink
  AliESDv0* fV0;                     //!placeholder for the current V0
  AliVParticle* fTrack;              //!the track to apply cuts on
  Double_t fTrackMass;               //!mass of the particle
  Double_t fTrackPt;                 //!track pt
  Double_t fTrackPhi;                //!track phi
  Double_t fTrackEta;                //!track eta
  Double_t fTrackWeight;             //!track weight
  Int_t fTrackLabel;                 //!track label, or its absolute value if FakesAreOK
  AliMCEvent* fMCevent;              //!mc event
  AliMCParticle* fMCparticle;        //!mc particle
  AliVEvent* fEvent;                 //!placeholder for current event
  AliESDtrack fTPCtrack;             //!placeholder for TPC only track to avoid new/delete on every track

  //PID
  AliESDpid fESDpid; //pid obj
  AliFlowBayesianPID *fBayesianResponse; //! Baysian response with all the TOF tuning (using fESDpid)
  PIDsource fPIDsource; //pid source
  TMatrixF* fTPCpidCuts; //tpc pid cuts
  TMatrixF* fTOFpidCuts; //tof pid cuts
  AliPID::EParticleType fParticleID; //alipid
  Double_t fParticleProbability; //desired prob for a particle type
  Bool_t fAllowTOFmismatchFlag; //allow TOFmismatch flag=1 in ESD
  Bool_t fRequireStrictTOFTPCagreement; //require stricter than TOFmismatch flag TOF-TPC agreement
  Bool_t fCutRejectElectronsWithTPCpid; //reject electrons with TPC pid

  // part added by F. Noferini
  static const Int_t fgkPIDptBin = 20; // pT bins for priors
  Float_t fC[fgkPIDptBin][5],fBinLimitPID[fgkPIDptBin]; // pt bin limit and priors
  Float_t fProbBayes; // bayesian probability
  Float_t fCurrCentr; // current centrality used for set the priors
  // end part added by F. Noferini
  Double_t fPtTOFPIDoff;
  
  //gain equalization and recentering for vzero
  TH1* fVZEROgainEqualization;     //! equalization histo
  TH2* fVZEROgainEqualizationCen;  //! equalization histo per centrality bin
  Bool_t fApplyRecentering;     // apply recentering of q-sub vectors in AliFlowEvent ?
  Bool_t fVZEROgainEqualizationPerRing;    // per ring vzero gain calibration
  Bool_t fDivSigma;                // divide by st.dev. after recentering
  Float_t fVZEROApol[4];           //! calibration info per ring
  Float_t fVZEROCpol[4];           //! calibration info per ring
  Bool_t fUseVZERORing[8];      // kTRUE means the ring is included
  static const Int_t fgkNumberOfVZEROtracks=64; //number of VZERO channels
  TArrayD*      fChi2A;                 // chi vs cent for vzero A ep_2
  TArrayD*      fChi2C;                 // chi vs cent for vzero C ep_2
  TArrayD*      fChi3A;                 // chi vs cent for vzero A ep_3
  TArrayD*      fChi3C;                 // chi vs cent for vzero C ep_3

  AliPIDResponse *fPIDResponse;            //! Pid reponse to manage Nsigma cuts
  Float_t fNsigmaCut2;                     // Number of sigma^2 (cut value) for TPC+TOF nsigma cut
    
  //TPC TOF nsigma Purity based cut functions
  TFile                 *fPurityFunctionsFile;       //! purity functions file
  TDirectory            *fPurityFunctionsList;     //! purity functions list
    
  TF2                   *fPurityFunction[180]; //TF2 purity functions
  
  Bool_t fCutITSclusterShared;          // cut fMaxITSClusterShared
  Int_t  fMaxITSclusterShared;          // fMaxITSclusterShared 
  Bool_t fCutITSChi2;                   // cut fMaxITSChi2
  Double_t  fMaxITSChi2;                // fMaxITSChi2
  Int_t         fRun;                   // run number
  
  ClassDef(AliFlowTrackCuts,20)
};

#endif



