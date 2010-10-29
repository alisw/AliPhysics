/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliFlowTrackESDCuts:
// A cut class for ESD, AOD and MC particles for the flow framework
// author: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#ifndef ALIFLOWTRACKCUTS_H
#define ALIFLOWTRACKCUTS_H

#include "AliFlowTrackSimpleCuts.h"
#include "AliESDtrackCuts.h"
#include "TMCProcess.h"

class TDirectory;
class AliVParticle;
class AliMCParticle;
class AliFlowTrack;
class AliMCEvent;
class AliVEvent;
class AliMultiplicity; 

class AliFlowTrackCuts : public AliFlowTrackSimpleCuts {

 public:
  AliFlowTrackCuts();
  AliFlowTrackCuts(const AliFlowTrackCuts& someCuts);
  AliFlowTrackCuts& operator=(const AliFlowTrackCuts& someCuts);
  virtual ~AliFlowTrackCuts();

  static AliFlowTrackCuts* GetStandardTPCOnlyTrackCuts();
  static AliFlowTrackCuts* GetStandardITSTPCTrackCuts2009(Bool_t selPrimaries=kTRUE);

  enum trackParameterType { kMC, kGlobal, kESD_TPConly, kESD_SPDtracklet };
  enum trackParameterMix  { kPure, kTrackWithMCkine, kTrackWithMCPID, kTrackWithMCpt };

  //setters (interface to AliESDtrackCuts)
  void SetMinNClustersTPC( Int_t a ) {fAliESDtrackCuts->SetMinNClustersTPC(a);}
  void SetMinNClustersITS( Int_t a ) {fAliESDtrackCuts->SetMinNClustersITS(a);}
  void SetClusterRequirementITS( AliESDtrackCuts::Detector det,
                                 AliESDtrackCuts::ITSClusterRequirement req = AliESDtrackCuts::kOff )
                                 { fAliESDtrackCuts->SetClusterRequirementITS(det,req); } 
  void SetMaxChi2PerClusterTPC( Float_t a ) {fAliESDtrackCuts->SetMaxChi2PerClusterTPC(a);}
  void SetMaxChi2PerClusterITS( Float_t a ) {fAliESDtrackCuts->SetMaxChi2PerClusterITS(a);}
  void SetRequireTPCRefit( Bool_t a ) {fAliESDtrackCuts->SetRequireTPCRefit(a);}
  void SetRequireTPCStandAlone( Bool_t a) {fAliESDtrackCuts->SetRequireTPCStandAlone(a);}
  void SetRequireITSRefit( Bool_t a ) {fAliESDtrackCuts->SetRequireITSRefit(a);}
  void SetRequireITSStandAlone( Bool_t a) {fAliESDtrackCuts->SetRequireITSStandAlone(a);}
  void SetAcceptKinkDaughters( Bool_t a ) {fAliESDtrackCuts->SetAcceptKinkDaughters(a);}
  void SetMaxDCAToVertexZ( Float_t a ) {fAliESDtrackCuts->SetMaxDCAToVertexZ(a);}
  void SetMaxDCAToVertexXY( Float_t a ) {fAliESDtrackCuts->SetMaxDCAToVertexXY(a);}
  void SetMaxDCAToVertexXYPtDep( const char* a ) {fAliESDtrackCuts->SetMaxDCAToVertexXYPtDep(a);}
  void SetRequireSigmaToVertex(Bool_t a) {fAliESDtrackCuts->SetRequireSigmaToVertex(a);}
  void SetMaxNsigmaToVertex(Float_t sigma=1e10) { fAliESDtrackCuts->SetMaxNsigmaToVertex(sigma); }
  void SetDCAToVertex2D( Bool_t a ) {fAliESDtrackCuts->SetDCAToVertex2D(a);}
  void SetEtaRange( Float_t r1, Float_t r2 ) { SetEtaMin(r1); SetEtaMax(r2); }
  void SetPtRange( Float_t r1, Float_t r2 ) { SetPtMin(r1); SetPtMax(r2); }
  void SetRequireCharge( Bool_t r ) {fRequireCharge=r;}
  void SetFakesAreOK( Bool_t b ) {fFakesAreOK=b;}
  void SetSPDtrackletDeltaPhiMax( Double_t m ) {fSPDtrackletDeltaPhiMax=m; fCutSPDtrackletDeltaPhi=kTRUE;}
  void SetSPDtrackletDeltaPhiMin( Double_t m ) {fSPDtrackletDeltaPhiMin=m; fCutSPDtrackletDeltaPhi=kTRUE;}

  Int_t GetMinNClustersTPC() const {return fAliESDtrackCuts->GetMinNClusterTPC();}
  Int_t GetMinNClustersITS() const {return fAliESDtrackCuts->GetMinNClustersITS();}
  AliESDtrackCuts::ITSClusterRequirement GetClusterRequirementITS( AliESDtrackCuts::Detector det ) const
                                 { return fAliESDtrackCuts->GetClusterRequirementITS(det); } 
  Float_t GetMaxChi2PerClusterTPC() const {return fAliESDtrackCuts->GetMaxChi2PerClusterTPC();}
  Float_t GetMaxChi2PerClusterITS() const {return fAliESDtrackCuts->GetMaxChi2PerClusterITS();}
  Bool_t GetRequireTPCRefit() const {return fAliESDtrackCuts->GetRequireTPCRefit();}
  Bool_t GetRequireTPCStandAlone() const {return fAliESDtrackCuts->GetRequireTPCStandAlone();}
  Bool_t GetRequireITSRefit() const {return fAliESDtrackCuts->GetRequireITSRefit();}
  Bool_t GetRequireITSStandAlone() const {return fAliESDtrackCuts->GetRequireITSStandAlone();}
  Bool_t GetAcceptKinkDaughters() const {return fAliESDtrackCuts->GetAcceptKinkDaughters();}
  Float_t GetMaxDCAToVertexZ() const {return fAliESDtrackCuts->GetMaxDCAToVertexZ();}
  Float_t GetMaxDCAToVertexXY() const {return fAliESDtrackCuts->GetMaxDCAToVertexXY();}
  const char* GetMaxDCAToVertexXYPtDep() const {return fAliESDtrackCuts->GetMaxDCAToVertexXYPtDep();}
  Bool_t GetRequireSigmaToVertex() const {return fAliESDtrackCuts->GetRequireSigmaToVertex();}
  Float_t GetMaxNsigmaToVertex() const {return fAliESDtrackCuts->GetMaxNsigmaToVertex(); }
  Bool_t GetDCAToVertex2D() const {return fAliESDtrackCuts->GetDCAToVertex2D();}
  void GetEtaRange( Float_t& r1, Float_t& r2 ) const { r1=GetEtaMin(); r2=GetEtaMax(); }
  void GetPtRange( Float_t& r1, Float_t& r2 ) const { r1=GetPtMin(); r2=GetPtMax(); }
  Bool_t GetRequireCharge() const {return fRequireCharge;}
  Bool_t GetFakesAreOK() const {return fFakesAreOK;}
  Double_t GetSPDtrackletDeltaPhiMax() const {return fSPDtrackletDeltaPhiMax;}
  Double_t GetSPDtrackletDeltaPhiMin() const {return fSPDtrackletDeltaPhiMin;}

  void SetQA(const char* dirname);
  TDirectory* GetQA() const {return fQA;}

  //MC stuff
  void SetMCprocessType( TMCProcess t ) { fMCprocessType = t; fCutMCprocessType=kTRUE; }
  TMCProcess GetMCprocessType() const { return fMCprocessType; }
  void SetMCisPrimary( Bool_t b ) { fMCisPrimary=b; fCutMCisPrimary=kTRUE; }
  Bool_t GetMCisPrimary() const {return fMCisPrimary;}

  void SetParamType(trackParameterType paramType) {fParamType=paramType;}
  trackParameterType GetParamType() const {return fParamType;}
  static const char* GetParamTypeName(trackParameterType type);
  void SetParamMix(trackParameterMix paramMix) {fParamMix=paramMix;}
  trackParameterMix GetParamMix() const {return fParamMix;}

  virtual Bool_t IsSelected(TObject* obj, Int_t id=-666);
  AliVParticle* GetTrack() const {return fTrack;}
  AliMCParticle* GetMCparticle() const {return fMCparticle;}
  AliFlowTrack* MakeFlowTrack() const;
  Bool_t IsPhysicalPrimary() const; 
  static Bool_t IsPhysicalPrimary(AliMCEvent* p, Int_t label); 
  
  void SetMCevent(AliMCEvent* mcEvent) {fMCevent=mcEvent;}
  AliMCEvent* GetMCevent() const {return fMCevent;}
  void SetEvent(AliVEvent* event) {fEvent=event;}
  AliVEvent* GetEvent() const {return fEvent;}
  Int_t GetNumberOfInputObjects() const;
  TObject* GetInputObject(Int_t i);

 protected:
  Bool_t PassesCuts(AliVParticle* track);
  Bool_t PassesCuts(AliFlowTrackSimple* track);
  Bool_t PassesCuts(AliMultiplicity* track, Int_t id);
  Bool_t PassesMCcuts();
  void HandleESDtrack(AliESDtrack* track);
  void HandleVParticle(AliVParticle* track);
  void DefineHistograms();

  //the cuts
  AliESDtrackCuts* fAliESDtrackCuts; //alianalysis cuts
  TDirectory* fQA;                   //qa histograms go here
  Bool_t fCutMCprocessType;          //do we cut on mc process type?
  TMCProcess fMCprocessType;         //mc process type
  Bool_t fCutMCPID;                  //cut on MC pid?
  Int_t fMCPID;                      //MC PID
  Bool_t fIgnoreSignInPID;           //when PID cut is set, pass also the antiparticle
  Bool_t fCutMCisPrimary;            //do we cut on primaryness?
  Bool_t fMCisPrimary;               //is MC primary
  Bool_t fRequireCharge;          //is charged?
  Bool_t fFakesAreOK;             //are fakes (negative labels) ok?
  Bool_t fCutSPDtrackletDeltaPhi; //are we cutting on the trcklet deltaphi?
  Double_t fSPDtrackletDeltaPhiMax; //maximal deltaphi for tracklets
  Double_t fSPDtrackletDeltaPhiMin; //minimal deltaphi for tracklets

  trackParameterType fParamType;     //parameter type tu cut on
  trackParameterMix fParamMix;       //parameter mixing
  Bool_t fCleanupTrack;              //check if we need to delete the track
  AliVParticle* fTrack;              //!the track to apply cuts on
  Double_t fTrackPhi;                //!track phi
  Double_t fTrackEta;                //!track eta
  Double_t fTrackWeight;             //!track weight
  Int_t fTrackLabel;                 //!track label, or its absolute value if FakesAreOK
  AliMCEvent* fMCevent;              //!mc event
  AliMCParticle* fMCparticle;        //!mc particle
  AliVEvent* fEvent;                 //!placeholder for current event

  ClassDef(AliFlowTrackCuts,3)
};

#endif


