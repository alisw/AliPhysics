/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliFlowTrackESDCuts:
// A cut class for ESD, AOD and MC particles for the flow framework
// author: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#ifndef ALIFLOWTRACKCUTS_H
#define ALIFLOWTRACKCUTS_H

#include <TMatrix.h>
#include "AliFlowTrackSimpleCuts.h"
#include "AliESDtrackCuts.h"
#include "TMCProcess.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliESDpid.h"

class TObjArray;
class AliVParticle;
class AliMCParticle;
class AliFlowTrack;
class AliMCEvent;
class AliVEvent;
class AliMultiplicity; 

class AliFlowTrackCuts : public AliFlowTrackSimpleCuts {

 public:
  AliFlowTrackCuts();
  AliFlowTrackCuts(const char* name);
  AliFlowTrackCuts(const AliFlowTrackCuts& someCuts);
  AliFlowTrackCuts& operator=(const AliFlowTrackCuts& someCuts);
  virtual ~AliFlowTrackCuts();

  static AliFlowTrackCuts* GetStandardTPCOnlyTrackCuts();
  static AliFlowTrackCuts* GetStandardITSTPCTrackCuts2009(Bool_t selPrimaries=kTRUE);

  enum trackParameterType { kMC, kGlobal, kESD_TPConly, kESD_SPDtracklet };
  enum trackParameterMix  { kPure, kTrackWithMCkine, kTrackWithMCPID, kTrackWithMCpt, kTrackWithPtFromFirstMother };
  enum PIDsource {kTPCpid, kTOFpid, kTPCTOFpid, kTOFbayesian};

  //setters (interface to AliESDtrackCuts)
  void SetMinNClustersTPC( Int_t a ) {InitESDcuts(); fAliESDtrackCuts->SetMinNClustersTPC(a);}
  void SetMinNClustersITS( Int_t a ) {InitESDcuts(); fAliESDtrackCuts->SetMinNClustersITS(a);}
  void SetClusterRequirementITS( AliESDtrackCuts::Detector det,
                                 AliESDtrackCuts::ITSClusterRequirement req = AliESDtrackCuts::kOff )
                                 { InitESDcuts(); fAliESDtrackCuts->SetClusterRequirementITS(det,req); } 
  void SetMaxChi2PerClusterTPC( Float_t a ) {fMaxChi2PerClusterTPC=a;fCutChi2PerClusterTPC=kTRUE;}
  void SetMinChi2PerClusterTPC( Float_t a ) {fMinChi2PerClusterTPC=a;fCutChi2PerClusterTPC=kTRUE;}
  void SetMaxChi2PerClusterITS( Float_t a ) {InitESDcuts(); fAliESDtrackCuts->SetMaxChi2PerClusterITS(a);}
  void SetRequireTPCRefit( Bool_t a ) {InitESDcuts(); fAliESDtrackCuts->SetRequireTPCRefit(a);}
  void SetRequireTPCStandAlone( Bool_t a) {InitESDcuts(); fAliESDtrackCuts->SetRequireTPCStandAlone(a);}
  void SetRequireITSRefit( Bool_t a ) {InitESDcuts(); fAliESDtrackCuts->SetRequireITSRefit(a);}
  void SetRequireITSStandAlone( Bool_t a) {InitESDcuts(); fAliESDtrackCuts->SetRequireITSStandAlone(a);}
  void SetAcceptKinkDaughters( Bool_t a ) {InitESDcuts(); fAliESDtrackCuts->SetAcceptKinkDaughters(a);}
  void SetMaxDCAToVertexZ( Float_t a ) {InitESDcuts(); fAliESDtrackCuts->SetMaxDCAToVertexZ(a);}
  void SetMaxDCAToVertexXY( Float_t a ) {InitESDcuts(); fAliESDtrackCuts->SetMaxDCAToVertexXY(a);}
  void SetMaxDCAToVertexXYPtDep( const char* a ) {InitESDcuts(); fAliESDtrackCuts->SetMaxDCAToVertexXYPtDep(a);}
  void SetRequireSigmaToVertex(Bool_t a) {InitESDcuts(); fAliESDtrackCuts->SetRequireSigmaToVertex(a);}
  void SetMaxNsigmaToVertex(Float_t sigma=1e10) {InitESDcuts(); fAliESDtrackCuts->SetMaxNsigmaToVertex(sigma); }
  void SetDCAToVertex2D( Bool_t a ) {InitESDcuts(); fAliESDtrackCuts->SetDCAToVertex2D(a);}
  void SetEtaRange( Float_t r1, Float_t r2 ) { SetEtaMin(r1); SetEtaMax(r2); }
  void SetPtRange( Float_t r1, Float_t r2 ) { SetPtMin(r1); SetPtMax(r2); }
  void SetRequireCharge( Bool_t r ) {fRequireCharge=r;}
  void SetFakesAreOK( Bool_t b ) {fFakesAreOK=b;}
  void SetSPDtrackletDeltaPhiMax( Double_t m ) {fSPDtrackletDeltaPhiMax=m; fCutSPDtrackletDeltaPhi=kTRUE;}
  void SetSPDtrackletDeltaPhiMin( Double_t m ) {fSPDtrackletDeltaPhiMin=m; fCutSPDtrackletDeltaPhi=kTRUE;}
  void SetIgnoreSignInPID( Bool_t b ) {fIgnoreSignInPID=b;}
  void SetIgnoreTPCzRange( Double_t min, Double_t max ) 
                         { fIgnoreTPCzRange=kTRUE; fIgnoreTPCzRangeMin=min; fIgnoreTPCzRangeMax=max; } 

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

  void SetQA(const char* dirname);
  TObjArray* GetQA() const {return fQA;}

  //MC stuff
  void SetCutMC( Bool_t b=kTRUE );
  void SetMCprocessType( TMCProcess t ) { fMCprocessType = t; fCutMCprocessType=kTRUE; SetCutMC();}
  void SetMCisPrimary( Bool_t b ) { fMCisPrimary=b; fCutMCisPrimary=kTRUE; SetCutMC();}
  void SetMCPID( Int_t pid ) { fMCPID=pid; fCutMCPID=kTRUE; SetCutMC(); }
  TMCProcess GetMCprocessType() const { return fMCprocessType; }
  Bool_t GetMCisPrimary() const {return fMCisPrimary;}
  Int_t GetMCPID() const {return fMCPID;}

  void SetParamType(trackParameterType paramType) {fParamType=paramType;}
  trackParameterType GetParamType() const {return fParamType;}
  static const char* GetParamTypeName(trackParameterType type);
  void SetParamMix(trackParameterMix paramMix) {fParamMix=paramMix;}
  trackParameterMix GetParamMix() const {return fParamMix;}

  virtual Bool_t IsSelected(TObject* obj, Int_t id=-666);
  virtual Bool_t IsSelectedMCtruth(TObject* obj, Int_t id=-666);
  AliVParticle* GetTrack() const {return fTrack;}
  AliMCParticle* GetMCparticle() const {return fMCparticle;}
  AliFlowTrack* MakeFlowTrack(int index=0) const;
  Bool_t IsPhysicalPrimary() const; 
  static Bool_t IsPhysicalPrimary(AliMCEvent* p, Int_t label, Bool_t requiretransported=kTRUE); 
  
  void SetMCevent(AliMCEvent* mcEvent) {fMCevent=mcEvent;}
  AliMCEvent* GetMCevent() const {return fMCevent;}
  void SetEvent(AliVEvent* event, AliMCEvent* mcEvent=NULL);
  AliVEvent* GetEvent() const {return fEvent;}
  Int_t GetNumberOfInputObjects() const;
  TObject* GetInputObject(Int_t i);
  void Clear(Option_t* option="");

  //PID
  void SetPID(AliPID::EParticleType pid, PIDsource s=kTPCTOFpid) {fAliPID=pid; fPIDsource=s; fCutPID=kTRUE; InitPIDcuts();}
  void SetTPCTOFpidCrossOverPt(Double_t pt) {fTPCTOFpidCrossOverPt=pt;}
  void SetTPCpidCuts(TMatrixF* mat) {fTPCpidCuts=new TMatrixF(*mat);}
  void SetTOFpidCuts(TMatrixF* mat) {fTOFpidCuts=new TMatrixF(*mat);}
  AliESDpid& GetESDpid() {return fESDpid;}

 protected:
  Bool_t PassesCuts(AliVParticle* track);
  Bool_t PassesCuts(AliFlowTrackSimple* track);
  Bool_t PassesCuts(AliMultiplicity* track, Int_t id);
  Bool_t PassesMCcuts();
  Bool_t PassesMCcuts(AliMCEvent* mcevent, Int_t label);
  Bool_t PassesTPCpidCut(AliESDtrack* track);
  Bool_t PassesTOFpidCut(AliESDtrack* track);  
  void HandleESDtrack(AliESDtrack* track);
  void HandleVParticle(AliVParticle* track);
  void DefineHistograms();
  void InitPIDcuts();
  void InitESDcuts() {if (!fAliESDtrackCuts) fAliESDtrackCuts=new AliESDtrackCuts();}
  // part added by F. Noferini
  Bool_t PassesTOFbayesianCut(AliESDtrack* track); 
  void SetPriors(); // set my favourite priors
  Int_t GetESDPdg(AliESDtrack *track,Option_t *option="bayesianTOF",Int_t ipart=2,Float_t cPi=-1.0,Float_t cKa=0.0,Float_t cPr=0.0); // 3sigma cut ipart=0(el),1(mu),2(pi),3(K),4(p)
  // end part added by F. Noferini

  //the cuts
  AliESDtrackCuts* fAliESDtrackCuts; //alianalysis cuts
  TObjArray* fQA;                   //qa histograms go here
  Bool_t fCutMC;                     //do we cut on MC?
  Bool_t fCutMCprocessType;          //do we cut on mc process type?
  TMCProcess fMCprocessType;         //mc process type
  Bool_t fCutMCPID;                  //cut on MC pid?
  Int_t fMCPID;                      //MC PID
  Bool_t fIgnoreSignInPID;           //when PID cut is set, pass also the antiparticle
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
  Bool_t fCutNClustersTPC;       //cut on clusters?
  Int_t fNClustersTPCMax;        //max tpc ncls
  Int_t fNClustersTPCMin;        //min tpc clusters  

  trackParameterType fParamType;     //parameter type tu cut on
  trackParameterMix fParamMix;       //parameter mixing
  AliVParticle* fTrack;              //!the track to apply cuts on
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
  PIDsource fPIDsource; //pid source
  TMatrixF* fTPCpidCuts; //tpc pid cuts
  TMatrixF* fTOFpidCuts; //tof pid cuts
  Double_t fTPCTOFpidCrossOverPt; //pt cross over for pid, below TPC is taken, above TOF
  AliPID::EParticleType fAliPID; //alipid

  // part added by F. Noferini
  static const Int_t fnPIDptBin = 20; // pT bins for priors
  Float_t fC[fnPIDptBin][5],fBinLimitPID[fnPIDptBin]; // pt bin limit and priors
  Float_t fProbBayes[5]; // bayesian probability
  // end part added by F. Noferini

  ClassDef(AliFlowTrackCuts,4)
};

#endif


