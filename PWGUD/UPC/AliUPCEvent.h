#ifndef ALIUPCEVENT_H
#define ALIUPCEVENT_H

//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//_____________________________________________________________________________

#define NTRG 24
#include "TObject.h"

class AliUPCEvent : public TObject
{
public:
  AliUPCEvent();
  virtual ~AliUPCEvent();

  void ClearEvent(void);

  //Setters
  void SetIsESD(Bool_t isESD=kTRUE);
  void SetIsMC(Bool_t isMC=kTRUE);
  void SetFlagBit(UChar_t ibit) {fFlags |= (1 << ibit);}
  void ResetFlagBits(void) {fFlags = 0;}

  void SetL0inputs(UInt_t inputs) {fL0inputs=inputs;}
  void SetTriggerClass(Int_t idx, Bool_t fired);
  void SetRecoPass(Char_t passId) {fRecoPass=passId;}
  void SetInputFileName(const char *s);
  void SetEventNumber(Long64_t evt) {fEvtNum=evt;}
  void SetRunNumber(Int_t run) {fRunNum=run;}

  void SetPrimaryVertex(Double_t *pos, Double_t chi2, Double_t *cov, Int_t ncontrib);
  void SetPrimaryVertexTitle(const char *s);
  void SetPrimaryVertexSPD(Double_t *pos, Double_t chi2, Double_t *cov, Int_t ncontrib);
  void SetPrimaryVertexSPDtitle(const char *s);
  void SetPrimaryVertexMC(Float_t *pos);

  void SetNumberOfTracklets(Int_t tcklets) {fNTracklets=tcklets;}
  void SetNSPDfiredInner(Int_t nInner) {fNSPDfiredInner=nInner;}
  void SetNSPDfiredOuter(Int_t nOuter) {fNSPDfiredOuter=nOuter;}

  void SetV0ADecision(Int_t decision) {fV0ADecision=decision;}
  void SetV0CDecision(Int_t decision) {fV0CDecision=decision;}
  void SetBBtriggerV0Cmask(UInt_t ibit);
  void SetBBFlagV0Cmask(UInt_t ibit);

  void SetZNCEnergy(Double_t energy) {fZNCEnergy = energy;}
  void SetZPCEnergy(Double_t energy) {fZPCEnergy = energy;}
  void SetZNAEnergy(Double_t energy) {fZNAEnergy = energy;}
  void SetZPAEnergy(Double_t energy) {fZPAEnergy = energy;}
  void SetZNCtdc(Bool_t tdc) {fZNCtdc = tdc;}
  void SetZPCtdc(Bool_t tdc) {fZPCtdc = tdc;}
  void SetZNAtdc(Bool_t tdc) {fZNAtdc = tdc;}
  void SetZPAtdc(Bool_t tdc) {fZPAtdc = tdc;}

  AliUPCTrack *AddTrack(void);
  AliUPCMuonTrack *AddMuonTrack(void);

  TParticle *AddMCParticle(void);

  //Getters
  Bool_t GetIsESD(void) const;
  Bool_t GetIsMC(void) const;
  Bool_t GetFlagBit(UChar_t ibit) const;
  UChar_t GetFlag(void) const { return fFlags; }

  UInt_t GetL0inputs(void) const { return fL0inputs; }
  Bool_t GetTriggerClass(Int_t idx) const;
  Char_t GetRecoPass(void) const { return fRecoPass; }
  TObjString *GetInputFileName(void) const { return fDataFilnam; }
  Long64_t GetEventNumber(void) const { return fEvtNum; }
  Int_t GetRunNumber(void) const { return fRunNum; }

  void GetPrimaryVertex(Double_t *pos, Double_t &chi2, Double_t *cov, Int_t &ncontrib) const;
  TObjString *GetPrimaryVertexTitle(void) const { return fVtxTitle; }
  void GetPrimaryVertexSPD(Double_t *pos, Double_t &chi2, Double_t *cov, Int_t &ncontrib) const;
  TObjString *GetPrimaryVertexSPDtitle(void) const { return fVtxSPDtitle; }
  void GetPrimaryVertexMC(Float_t *pos) const;

  Int_t GetNumberOfTracklets(void) const { return fNTracklets; }
  Int_t GetNSPDfiredInner(void) const { return fNSPDfiredInner; }
  Int_t GetNSPDfiredOuter(void) const { return fNSPDfiredOuter; }

  Int_t GetV0ADecision(void) const { return fV0ADecision; }
  Int_t GetV0CDecision(void) const { return fV0CDecision; }
  UInt_t GetBBtriggerV0C(void) const { return fBBtriggerV0C; }
  UInt_t GetBBFlagV0C(void) const { return fBBFlagV0C; }
  Int_t GetNV0ChitsOffline(void) const;
  Int_t GetNV0ChitsOnline(void) const;

  Double_t GetZNCEnergy(void) const { return fZNCEnergy; }
  Double_t GetZPCEnergy(void) const { return fZPCEnergy; }
  Double_t GetZNAEnergy(void) const { return fZNAEnergy; }
  Double_t GetZPAEnergy(void) const { return fZPAEnergy; }
  Bool_t GetZNCtdc(void) const { return fZNCtdc; }
  Bool_t GetZPCtdc(void) const { return fZPCtdc; }
  Bool_t GetZNAtdc(void) const { return fZNAtdc; }
  Bool_t GetZPAtdc(void) const { return fZPAtdc; }
  Bool_t GetAllZDCtdc(void) const;

  Int_t GetNumberOfTracks(void) const { return fNtracks; }
  AliUPCTrack *GetTrack(Int_t iTrack) const;
  Int_t GetNumberOfMuonTracks(void) const { return fNmuons; }
  AliUPCMuonTrack *GetMuonTrack(Int_t iTrack) const;

  Int_t GetNumberOfMCParticles(void) const { return fNmc; }
  TParticle *GetMCParticle(Int_t iMC) const;

protected:
  AliUPCEvent(const AliUPCEvent &o); // not implemented
  AliUPCEvent &operator=(const AliUPCEvent &o); // not implemented

  UChar_t fFlags; // event flags bits: 0 = ESD, 1 = MC
  UInt_t fL0inputs; // L0 trigger inputs
  Bool_t fTrgClasses[NTRG]; // fired trigger classes
  Char_t fRecoPass; // reconstruction pass identifier
  TObjString *fDataFilnam; //-> input file name and path
  Long64_t fEvtNum; // event number in input file
  Int_t fRunNum; // run number
  Double_t fVtxPos[3]; // default primary vertex position
  Double_t fVtxChi2perNDF; // chi2/ndf of vertex fit
  Double_t fVtxCov[6]; // vertex covariance matrix
  Int_t fVtxNContributors; // # of tracklets/tracks used for the estimate
  TObjString *fVtxTitle; //-> title of default primary vertex
  Double_t fVtxSPDpos[3]; // SPD primary vertex position
  Double_t fVtxSPDchi2perNDF; // chi2/ndf of vertex fit
  Double_t fVtxSPDcov[6]; // vertex covariance matrix
  Int_t fVtxSPDnContributors; // # of tracklets/tracks used for the estimate
  TObjString *fVtxSPDtitle; //-> title of SPD primary vertex
  Float_t fVtxMCpos[3]; // MC primary vertex position
  Int_t fNTracklets; // number of SPD tracklets
  Int_t fNSPDfiredInner; // number of fired SPD FO chips, inner layer
  Int_t fNSPDfiredOuter; // number of fired SPD FO chips, outer layer
  Int_t fV0ADecision; // V0A decision, set by enumeration: kV0Invalid = -1, kV0Empty = 0, kV0BB, kV0BG, kV0Fake
  Int_t fV0CDecision; // V0C decision
  UInt_t fBBtriggerV0C; // offline beam-beam flags in V0C one bit per cell
  UInt_t fBBFlagV0C; // online beam-beam flags in V0C one bit per cell
  Double_t fZNCEnergy; // reconstructed energy in the neutron ZDC, C-side
  Double_t fZPCEnergy; // reconstructed energy in the proton ZDC, C-side
  Double_t fZNAEnergy; // reconstructed energy in the neutron ZDC, A-side
  Double_t fZPAEnergy; // reconstructed energy in the proton ZDC, A-side
  Bool_t fZNCtdc; // ZDC TDC data, NC
  Bool_t fZPCtdc; // ZDC TDC data, PC
  Bool_t fZNAtdc; // ZDC TDC data, NA
  Bool_t fZPAtdc; // ZDC TDC data, PA
  TClonesArray *fUPCTracks; //-> array of central upc tracks
  Int_t fNtracks; // number of central upc tracks in event
  TClonesArray *fUPCMuonTracks; //-> array of muon upc tracks
  Int_t fNmuons; // number of muon upc tracks in event
  TClonesArray *fMCParticles; // array of MC particles
  Int_t fNmc; // number of mc particles in event

  static TClonesArray *fgUPCTracks; // array of central upc tracks
  static TClonesArray *fgUPCMuonTracks; // array of muon upc tracks
  static TClonesArray *fgMCParticles; // array of MC particles

  ClassDef(AliUPCEvent,1)
};

#endif





















