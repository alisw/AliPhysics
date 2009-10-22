#ifndef ALIANALYSISTASKJETCORRECTIONS_H
#define ALIANALYSISTASKJETCORRECTIONS_H

#include "AliAnalysisTaskSE.h"
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// ***
// corrections to jet energy by sona.pochybova@cern.ch 
// ***

class AliJetFinder;
class AliESDEvent;
class AliAODEvent;
class AliAODJet;
class AliGenPythiaEventHeader;
class AliAODPid;

class TList;
class TArrayD;
class TChain;
class TH1;
class TH2;
class TH1F;
class TH2F;
class TH2I;
class TH3D;
class TTree;
class TProfile;
class TLorentzVector;
class TVector3;
class TParticle;

class AliAnalysisTaskJetCorrections : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskJetCorrections();
  AliAnalysisTaskJetCorrections(const char * name);
  virtual ~AliAnalysisTaskJetCorrections() {;}
  
  //Implementation of interface methods
  virtual Bool_t Notify(); 
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() { Init(); };
  virtual void UserExec(Option_t * option);
  virtual void Terminate(Option_t * option);

  virtual void SetAODInput(Bool_t b){fUseAODInput = b;}

  virtual void SetBranchGen(const char* c){fBranchGen = c;}
  virtual void SetBranchRec(const char* c){fBranchRec = c;}

  virtual Double_t SetR(Double_t b){fR = b; return fR;} 

  virtual void GetThrustAxis(TVector3 &n01, TVector3 * pTrack,const Int_t &nTracks);
 private:
  AliAnalysisTaskJetCorrections(const AliAnalysisTaskJetCorrections&);
  AliAnalysisTaskJetCorrections& operator = (const AliAnalysisTaskJetCorrections&);

  enum {kMaxJets = 6};
  enum {kMaxEvents = 10};
  enum {kJets = 3};
  enum {kTracks = 1000};

  AliAODEvent  *fAOD; // where we take the jets from can be input or output AOD
  
  TString       fBranchRec;  // AOD branch name for reconstructed
  TString       fBranchGen;  // AOD brnach for genereated
  
  Bool_t        fUseAODInput; // use aod input not output
  Double_t fR; // radius
  TList * fList; // output list
 
  Int_t fGlobVar; // switch enabling checking out just one event
  Double_t fXsection; // cross-section 


  TH1F * fhEGen; // generated energy histogram
  TH1F * fhERec; // reconstructed energy histogram
  TH1F * fhEGenRest; // generated energy in the rest frame of three-jet events
  TH1F * fhERecRest; // generated energy in the rest frame of three-jet events
  TH1F * fhEsumGenRest; //generated summed energy in the rest frame of three-jet events
  TH1F * fhEsumRecRest; // reconstructed summed energy in the rest frame of three-jet events

  TH2F * fhE2vsE1Gen; // leading vs next-to leading jet energy, generated
  TH2F * fhE2vsE1Rec; // leading vs next-to leading jet energy, reconstruted
  TH2F * fhE2E1vsEsumGen; // E2,E1 diference vs. summed energy, generated
  TH2F * fhE2E1vsEsumRec; // E2, E1 difference vs summed jet energy, reconstructed
  TH2F * fhE2E1vsE1Gen; // E2, E1 difference vs E1, generated
  TH2F * fhE2E1vsE1Rec; // E2, E1 diference vs E1, reeconstructed
  TH2F * fhE2E1vsdPhiGen; // E2, E1 difference vs dPhi, generated
  TH2F * fhE2E1vsdPhiRec; // E2, E1 difference vs dPhi, reconstrted

  TH2F * fhTrackBalance2; // track balance in 2-jet events
  TH2F * fhTrackBalance3; // track balance in 3-jet events

  TH2F * fhEt1Et22; // Et1,2 in back-to-back cones in 2-jet events
  TH2F * fhEt1Et23; // Et1,2 in back-to-back cones in 3-jet events

  TProfile * fhECorrJet10[3]; // corr. factor for jets matched within dR = 1.
  TProfile * fhECorrJet05[3]; // corr. factor for jets matched within dR = .5
  TProfile * fhECorrJet01[3]; // corr. factor for jets matched within dR = .1
  TProfile * fhECorrJet001[3]; // corr. factor for jets matched within dR = .01
  
  TProfile * fhdEvsErec10[3]; // energy difference vs rec. energy, dR = 1.
  TProfile * fhdEvsErec05[3]; // energy difference vs rec. energy, dR = .5
  TProfile * fhdEvsErec01[3]; // energy difference vs rec. energy, dR = .1
  TProfile * fhdEvsErec001[3]; // energy difference vs rec. energy, dR = .01

  TH2F * fhdPhidEta10[3]; // dPhi vs dEta of particles, dR = 1.
  TH2F * fhdPhidEta05[3]; // dPhi vs dEta of particles, dR = .5
  TH2F * fhdPhidEta01[3]; // dPhi vs dEta of particles, dR = .1
  TH2F * fhdPhidEta001[3]; // dPhi vs dEta of particles, dR = .01

  TH2F * fhdPhidEtaPt10[3]; // Pt weighted dPhi vs dEta of particles, dR = 1.
  TH2F * fhdPhidEtaPt05[3]; // Pt weighted dPhi vs dEta of particles, dR = .5
  TH2F * fhdPhidEtaPt01[3]; // Pt weighted dPhi vs dEta of particles, dR = .1
  TH2F * fhdPhidEtaPt001[3]; // Pt weighted dPhi vs dEta of particles, dR = .01
  ClassDef(AliAnalysisTaskJetCorrections, 1)
};

#endif
