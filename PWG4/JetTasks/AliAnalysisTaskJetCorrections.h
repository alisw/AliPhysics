#ifndef ALIANALYSISTASKJETCORRECTIONS_H
#define ALIANALYSISTASKJETCORRECTIONS_H

#include "AliAnalysisTaskSE.h"

//
// corrections to jet energy by sona 
// 

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
 
  Int_t fGlobVar; //
  Double_t fXsection; //


  TH1F * fhEGen; //
  TH1F * fhERec; //
  TH1F * fhEGenRest; //
  TH1F * fhERecRest; //
  TH1F * fhEsumGenRest; //
  TH1F * fhEsumRecRest; //

  TH2F * fhE2vsE1Gen; //
  TH2F * fhE2vsE1Rec; //
  TH2F * fhE2E1vsEsumGen; //
  TH2F * fhE2E1vsEsumRec; //
  TH2F * fhE2E1vsE1Gen; //
  TH2F * fhE2E1vsE1Rec; //
  TH2F * fhE2E1vsdPhiGen; //
  TH2F * fhE2E1vsdPhiRec; //

  TH2F * fhTrackBalance2; //
  TH2F * fhTrackBalance3; //

  TH2F * fhEt1Et22; //
  TH2F * fhEt1Et23; //

  TProfile * fhECorrJet10[3]; //
  TProfile * fhECorrJet05[3]; //
  TProfile * fhECorrJet01[3]; //
  TProfile * fhECorrJet001[3]; //
  
  TProfile * fhdEvsErec10[3]; //
  TProfile * fhdEvsErec05[3]; //
  TProfile * fhdEvsErec01[3]; //
  TProfile * fhdEvsErec001[3]; //

  TH2F * fhdPhidEta10[3]; //
  TH2F * fhdPhidEta05[3]; //
  TH2F * fhdPhidEta01[3]; //
  TH2F * fhdPhidEta001[3]; //

  TH2F * fhdPhidEtaPt10[3]; //
  TH2F * fhdPhidEtaPt05[3]; //
  TH2F * fhdPhidEtaPt01[3]; //
  TH2F * fhdPhidEtaPt001[3]; //

  ClassDef(AliAnalysisTaskJetCorrections, 1)
};

#endif
