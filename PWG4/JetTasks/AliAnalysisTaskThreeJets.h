#ifndef ALIANALYSISTASKTHREEJETS_H
#define ALIANALYSISTASKTHREEJETS_H

#include "AliAnalysisTaskSE.h"

// 
// Task for three jet ana by sona 
// 

//class AliJetHistos;
//class AliJetCorrHistos;
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

class AliAnalysisTaskThreeJets : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskThreeJets();
  AliAnalysisTaskThreeJets(const char * name);
  virtual ~AliAnalysisTaskThreeJets() {;}
  
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

  virtual void FillTopology(TH2F * Dalitz, TH1F * fhMu34, TH1F * fhMu45, TH1F * fhMu35, Double_t * x, TVector3 * pRest, Double_t xsection);

  virtual Double_t SetR(Double_t b){fR = b; return fR;} 
  virtual Bool_t IsPrimChar(TParticle* aParticle, Int_t aTotalPrimaries, Bool_t adebug);

  Double_t Exponent(Double_t x,const Double_t * const par) const;
  Double_t Exponent2(Double_t x,const Double_t * const par) const;
  Double_t Gauss(Double_t x,const Double_t * const par) const;
  Double_t Total(Double_t x,const Double_t * const par) const;

 private:
  AliAnalysisTaskThreeJets(const AliAnalysisTaskThreeJets&);
  AliAnalysisTaskThreeJets& operator = (const AliAnalysisTaskThreeJets&);

  enum {kMaxJets = 6};
  enum {kMaxEvents = 10};
  enum {kJets = 3};
  enum {kTracks = 1000};

  AliAODEvent  *fAOD; // where we take the jets from can be input or output AOD
  
  TString       fBranchRec;  // AOD branch name for reconstructed
  TString       fBranchGen;  // AOD brnach for genereated
  
  Bool_t        fUseAODInput; // read the AOD from the input no from the output

  Double_t fR; // radius
  TList * fList; // output list

  Int_t fGlobVar; // globvar
  Double_t fXsection; // xsectio

TH2F * fhX3X4Rec; // Dalitz variables, reconstructed
  TH2F * fhX3X4Gen; // Dalitz variables, generated
  
  TH1F * fhMu35Rec; // scaled masses, 35, reconstructed
  TH1F * fhMu34Rec; // scaled masses, 34, reconstructed
  TH1F * fhMu45Rec; // scaled masses, 45, reconstructed
  
  TH1F * fhMu35Gen; // scaled masses, 35, generated
  TH1F * fhMu34Gen; // scaled masses, 34, generated
  TH1F * fhMu45Gen; // scaled masses, 45, generated

  TH1I * fhInOut; // number of gen. vs number of rec.
  TH1F * fhThrustRec2; // thrust for reco 2-jet events
  TH1F * fhThrustRec3; // thrust for reco 3-jet events

  TH1F * fhThrustGen2; // thrust for gen 2-jet events
  TH1F * fhThrustGen3; // thrust for gen 3-jet events

  TH1F * fhCGen2; // C-variable for gen 2-jets
  TH1F * fhCGen3; // C-variable for gen 3-jets

  TH1F * fhSGen2; // Sphericity for gen 2-jets
  TH1F * fhSGen3; // Sphericity for gen 3-jets

  TH1F * fhAGen2; // A-variable for gen 2-jets
TH1F * fhAGen3; // A-variable for gen 3-jets

  TH1F * fhCRec2; // C-variable for reco 2-jets
  TH1F * fhCRec3; // C-variable for reco 3-jets

  TH1F * fhSRec2; // Sphericity for reco 2-jets
  TH1F * fhSRec3; // Sphericity for reco 3-jets

  TH1F * fhARec2; // A-variable for reco 2-jets
  TH1F * fhARec3; // A-variable for reco 3-jets

  TH2F * fhX3; // dX3 vs X3 rec
  TH2F * fhX4; // dX4 vs X4 rec
  TH2F * fhX5; // dX5 vs X5 rec

  TProfile * fhXSec; // cross-section vs PtHard
  TH2F * fhX3X4Rec60; // Dalitz plane for Esum < 60, reco
  TH2F * fhX3X4Rec60100; // Dalitz plane for 60 < Esum < 100, reco
  TH2F * fhX3X4Rec100; // Dalitz plane for Esum > 100, reco^M 
  TH2F * fhX3X4Gen60; // Dalitz plane for Esum < 60, gen
  TH2F * fhX3X4Gen60100; // Dalitz plane for 60 < Esum < 100, gen
  TH2F * fhX3X4Gen100; // Dalitz plane for Esum > 100, gen

  TH2F * fhdPhiThrustGen; // energy distribution with rspct to thrust axis, gen, 2-jets
  TH2F * fhdPhiThrustGenALL; // energy distribution with rspct to thrust axis, gen 3-jets
  TH2F * fhdPhiThrustRec; // energy distribution with rspct to thrust axis, reco, 2-jets
  TH2F * fhdPhiThrustRecALL; // energy distribution with rspct to thrust axis, reco, 3-jets

  ClassDef(AliAnalysisTaskThreeJets, 1)
};

#endif
