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

/*   virtual void SetUseBkg(Bool_t b) {fUseBkg = b;} */
/*   virtual void SetUseSimTPC(Bool_t b) {fUseFastTPC = b;} */

  virtual void FillTopology(TH2F * Dalitz, TH1F * fhMu34, TH1F * fhMu45, TH1F * fhMu35, Double_t * x, TVector3 * pRest, Double_t xsection);

  virtual Double_t SetR(Double_t b){fR = b; return fR;} 
//  TArrayD * GetThrustParamMC(AliMCEvent* mcEvent, Int_t  nstudymin, Double_t ptcutoff, Double_t etacutoff, Bool_t chom, TArrayD * evsh);
  virtual void GetThrustAxis(TVector3 &n01, TVector3 * p_track,const Int_t &nTracks);
  virtual void GetEventShapes(TVector3 &n01, TVector3 * pTrack, Int_t nTracks, Double_t * eventShapes);
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

  TH2F * fhX3X4Rec; //
  TH2F * fhX3X4Gen; //
  
  TH1F * fhMu35Rec; //
  TH1F * fhMu34Rec; //
  TH1F * fhMu45Rec; //
  
  TH1F * fhMu35Gen; //
  TH1F * fhMu34Gen; //
  TH1F * fhMu45Gen; //

  TH1I * fhInOut; //
  TH1F * fhThrustRec2; //
  TH1F * fhThrustRec3; // 

  TH1F * fhThrustGen2; //
  TH1F * fhThrustGen3; // 

  TH1F * fhCGen2; //
  TH1F * fhCGen3; // 

  TH1F * fhSGen2; //
  TH1F * fhSGen3; // 

  TH1F * fhAGen2; // 
  TH1F * fhAGen3; //

  TH1F * fhCRec2; // 
  TH1F * fhCRec3; // 

  TH1F * fhSRec2; //
  TH1F * fhSRec3; //

  TH1F * fhARec2; //
  TH1F * fhARec3; //

  TH2F * fhX3; //
  TH2F * fhX4; //
  TH2F * fhX5; //

  TProfile * fhXSec; //
  TH2F * fhX3X4Rec60; //
  TH2F * fhX3X4Rec60100; //
  TH2F * fhX3X4Rec100; //
  TH2F * fhX3X4Gen60; //
  TH2F * fhX3X4Gen60100; //
  TH2F * fhX3X4Gen100; //

  TH2F * fhdPhiThrustGen; //
  TH2F * fhdPhiThrustGenALL; //
  TH2F * fhdPhiThrustRec; //
  TH2F * fhdPhiThrustRecALL; //

  ClassDef(AliAnalysisTaskThreeJets, 1)
};

#endif
