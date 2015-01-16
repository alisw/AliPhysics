#ifndef ALIANALYSISTASKLAMBDABAYES_H
#define ALIANALYSISTASKLAMBDABAYES_H

// ROOT includes
#include <TObject.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TProfile.h>
#include "AliPIDCombined.h"
#include "TF1.h"

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliAODEvent.h>
#include "AliPIDperfContainer.h"
#include "AliPIDperfCut.h"

class TH1F;
class TH2F;
class AliESDtrackCuts;

class AliAnalysisTaskLambdaBayes : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskLambdaBayes();
  AliAnalysisTaskLambdaBayes(const char *name);

  virtual ~AliAnalysisTaskLambdaBayes();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *); 
  
  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }     
  Double_t GetMinPt() { return fMinPt; }   

  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetMinPt(Double_t value) {fMinPt = value;}   

  virtual void SetMC(Bool_t flag = kTRUE){fIsMC = flag;};
  virtual void SetQA(Bool_t flag = kTRUE){fQAsw = flag;};

  void SetTPCclusterN(Int_t ncl){fNcluster=ncl;};

  void SetCorrEP(Bool_t flag=kTRUE) {fToEP = flag;};

  static const Int_t nPtBin = 13;  //! # pt ks bins

  void SetFilterBit(Int_t fb){fFilterBit=fb;};
  Int_t GetFilterBit() const {return fFilterBit;};

  void SetTypeCollisions(Int_t type){fTypeCol = type;}; // 0=pp, 1=pPb, 2=PbPb
  Int_t GetTypeCollisions() const {return fTypeCol;}; 

  void SetPIDuserCut(AliPIDperfCut *usercut){fPIDuserCut = usercut;};

  void SetRefSpecies(Int_t species) {fSpeciesRef = species;}

 private:
  AliAnalysisTaskLambdaBayes(const AliAnalysisTaskLambdaBayes &old); 
  AliAnalysisTaskLambdaBayes& operator=(const AliAnalysisTaskLambdaBayes &source); 

  Int_t PassesAODCuts(AliAODv0 *myV0, AliAODEvent *tAOD,Int_t specie);
  Int_t FindDaugheterIndex(AliAODTrack *trk);

  virtual Float_t GetVertex(AliAODEvent* aod) const;
  virtual void Analyze(AliAODEvent* aodEvent); 
  virtual Int_t IsChannelValid(Float_t etaAbs);
  virtual void SelectLambda();

  Double_t     fVtxCut;             // Vtx cut on z position in cm
  Double_t     fEtaCut;             // Eta cut used to select particles
  Double_t     fMinPt;              // Min pt - for histogram limits

  Bool_t fIsMC; // if MC
  Bool_t fQAsw;   // if QA

  static const Int_t nCentrBin = 9; //! # cenrality bins

  //
  // Cuts and options
  //

  Int_t fNcluster;           // Numer of TPC cluster required
  Int_t fFilterBit;          // filter bit to be used
  TList *fList;              //! List for output objects
  TList *fList2;             //! List for output objects
  TList *fList3;             //! List for output objects

  //
  Float_t fCentrality;  //! current centrality for the tree
  Float_t fPsi;         //! current Psi from TPC

  Float_t fPtLambdaC;   //! variable to fill the tree
  Float_t fPhiLambdaC;  //! variable to fill the tree
  Float_t fEtaLambda;  //! variable to fill the tree
  Float_t fMassV0; //! variable to fill the tree
  Float_t fPtKp;  //! variable to fill the tree
  Float_t fPhiKp; //! variable to fill the tree
  Float_t fEtaKp; //! variable to fill the tree
  Float_t fPtKn;  //! variable to fill the tree
  Float_t fPhiKn; //! variable to fill the tree
  Float_t fEtaKn; //! variable to fill the tree
  Int_t fPidKp;   //! variable to fill the tree
  Int_t fPidKn;   //! variable to fill the tree

  TH2F *hMatching[4]; //! matching (all, pi, k ,p)
  TH2F *hTracking[4]; //! tracking (all, pi, k ,p)
 
  TH2F *fTOFTPCsignal;   //! histo with tof signal
  TH1F *fCombsignal;   //! histo with tof signal

  static Float_t fPtLambdaMin[nPtBin];// ptmin bin
  static Float_t fPtLambdaMax[nPtBin];// ptmax bin

  // QA plots
  TH2F *fPiTPC[nCentrBin];//! TPC dE/dx plot
  TH2F *fKaTPC[nCentrBin];//! TPC dE/dx plot
  TH2F *fPrTPC[nCentrBin];//! TPC dE/dx plot
  TH2F *fElTPC[nCentrBin];//! TPC dE/dx plot

  TH2F *fPiTOF[nCentrBin];//! TPC dE/dx plot
  TH2F *fKaTOF[nCentrBin];//! TPC dE/dx plot
  TH2F *fPrTOF[nCentrBin];//! TPC dE/dx plot
  TH2F *fElTOF[nCentrBin];//! TPC dE/dx plot

  AliESDtrackCuts *fCutsDaughter;

  AliPIDCombined *fPIDCombined;  //! PID combined object
  AliPIDperfContainer* fContPid;      //! results for positive
  AliPIDperfContainer* fContPid2;     //! results for negative
  AliPIDperfContainer* fContUser;     //! results for positive user cut
  AliPIDperfContainer* fContUser2;    //! results for negative user cut

  Int_t fNLambda; //! number of Lambda in my private selection
  Float_t fPhiLambda[1000]; //! phi of Lambda in my private selection
  Float_t fPtLambda[1000];//! pt of Lambda in my private selection
  Int_t fNpPos; //! number of positive pions for Lambda selection
  Int_t fNpNeg; //! number of negative pions for Lambda selection
  Int_t fIPPos[1000]; //! position in the AOD stack of positive pions candidates
  Int_t fIPNeg[1000]; //! position in the AOD stack of negative pions candidates
  Int_t fIpP[1000]; //! position in the AOD stack of positive pions for Lambda
  Int_t fIpN[1000]; //! position in the AOD stack of negative pions for Lambda
  Float_t fMassLambda[1000]; //! Lambda mass
  
  TH1F *fHmismTOF; //! TOF mismatch distribution
  TH1D *fHchannelTOFdistr; //! TOF channel distance w.r.t. IP

  Int_t fTypeCol; // type of collision system (0=pp, 1=pPb, 2=PbPb)
  AliPIDperfCut *fPIDuserCut;  // pid user cut to be cheked

  Bool_t fToEP;        // correlate to TPC EP
  Int_t fSpeciesRef;    // reference species (default = 4, proton)

  ClassDef(AliAnalysisTaskLambdaBayes, 5);    //Analysis task for Bayesian (Lambda)
};

#endif
