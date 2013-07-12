#ifndef ALIANALYSISTASKPHIBAYES_H
#define ALIANALYSISTASKPHIBayes_H

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

class TH1F;
class TH2F;
class AliESDtrackCuts;

class AliAnalysisTaskPhiBayes : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPhiBayes();
  AliAnalysisTaskPhiBayes(const char *name);

  virtual ~AliAnalysisTaskPhiBayes();

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

  static const Int_t nPtBin = 13;  //! # pt phi bins

  void SetFilterBit(Int_t fb){fFilterBit=fb;};
  Int_t GetFilterBit() const {return fFilterBit;};

 private:
  AliAnalysisTaskPhiBayes(const AliAnalysisTaskPhiBayes &old); 
  AliAnalysisTaskPhiBayes& operator=(const AliAnalysisTaskPhiBayes &source); 

  virtual Float_t GetVertex(AliAODEvent* aod) const;
  virtual void Analyze(AliAODEvent* aodEvent); 
  virtual Int_t IsChannelValid(Float_t etaAbs);

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

  Float_t fPtPhi;   //! variable to fill the tree
  Float_t fPhiPhi;  //! variable to fill the tree
  Float_t fEtaPhi;  //! variable to fill the tree
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

  static Float_t fPtPhiMin[nPtBin];// ptmin bin
  static Float_t fPtPhiMax[nPtBin];// ptmax bin

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

  TH1F *fHmismTOF;         //! TOF mismatch distribution
  TH1D *fHchannelTOFdistr; //! TOF channel distance w.r.t. IP

  ClassDef(AliAnalysisTaskPhiBayes, 1);    //Analysis task for bayesian (K0s)
};

#endif
