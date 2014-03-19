#ifndef ALIANALYSISTASKCOMBINHF_H
#define ALIANALYSISTASKCOMBINHF_H

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */ 

//*************************************************************************
// Class AliAnalysisTaskCombinHF
// AliAnalysisTaskSE to build D meson candidates by combining tracks
//  background is computed LS and track rotations is
// Authors: F. Prino, A. Rossi
/////////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliNormalizationCounter.h"
#include "AliRDHFCuts.h"

class AliAnalysisTaskCombinHF : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskCombinHF();
  AliAnalysisTaskCombinHF(Int_t meson, AliRDHFCuts* analysiscuts);
  virtual ~AliAnalysisTaskCombinHF();

  virtual void UserCreateOutputObjects();
  virtual void Init(){};
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetReadMC(Bool_t read){fReadMC=read;}
  void SelectPromptD(){fPromptFeeddown=kPrompt;}
  void SelectFeeddownD(){fPromptFeeddown=kFeeddown;}
  void SelectPromptAndFeeddownD(){fPromptFeeddown=kBoth;}
  void SetGoUpToQuark(Bool_t opt){fGoUpToQuark=opt;}

  void SetTrackCuts(AliESDtrackCuts* cuts){
    if(fTrackCutsAll) delete fTrackCutsAll;
    fTrackCutsAll=new AliESDtrackCuts(*cuts);
  }
  void SetPionTrackCuts(AliESDtrackCuts* cuts){
    if(fTrackCutsPion) delete fTrackCutsPion;
    fTrackCutsPion=new AliESDtrackCuts(*cuts);
  }
  void SetKaonTrackCuts(AliESDtrackCuts* cuts){
    if(fTrackCutsKaon) delete fTrackCutsKaon;
    fTrackCutsKaon=new AliESDtrackCuts(*cuts);
  }
  void SetPIDHF(AliAODPidHF* pid){
    if(fPidHF) delete fPidHF;
    fPidHF=new AliAODPidHF(*pid);
  }
  void SetRDHFCuts(AliRDHFCuts* cuts){
    fAnalysisCuts=cuts;
  }
  void SetFilterMask(UInt_t mask=16){fFilterMask=mask;} 
  void SetAnalysisLevel(Int_t level){fFullAnalysis=level;}
  void ConfigureRotation(Int_t n, Double_t phimin, Double_t phimax){
    fNRotations=n;
    fMinAngleForRot=phimin;
    fMaxAngleForRot=phimax;
  }
  Bool_t IsTrackSelected(AliAODTrack* track);
  Bool_t IsKaon(AliAODTrack* track);
  Bool_t IsPion(AliAODTrack* track);
  Bool_t SelectAODTrack(AliAODTrack *track, AliESDtrackCuts *cuts);
  Bool_t FillHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, TClonesArray *arrayMC, Int_t* dgLabels);
  void FillLSHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, Int_t charge);
  void FillGenHistos(TClonesArray* arrayMC);
  Bool_t CheckAcceptance(TClonesArray* arrayMC, Int_t nProng, Int_t *labDau);
  enum EMesonSpecies {kDzero, kDplus, kDstar, kDs};
  enum EPrompFd {kNone,kPrompt,kFeeddown,kBoth};

 private:

  AliAnalysisTaskCombinHF(const AliAnalysisTaskCombinHF &source);
  AliAnalysisTaskCombinHF& operator=(const AliAnalysisTaskCombinHF& source); 

  TList   *fOutput; //! list send on output slot 0
  TH1F *fHistNEvents;         //!hist. for No. of events
  TH1F *fHistTrackStatus;     //!hist. of status of tracks
  TH1F *fHistCheckOrigin;     //!hist. of origin (c/b) of D meson
  TH1F *fHistCheckOriginSel;  //!hist. of origin (c/b) of D meson
  TH1F *fHistCheckDecChan;    //!hist. of decay channel of D meson
  TH1F *fHistCheckDecChanAcc; //!hist. of decay channel of D meson in acc.
  TH2F *fPtVsYGen;        //! hist. of Y vs. Pt generated (all D)
  TH2F *fPtVsYGenLimAcc;  //! hist. of Y vs. Pt generated (|y|<0.5)
  TH2F *fPtVsYGenAcc;     //! hist. of Y vs. Pt generated (D in acc)
  TH2F *fPtVsYReco;       //! hist. of Y vs. Pt generated (Reco D)
  TH3F *fMassVsPtVsY;     //! hist. of Y vs. Pt vs. Mass (all cand)
  TH3F *fMassVsPtVsYRot;   //! hist. of Y vs. Pt vs. Mass (rotations)
  TH3F *fMassVsPtVsYLSpp;  //! hist. of Y vs. Pt vs. Mass (like sign ++)
  TH3F *fMassVsPtVsYLSmm;  //! hist. of Y vs. Pt vs. Mass (like sign --)
  TH3F *fMassVsPtVsYSig;   //! hist. of Y vs. Pt vs. Mass (signal)
  TH3F *fMassVsPtVsYRefl;  //! hist. of Y vs. Pt vs. Mass (reflections)
  TH3F *fMassVsPtVsYBkg;   //! hist. of Y vs. Pt vs. Mass (background)
  TH1F *fNSelected;        //! hist. of n. of selected D+
  TH1F *fNormRotated;      //! hist. rotated/selected D+
  TH1F *fDeltaMass;        //! hist. mass difference after rotations
  THnSparse *fDeltaMassFullAnalysis; //! hist. mass difference after rotations with more details

  UInt_t fFilterMask; // FilterMask
  AliESDtrackCuts* fTrackCutsAll; // track selection
  AliESDtrackCuts* fTrackCutsPion; // pion track selection
  AliESDtrackCuts* fTrackCutsKaon; // kaon track selection
  AliAODPidHF* fPidHF; // PID configuration
  AliRDHFCuts *fAnalysisCuts; // Cuts for candidates
  Double_t fMinMass; // minimum value of invariant mass
  Double_t fMaxMass; // maximum value of invariant mass

  Double_t fEtaAccCut; // eta limits for acceptance step
  Double_t fPtAccCut; // pt limits for acceptance step


  Int_t fNRotations; // number of rotations
  Double_t fMinAngleForRot; // minimum angle for track rotation
  Double_t fMaxAngleForRot; // maximum angle for track rotation
  Double_t fMinAngleForRot3; // minimum angle for track rotation (3rd prong)
  Double_t fMaxAngleForRot3; // maximum angle for track rotation (3rd prong)

  AliNormalizationCounter *fCounter;//!Counter for normalization

  Int_t fMeson;          // mesonSpecies (see enum)
  Bool_t  fReadMC;       //  flag for access to MC
  Int_t fPromptFeeddown; // flag to select prompt (1), feeddown (2) or all (3) 
  Bool_t fGoUpToQuark;   // flag for definition of c,b origin
  Int_t fFullAnalysis;   // flag to set analysis level (0 is the fastest)

  ClassDef(AliAnalysisTaskCombinHF,2); // D+ task from AOD tracks
};

#endif
