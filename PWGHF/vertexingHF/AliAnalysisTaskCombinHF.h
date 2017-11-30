#ifndef ALIANALYSISTASKCOMBINHF_H
#define ALIANALYSISTASKCOMBINHF_H

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///*************************************************************************
/// \class Class AliAnalysisTaskCombinHF
/// \brief AliAnalysisTaskSE to build D meson candidates by combining tracks
///  background is computed LS and track rotations is
/// \author Authors: F. Prino, A. Rossi
//////////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TH3F.h>
#include <TObjString.h>
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
  virtual void FinishTaskOutput();
  
  void SetReadMC(Bool_t read){fReadMC=read;}

  void SetEventMixingWithCuts(Double_t maxDeltaVz, Double_t maxDeltaMult){
    fDoEventMixing=2; fMaxzVertDistForMix=maxDeltaVz; fMaxMultDiffForMix=maxDeltaMult;
  }
  void SetEventMixingWithPools(){fDoEventMixing=1;}
  void SetEventMixingOff(){fDoEventMixing=0;}
  void SetNumberOfEventsForMixing(Int_t minn){fNumberOfEventsForMixing=minn;}

  void ConfigureZVertPools(Int_t nPools, Double_t*  zVertLimits);
  void ConfigureMultiplicityPools(Int_t nPools, Double_t*  multLimits);
  void SelectPromptD(){fPromptFeeddown=kPrompt;}
  void SelectFeeddownD(){fPromptFeeddown=kFeeddown;}
  void SelectPromptAndFeeddownD(){fPromptFeeddown=kBoth;}
  void SetGoUpToQuark(Bool_t opt){fGoUpToQuark=opt;}
  void SetKeepNegIDtracks(Bool_t nid){fKeepNegID=nid;}//set it to kTRUE only if you know what you are doing
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
  void SetCutOnKKInvMass(Double_t cut){
    fPhiMassCut=cut;
  }
  void SetCutOnCos3PiKPhiRFrame(Double_t cut){
    fCutCos3PiKPhiRFrame=cut;
  }
  void SetCutCosPiDsLabFrame(Double_t cut){
    fCutCosPiDsLabFrame=cut;
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
  void ConfigureRotation3rdProng(Int_t n, Double_t phimin, Double_t phimax){
    fNRotations3=n;
    fMinAngleForRot3=phimin;
    fMaxAngleForRot3=phimax;
  }
  void SetMassWindow(Double_t minMass, Double_t maxMass){fMinMass=minMass; fMaxMass=maxMass;}
  void SetMaxPt(Double_t maxPt){fMaxPt=maxPt;}
  void SetPtBinWidth(Double_t binw){fPtBinWidth=binw;}
  void SetEtaAccCut(Double_t etacut){fEtaAccCut=etacut;}
  void SetPtAccCut(Double_t ptcut){fPtAccCut=ptcut;}
  void SetMultiplicityRange(Double_t mmin=-0.5, Double_t mmax=199.5){
    fMinMultiplicity=mmin;
    fMaxMultiplicity=mmax;
  }

  void SetPIDstrategy(Int_t strat){fPIDstrategy=strat;}
  void SetMaxPforIDPion(Double_t maxpIdPion){fmaxPforIDPion=maxpIdPion;}
  void SetMaxPforIDKaon(Double_t maxpIdKaon){fmaxPforIDKaon=maxpIdKaon;}
  void SetPIDselCaseZero(Int_t strat){fPIDselCaseZero=strat;}
  void SetBayesThres(Double_t thresKaon, Double_t thresPion){
    fBayesThresKaon=thresKaon;
    fBayesThresPion=thresPion;
  }
  
  Bool_t IsTrackSelected(AliAODTrack* track);
  Bool_t IsKaon(AliAODTrack* track);
  Bool_t IsPion(AliAODTrack* track);
  Bool_t SelectAODTrack(AliAODTrack *track, AliESDtrackCuts *cuts);
  
  Bool_t FillHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, TClonesArray *arrayMC, Int_t* dgLabels);
  void FillLSHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, Int_t charge);
  void FillMEHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau);
  void FillMEHistosLS(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, Int_t charge);
  void FillGenHistos(TClonesArray* arrayMC, Bool_t isEvSel);
  Bool_t CheckAcceptance(TClonesArray* arrayMC, Int_t nProng, Int_t *labDau); 
  Int_t GetPoolIndex(Double_t zvert, Double_t mult);
  void ResetPool(Int_t poolIndex);
  void DoMixingWithPools(Int_t poolIndex);
  void DoMixingWithCuts();
  Bool_t CanBeMixed(Double_t zv1, Double_t zv2, Double_t mult1, Double_t mult2);
  enum EMesonSpecies {kDzero, kDplus, kDstar, kDs};
  enum EPrompFd {kNone,kPrompt,kFeeddown,kBoth};
  enum EPIDstrategy {knSigma, kBayesianMaxProb, kBayesianThres};
  
private:
  
  AliAnalysisTaskCombinHF(const AliAnalysisTaskCombinHF &source);
  AliAnalysisTaskCombinHF& operator=(const AliAnalysisTaskCombinHF& source);
  Double_t ComputeInvMassKK(AliAODTrack* tr1, AliAODTrack* tr2) const;

  TList   *fOutput; //!<! list send on output slot 0
  TH1F *fHistNEvents;         //!<!hist. for No. of events
  TH2F *fHistEventMultZv;       //!<!hist. of evnt Mult vs. Zv for all events
  TH2F *fHistEventMultZvEvSel;  //!<!hist. of evnt Mult vs. Zv for selected ev
  TH1F *fHistTrackStatus;     //!<!hist. of status of tracks
  TH3F* fHistTrackEtaMultZv;  // track distribution vs. era z vertex and mult
  TH1F *fHistCheckOrigin;     //!<!hist. of origin (c/b) of D meson
  TH1F *fHistCheckOriginSel;  //!<!hist. of origin (c/b) of D meson
  TH1F *fHistCheckDecChan;    //!<!hist. of decay channel of D meson
  TH1F *fHistCheckDecChanAcc; //!<!hist. of decay channel of D meson in acc.
  TH3F *fPtVsYVsMultGen;        //!<! hist. of Y vs. Pt vs. Mult generated (all D)
  TH3F *fPtVsYVsMultGenLargeAcc; //!<! hist. of Y vs. Pt vs. Mult generated (|y|<0.9)
  TH3F *fPtVsYVsMultGenLimAcc;  //!<! hist. of Y vs. Pt vs. Mult generated (|y|<0.5)
  TH3F *fPtVsYVsMultGenAcc;     //!<! hist. of Y vs. Pt vs. Mult generated (D in acc)
  TH3F *fPtVsYVsMultGenAccEvSel; //!<! hist. of Y vs. Pt vs. Mult generated (D in acc, sel ev.)
  TH3F *fPtVsYVsMultReco;       //!<! hist. of Y vs. Pt vs. Mult generated (Reco D)
  TH3F *fMassVsPtVsY;     //!<! hist. of Y vs. Pt vs. Mass (all cand)
  TH3F *fMassVsPtVsYRot;   //!<! hist. of Y vs. Pt vs. Mass (rotations)
  TH3F *fMassVsPtVsYLSpp;  //!<! hist. of Y vs. Pt vs. Mass (like sign ++)
  TH3F *fMassVsPtVsYLSmm;  //!<! hist. of Y vs. Pt vs. Mass (like sign --)
  TH3F *fMassVsPtVsYSig;   //!<! hist. of Y vs. Pt vs. Mass (signal)
  TH3F *fMassVsPtVsYRefl;  //!<! hist. of Y vs. Pt vs. Mass (reflections)
  TH3F *fMassVsPtVsYBkg;   //!<! hist. of Y vs. Pt vs. Mass (background)
  TH1F *fNSelected;        //!<! hist. of n. of selected D+
  TH1F *fNormRotated;      //!<! hist. rotated/selected D+
  TH1F *fDeltaMass;        //!<! hist. mass difference after rotations
  THnSparse *fDeltaMassFullAnalysis; //!<! hist. mass difference after rotations with more details
  TH3F *fMassVsPtVsYME;   //!<! hist. of Y vs. Pt vs. Mass (mixedevents)
  TH3F *fMassVsPtVsYMELSpp;   //!<! hist. of Y vs. Pt vs. Mass (mixedevents)
  TH3F *fMassVsPtVsYMELSmm;   //!<! hist. of Y vs. Pt vs. Mass (mixedevents)
  TH2F* fEventsPerPool;   //!<! hist with number of events per pool  
  TH2F* fMixingsPerPool;    //!<! hist with number of mixings per pool  
  UInt_t fFilterMask; /// FilterMask
  AliESDtrackCuts* fTrackCutsAll; //// track selection
  AliESDtrackCuts* fTrackCutsPion; /// pion track selection
  AliESDtrackCuts* fTrackCutsKaon; /// kaon track selection
  Double_t fPhiMassCut;   /// cut on the KK inv mass for phi selection
  Double_t fCutCos3PiKPhiRFrame; // cut on the Ds decay angles
  Double_t fCutCosPiDsLabFrame;  // cut on the Ds decay angles
  AliAODPidHF* fPidHF; /// PID configuration
  AliRDHFCuts *fAnalysisCuts; /// Cuts for candidates
  
  Double_t fMinMass; /// minimum value of invariant mass
  Double_t fMaxMass; /// maximum value of invariant mass
  Double_t fMaxPt;   /// maximum pT value for inv. mass histograms
  Double_t fPtBinWidth; /// width of pt bin (GeV/c)
  Double_t fEtaAccCut; /// eta limits for acceptance step
  Double_t fPtAccCut; /// pt limits for acceptance step
  
  Int_t fNRotations; /// number of rotations
  Double_t fMinAngleForRot; /// minimum angle for track rotation
  Double_t fMaxAngleForRot; /// maximum angle for track rotation
  Int_t fNRotations3; /// number of rotations (3rd prong)
  Double_t fMinAngleForRot3; /// minimum angle for track rotation (3rd prong)
  Double_t fMaxAngleForRot3; /// maximum angle for track rotation (3rd prong)
  
  AliNormalizationCounter *fCounter;//!<!Counter for normalization
  
  Int_t fMeson;          /// mesonSpecies (see enum)
  Bool_t  fReadMC;       ///  flag for access to MC
  Int_t fPromptFeeddown; /// flag to select prompt (1), feeddown (2) or all (3)
  Bool_t fGoUpToQuark;   /// flag for definition of c,b origin
  Int_t fFullAnalysis;   /// flag to set analysis level (0 is the fastest)
  
  Int_t    fPIDstrategy;   /// knSigma, kBayesianMaxProb, kBayesianThres
  Double_t fmaxPforIDPion; /// flag for upper p limit for id band for pion
  Double_t fmaxPforIDKaon; /// flag for upper p limit for id band for kaon
  Bool_t   fKeepNegID;    /// flag to keep also track with negative ID (default kFALSE, change it only if you know what you are doing)
  Int_t    fPIDselCaseZero;  /// flag to change PID strategy
  Double_t fBayesThresKaon;  /// threshold for kaon identification via Bayesian PID
  Double_t fBayesThresPion;  /// threshold for pion identification via Bayesian PID

  Int_t fDoEventMixing; /// flag for event mixing
  Int_t  fNumberOfEventsForMixing; /// maximum number of events to be used in event mixing
  Double_t fMaxzVertDistForMix; /// cut on zvertex distance for event mixing with cuts
  Double_t fMaxMultDiffForMix; /// cut on multiplicity difference for event mixing with cuts
  Int_t fNzVertPools; /// number of pools in z vertex for event mixing
  Int_t fNzVertPoolsLimSize; /// number of pools in z vertex for event mixing +1
  Double_t* fzVertPoolLims; //[fNzVertPoolsLimSize] limits of the pools in zVertex
  Int_t fNMultPools; /// number of pools in multiplicity for event mixing
  Int_t fNMultPoolsLimSize; /// number of pools in multiplicity for event mixing +1
  Double_t* fMultPoolLims; //[fNMultPoolsLimSize] limits of the pools in multiplicity
  Int_t  fNOfPools; /// number of pools
  TTree** fEventBuffer;   //!<! structure for event mixing
  TObjString* fEventInfo;  /// unique event Id for event mixing checks
  Double_t fVtxZ;         /// zVertex
  Double_t fMultiplicity; /// multiplicity
  Double_t fMinMultiplicity;  /// lower limit for multiplcities in MC histos
  Double_t fMaxMultiplicity; /// upper limit for multiplcities in MC histos
  TObjArray* fKaonTracks; /// array of kaon-compatible tracks (TLorentzVectors)
  TObjArray* fPionTracks; /// array of pion-compatible tracks (TLorentzVectors)

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCombinHF,13); /// D0D+ task from AOD tracks
  /// \endcond
};

#endif
