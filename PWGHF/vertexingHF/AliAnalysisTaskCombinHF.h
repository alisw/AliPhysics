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
  void UseOnlySignalInMC(Bool_t opt){fSignalOnlyMC=opt;}
  void UseMBTrigMaskInMC(){fEnforceMBTrigMaskInMC=kTRUE;}
  void UseTrigMaskFromCutFileInMC(){fEnforceMBTrigMaskInMC=kFALSE;}
  void SetPtHardRange(double pmin, double pmax){
    fSelectPtHardRange=kTRUE; fMinPtHard=pmin; fMaxPtHard=pmax;
  }
  void SetRejectGeneratedEventsWithPileup(Bool_t opt=kTRUE){
    fRejectGeneratedEventsWithPileup=opt;
  }
  void SetRejectSignalsFromOOBPileupEvents(Bool_t opt=kTRUE){
    fRejectSignalsFromOOBPileupEvents=opt;
  }
  void SetTrackletEta1MultiplicityEstimatorForMC(){fMultEstimMC=0;}
  void SetTrackletFullEtaMultiplicityEstimatorForMC(){fMultEstimMC=1;}
  void SetCentralityPercMultiplicityEstimatorForMC(){fMultEstimMC=2;}
  void SetTPCClustersMultiplicityEstimatorForMC(){fMultEstimMC=3;}
  
  void SetEventMixingWithCuts(Double_t maxDeltaVz, Double_t maxDeltaMult){
    fDoEventMixing=2; fMaxzVertDistForMix=maxDeltaVz; fMaxMultDiffForMix=maxDeltaMult;
  }
  void SetEventMixingWithPools(){fDoEventMixing=1;}
  void SetEventMixingOff(){fDoEventMixing=0;}
  void SetNumberOfEventsForMixing(Int_t minn){fNumberOfEventsForMixing=minn;}

  void ConfigureZVertPools(Int_t nPools, Double_t*  zVertLimits);
  void ConfigureMultiplicityPools(Int_t nPools, Double_t*  multLimits);
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
  void SetMinNumTPCClsForPID(Int_t cut=0.) {fCutTPCSignalN = cut;}
  void SetCutOnCosThetaStar(Double_t cut){
    if(cut>0 && cut<1) fApplyCutCosThetaStar=kTRUE;
    else fApplyCutCosThetaStar=kFALSE;
    fCutCosThetaStar=cut;
  }
  void EnableHistosVsCosThetaStar(Bool_t opt){
    fFillHistosVsCosThetaStar=opt;
  }
  void SetUseDzeroTopologicalCuts(Bool_t opt){
    fUseDzeroTopologicalCuts=opt;
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
  void SetMultiplicityRange(Double_t mmin=-0.5, Double_t mmax=199.5, Int_t nbins=200){
    fMinMultiplicity=mmin;
    fMaxMultiplicity=mmax;
    fNumOfMultBins=nbins;
  }

  void SetPIDstrategy(Int_t strat){fPIDstrategy=strat;}
  void SetMaxPforIDPion(Double_t maxpIdPion){fmaxPforIDPion=maxpIdPion;}
  void SetMaxPforIDKaon(Double_t maxpIdKaon){fmaxPforIDKaon=maxpIdKaon;}
  void SetPIDselCaseZero(Int_t strat){fPIDselCaseZero=strat;}
  void SetBayesThres(Double_t thresKaon, Double_t thresPion, Double_t thresProton=0.4){
    fBayesThresKaon=thresKaon;
    fBayesThresPion=thresPion;
    fBayesThresProton=thresProton;
  }
  
  Bool_t IsTrackSelected(AliAODTrack* track);
  Bool_t IsKaon(AliAODTrack* track);
  Bool_t IsPion(AliAODTrack* track);
  Bool_t IsProton(AliAODTrack* track);
  Bool_t SelectAODTrack(AliAODTrack *track, AliESDtrackCuts *cuts);
  
  Bool_t FillHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, TClonesArray *arrayMC, AliAODMCHeader *mcHeader, Int_t* dgLabels);
  void FillLSHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, Int_t charge);
  void FillMEHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau);
  void FillMEHistosLS(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, Int_t charge);
  void FillGenHistos(TClonesArray* arrayMC, AliAODMCHeader *mcHeader, Bool_t isEvSel);
  Bool_t CheckAcceptance(TClonesArray* arrayMC, Int_t nProng, Int_t *labDau); 
  Int_t GetPoolIndex(Double_t zvert, Double_t mult);
  void ResetPool(Int_t poolIndex);
  void DoMixingWithPools(Int_t poolIndex);
  void DoMixingWithCuts();
  Bool_t CanBeMixed(Double_t zv1, Double_t zv2, Double_t mult1, Double_t mult2);
  enum EMesonSpecies {kDzero, kDplus, kDstar, kDs, kJpsi, kEtac};
  enum EPrompFd {kNone,kPrompt,kFeeddown,kBoth};
  enum EPIDstrategy {knSigma, kBayesianMaxProb, kBayesianThres};
  
private:
  
  AliAnalysisTaskCombinHF(const AliAnalysisTaskCombinHF &source);
  AliAnalysisTaskCombinHF& operator=(const AliAnalysisTaskCombinHF& source);
  Double_t ComputeInvMassKK(AliAODTrack* tr1, AliAODTrack* tr2) const;
  Double_t ComputeInvMassKK(TLorentzVector* tr1, TLorentzVector* tr2) const;
  Double_t CosPiKPhiRFrame(TLorentzVector* dauK1, TLorentzVector* dauK2, TLorentzVector* daupi) const;
  Double_t CosPiDsLabFrame(TLorentzVector* dauK1, TLorentzVector* dauK2, TLorentzVector* daupi) const;

  TList *fOutput;                       //!<! list with output histograms
  TList *fListCuts;                     //!<! list with cut values 
  TH1F *fHistNEvents;                   //!<! hist. for N of events
  TH1F *fHistNEventsMCCharmInj;         //!<! hist. for N of events with ccbar
  TH1F *fHistNEventsMCBeautyInj;        //!<! hist. for N of events with bbbar
  TH2F *fHistEventMultCent;             //!<! hist. for evnt Mult vs. centrality (all)
  TH2F *fHistEventMultCentEvSel;        //!<! hist. for evnt Mult vs. centrality (sel)
  TH2F *fHistEventMultZv;               //!<! hist. of evnt Mult vs. Zv for all events
  TH2F *fHistEventMultZvEvSel;          //!<! hist. of evnt Mult vs. Zv for selected ev
  TH2F *fHistEventTrackletCent;             //!<! hist. for evnt Tracklet vs. centrality (all)
  TH2F *fHistEventTrackletCentEvSel;        //!<! hist. for evnt Tracklet vs. centrality (sel)
  TH2F *fHistEventTrackletZv;               //!<! hist. of evnt Tracklet vs. Zv for all events
  TH2F *fHistEventTrackletZvEvSel;          //!<! hist. of evnt Tracklet vs. Zv for selected ev
  TH1F *fHistXsecVsPtHard;              //!<! hist. of xsec vs pthard (MC)
  TH1F *fHistTrackStatus;               //!<! hist. of status of tracks
  TH3F *fHistTrackEtaMultZv;            //!<! track distribution vs. eta z vertex and mult
  TH3F *fHistTrackEtaTrackletZv;            //!<! track distribution vs. eta z vertex and mult
  TH1D* fHistTrackSelSteps;             //!<! track cuts statistics
  TH2F *fHistSelTrackPhiPt;             //!<! track distribution vs. phi and pt
  TH2F *fHistSelTrackChi2ClusPt;        //!<! track chi2 distribution vs. pt
  TH2F *fHistSelTrackDCAxyPt;           //!<! impact patamter histos
  TH2F *fHistSelTrackFineDCAxyPt;           //!<! impact patamter histos
  TH2F *fHistSelTrackDCAzPt;            //!<! impact patamter histos
  TH2F *fHistSelTrackDCAxyPtAfterProp;  //!<! impact patamter histos
  TH2F *fHistSelTrackFineDCAxyPtAfterProp;  //!<! impact patamter histos
  TH2F *fHistSelTrackDCAzPtAfterProp;   //!<! impact patamter histos
  TH2F *fHistCheckOrigin;               //!<! hist. of origin (c/b) of D meson (gen)
  TH2F *fHistCheckOriginRecoD;          //!<! hist. of origin (c/b) of D meson (reco)
  TH2F *fHistCheckOriginRecoVsGen;      //!<! hist. of origin (c/b) of D meson
  TH1F *fHistCheckDecChan;              //!<! hist. of decay channel of D meson
  TH1F *fHistCheckDecChanAcc;           //!<! hist. of decay channel of D meson in acc.
  TH3F *fPtVsYVsMultGenPrompt;          //!<! hist. of Y vs. Pt vs. Mult generated (all D)
  TH3F *fPtVsYVsMultGenLargeAccPrompt;  //!<! hist. of Y vs. Pt vs. Mult generated (|y|<0.9)
  TH3F *fPtVsYVsMultGenLimAccPrompt;    //!<! hist. of Y vs. Pt vs. Mult generated (|y|<0.5)
  TH3F *fPtVsYVsMultGenAccPrompt;       //!<! hist. of Y vs. Pt vs. Mult generated (D in acc)
  TH3F *fPtVsYVsMultGenAccEvSelPrompt;  //!<! hist. of Y vs. Pt vs. Mult generated (D in acc, sel ev.)
  TH3F *fPtVsYVsMultRecoPrompt;         //!<! hist. of Y vs. Pt vs. Mult generated (Reco D)
  TH3F *fPtVsPhiVsMultGenPrompt;        //!<! hist. of Phi vs. Pt vs. Mult generated (all D)
  TH3F *fPtVsPhiVsMultGenLimAccPrompt;  //!<! hist. of Phi vs. Pt vs. Mult generated (all D)
  TH3F *fPtVsPhiVsMultGenAccPrompt;     //!<! hist. of Phi vs. Pt vs. Mult generated (all D)
  TH3F *fPtVsPhiVsMultRecoPrompt;       //!<! hist. of Y vs. Pt vs. Mult generated (Reco D)
  TH3F *fPtVsYVsMultGenFeeddw;          //!<! hist. of Y vs. Pt vs. Mult generated (all D)
  TH3F *fPtVsYVsMultGenLargeAccFeeddw;  //!<! hist. of Y vs. Pt vs. Mult generated (|y|<0.9)
  TH3F *fPtVsYVsMultGenLimAccFeeddw;    //!<! hist. of Y vs. Pt vs. Mult generated (|y|<0.5)
  TH3F *fPtVsYVsMultGenAccFeeddw;       //!<! hist. of Y vs. Pt vs. Mult generated (D in acc)
  TH3F *fPtVsYVsMultGenAccEvSelFeeddw;  //!<! hist. of Y vs. Pt vs. Mult generated (D in acc, sel ev.)
  TH3F *fPtVsYVsMultRecoFeeddw;         //!<! hist. of Y vs. Pt vs. Mult generated (Reco D)
  TH3F *fPtVsPhiVsMultGenFeeddw;        //!<! hist. of Phi vs. Pt vs. Mult generated (all D)
  TH3F *fPtVsPhiVsMultGenLimAccFeeddw;  //!<! hist. of Phi vs. Pt vs. Mult generated (all D)
  TH3F *fPtVsPhiVsMultGenAccFeeddw;     //!<! hist. of Phi vs. Pt vs. Mult generated (all D)
  TH3F *fPtVsPhiVsMultRecoFeeddw;       //!<! hist. of Phi vs. Pt vs. Mult generated (all D)
  TH3F *fPtVsYVsPtBGenFeeddw;           //!<! hist. of Y vs. Pt vs. PtB generated (all D)
  TH3F *fPtVsYVsPtBGenLargeAccFeeddw;   //!<! hist. of Y vs. Pt vs. PtB generated (|y|<0.9)
  TH3F *fPtVsYVsPtBGenLimAccFeeddw;     //!<! hist. of Y vs. Pt vs. PtB generated (|y|<0.5)
  TH3F *fPtVsYVsPtBGenAccFeeddw;        //!<! hist. of Y vs. Pt vs. PtB generated (D in acc)
  TH3F *fPtVsYVsPtBGenAccEvSelFeeddw;   //!<! hist. of Y vs. Pt vs. PtB generated (D in acc, sel ev.)
  TH3F *fPtVsYVsPtBRecoFeeddw;          //!<! hist. of Y vs. Pt vs. PtB generated (Reco D)
  TH3F *fMassVsPtVsY;      //!<! hist. of Y vs. Pt vs. Mass (all cand)
  TH3F *fMassVsPtVsYRot;   //!<! hist. of Y vs. Pt vs. Mass (rotations)
  TH3F *fMassVsPtVsYLSpp;  //!<! hist. of Y vs. Pt vs. Mass (like sign ++)
  TH3F *fMassVsPtVsYLSmm;  //!<! hist. of Y vs. Pt vs. Mass (like sign --)
  TH3F *fMassVsPtVsYSig;   //!<! hist. of Y vs. Pt vs. Mass (signal)
  TH3F *fMassVsPtVsYRefl;  //!<! hist. of Y vs. Pt vs. Mass (reflections)
  TH3F *fMassVsPtVsYBkg;   //!<! hist. of Y vs. Pt vs. Mass (background)
  TH1F *fBMohterPtGen;     //!<! hist. of beauty mother pt
  TH1F *fNSelected;        //!<! hist. of n. of selected D+
  TH1F *fNormRotated;      //!<! hist. rotated/selected D+
  TH1F *fDeltaMass;        //!<! hist. mass difference after rotations
  THnSparse *fDeltaMassFullAnalysis; //!<! hist. mass difference after rotations with more details
  TH3F *fMassVsPtVsYME;   //!<! hist. of Y vs. Pt vs. Mass (mixedevents)
  TH3F *fMassVsPtVsYMELSpp;   //!<! hist. of Y vs. Pt vs. Mass (mixedevents)
  TH3F *fMassVsPtVsYMELSmm;   //!<! hist. of Y vs. Pt vs. Mass (mixedevents)
  TH2F* fEventsPerPool;   //!<! hist with number of events per pool  
  TH2F* fMixingsPerPool;    //!<! hist with number of mixings per pool
  TH3F *fMassVsPtVsCosthSt;         //!<! hist. of Pt vs. Mass vs. cos(th*) (all cand)
  TH3F *fMassVsPtVsCosthStRot;      //!<! hist. of Pt vs. Mass vs. cos(th*) (rotations)
  TH3F *fMassVsPtVsCosthStLSpp;     //!<! hist. of Pt vs. Mass vs. cos(th*) (like sign ++)
  TH3F *fMassVsPtVsCosthStLSmm;     //!<! hist. of Pt vs. Mass vs. cos(th*) (like sign --)
  TH3F *fMassVsPtVsCosthStSig;      //!<! hist. of Pt vs. Mass vs. cos(th*) (signal)
  TH3F *fMassVsPtVsCosthStRefl;     //!<! hist. of Pt vs. Mass vs. cos(th*) (reflections)
  TH3F *fMassVsPtVsCosthStBkg;      //!<! hist. of Pt vs. Mass vs. cos(th*) (background)
  TH3F *fMassVsPtVsCosthStME;       //!<! hist. of Pt vs. Mass vs. cos(th*) (mixedevents)
  TH3F *fMassVsPtVsCosthStMELSpp;   //!<! hist. of Pt vs. Mass vs. cos(th*) (mixedevents)
  TH3F *fMassVsPtVsCosthStMELSmm;   //!<! hist. of Pt vs. Mass vs. cos(th*) (mixedevents)
  TH2F *fHistonSigmaTPCPion;        //!<! hist. of nSigmaTPC pion
  TH2F *fHistonSigmaTPCPionGoodTOF; //!<! hist. of nSigmaTPC pion
  TH2F *fHistonSigmaTOFPion;        //!<! hist. of nSigmaTOF pion
  TH2F *fHistonSigmaTPCKaon;        //!<! hist. of nSigmaTPC kaon
  TH2F *fHistonSigmaTPCKaonGoodTOF; //!<! hist. of nSigmaTPC kaon
  TH2F *fHistonSigmaTOFKaon;        //!<! hist. of nSigmaTOF kaon
  TH2F *fHistonSigmaTPCProton;        //!<! hist. of nSigmaTPC proton
  TH2F *fHistonSigmaTPCProtonGoodTOF; //!<! hist. of nSigmaTPC proton
  TH2F *fHistonSigmaTOFProton;        //!<! hist. of nSigmaTOF proton
  TH3F *fHistoPtKPtPiPtD;           //!<! hist. for propagation of tracking unc
  TH3F *fHistoPtKPtPiPtDSig;        //!<! hist. for propagation of tracking unc
  TH1F *fHistd0xd0;                 //!<! hist. to check topological cuts
  TH1F *fHistCosPoint;              //!<! hist. to check topological cuts
  TH1F *fHistCosPointXY;            //!<! hist. to check topological cuts
  TH1F *fHistDecLen;                //!<! hist. to check topological cuts
  TH1F *fHistNormDecLenXY;          //!<! hist. to check topological cuts
  UInt_t fFilterMask; /// FilterMask
  AliESDtrackCuts* fTrackCutsAll; //// track selection
  AliESDtrackCuts* fTrackCutsPion; /// pion track selection
  AliESDtrackCuts* fTrackCutsKaon; /// kaon track selection
  Int_t fCutTPCSignalN;   /// min. value of number of TPC clusters for PID, cut if !=0 
  Bool_t fFillHistosVsCosThetaStar; /// flag to control cos(theta*) cut
  Bool_t fApplyCutCosThetaStar; /// flag to control cos(theta*) cut
  Double_t fCutCosThetaStar;    /// cos(theta*) cut
  Bool_t fUseDzeroTopologicalCuts;  /// flag to eanble D0 topological cuts
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
  
  Int_t fMeson;                    /// mesonSpecies (see enum)
  Double_t fMassMeson;             /// mass of the selected meson
  Bool_t  fReadMC;                 ///  flag for access to MC
  Bool_t  fEnforceMBTrigMaskInMC;  /// if true force the MC to use
  Bool_t fGoUpToQuark;             /// flag for definition of c,b origin
  Int_t fFullAnalysis;             /// flag to set analysis level (0 is the fastest)
  Bool_t fSignalOnlyMC;            /// flag to speed up the MC 
  Bool_t   fSelectPtHardRange;     /// flag to select specific phard range in MC
  Double_t fMinPtHard;             /// minimum pthard
  Double_t fMaxPtHard;             /// maximum pthard
  Bool_t fRejectGeneratedEventsWithPileup;  /// reject events with generated pileup
  Bool_t fRejectSignalsFromOOBPileupEvents; /// reject signals from OOB pileup
  
  Int_t    fPIDstrategy;           /// knSigma, kBayesianMaxProb, kBayesianThres
  Double_t fmaxPforIDPion;         /// flag for upper p limit for id band for pion
  Double_t fmaxPforIDKaon;         /// flag for upper p limit for id band for kaon
  Bool_t   fKeepNegID;             /// flag to keep also track with negative ID (default kFALSE, change it only if you know what you are doing)
  Int_t    fPIDselCaseZero;        /// flag to change PID strategy
  Double_t fBayesThresKaon;        /// threshold for kaon identification via Bayesian PID
  Double_t fBayesThresPion;        /// threshold for pion identification via Bayesian PID
  Double_t fBayesThresProton;      /// threshold for proton identification via Bayesian PID
  Int_t  fOrigContainer[200000];   /// container for checks

  Int_t fDoEventMixing;            /// flag for event mixing
  Int_t  fNumberOfEventsForMixing; /// maximum number of events to be used in event mixing
  Double_t fMaxzVertDistForMix;    /// cut on zvertex distance for event mixing with cuts
  Double_t fMaxMultDiffForMix;     /// cut on multiplicity difference for event mixing with cuts
  Int_t fNzVertPools;              /// number of pools in z vertex for event mixing
  Int_t fNzVertPoolsLimSize;       /// number of pools in z vertex for event mixing +1
  Double_t* fzVertPoolLims;        //[fNzVertPoolsLimSize] limits of the pools in zVertex
  Int_t fNMultPools;               /// number of pools in multiplicity for event mixing
  Int_t fNMultPoolsLimSize;        /// number of pools in multiplicity for event mixing +1
  Double_t* fMultPoolLims;         //[fNMultPoolsLimSize] limits of the pools in multiplicity
  Int_t  fNOfPools;                /// number of pools
  TTree** fEventBuffer;            //!<! structure for event mixing
  TObjString* fEventInfo;          /// unique event Id for event mixing checks
  Double_t fVtxZ;                  /// zVertex
  Double_t fMultiplicityEM;        /// multiplicity for ev mix pools
  Double_t fMultiplicityMC;        /// multiplicity for MC efficiencies
  Int_t fMultEstimMC;              /// multiplicity estimator (0=tracklets in eta<1)
  Int_t fNumOfMultBins;            /// number of bins for multiplcities in MC histos
  Double_t fMinMultiplicity;       /// lower limit for multiplcities in MC histos
  Double_t fMaxMultiplicity;       /// upper limit for multiplcities in MC histos
  TObjArray* fKaonTracks;          /// array of kaon-compatible tracks (TLorentzVectors)
  TObjArray* fPionTracks;          /// array of pion-compatible tracks (TLorentzVectors)
    
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCombinHF,41); /// D0D+ task from AOD tracks
  /// \endcond
};

#endif
