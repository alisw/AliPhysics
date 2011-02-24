#ifndef ALIANALYSISTASKSEDSTARJETS_H
#define ALIANALYSISTASKSEDSTARJETS_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */ 

//-----------------------------------------------------------------------
// Author : A. Grelli, UTRECHT
//-----------------------------------------------------------------------


#include <TH2F.h>
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"

class TH2F;
class TH1I;
class TParticle ;
class TFile ;
class TClonesArray ;
class AliCFManager;
class AliAODRecoDecay;
class AliAODRecoDecayHF2Prong;
class AliAODMCParticle;


class AliAnalysisTaskSEDStarJets : public AliAnalysisTaskSE {
  
 public:
  
  AliAnalysisTaskSEDStarJets();
  AliAnalysisTaskSEDStarJets(const Char_t* name);
  AliAnalysisTaskSEDStarJets& operator= (const AliAnalysisTaskSEDStarJets& c);
  AliAnalysisTaskSEDStarJets(const AliAnalysisTaskSEDStarJets& c);
  virtual ~AliAnalysisTaskSEDStarJets();
  
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  
  // User functions

  Double_t GetInvariantMass(TLorentzVector LorentzTrack1, TLorentzVector LorentzTrack2);
  Double_t GetInvariantMassDStar(TLorentzVector LorentzTrack3,TLorentzVector LorentzTrack4);
 
  //side band background eval
  void     SideBandBackground(Double_t finvM, Double_t finvMDStar, Double_t fejet, Double_t ejet, Int_t nJets);
  
  // inizializations
  Bool_t   DefineHistoFroAnalysis();
  
  //MC values for D0 and D*
  
  Bool_t   DstarInMC(AliAODMCParticle* const mcPart, TClonesArray* mcArray);
  Bool_t   EvaluateIfD0toKpi(AliAODMCParticle* neutralDaugh, TClonesArray* mcArray)const;

  // Alternative cut method
  Bool_t   EvaluateCutOnPiD0pt(AliAODRecoDecayHF2Prong* const vtx, AliAODTrack* const aodTrack);
  // set minimum ITS clusters for the analysis
  void     SetMinITSClusters(Int_t minITSClusters) {fMinITSClusters = minITSClusters;}
  Int_t    GetMinITSClusters() const {return fMinITSClusters;}

  //set the analysis type D*+ or D*-
  void     SetAnalType(Bool_t computeD0) {fComputeD0 = computeD0;}
  Bool_t   GetAnalType() const {return fComputeD0;}

  // set MC usage
  void    SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t  GetMC() const {return fUseMCInfo;}

  // set cut type
  void     SetCutType(Bool_t topologicalCut) {ftopologicalCut = topologicalCut;}
  Bool_t   GetCutType() const {return ftopologicalCut;}
  
 protected:
  
  Int_t  fCountReco;             //  Reco particle found that satisfy cuts
  Int_t  fCountRecoAcc;          //  Reco particle found that satisfy cuts in requested acceptance
  Int_t  fCountRecoITSClusters;  //  Reco particle found that satisfy cuts in n. of ITS clusters
  Int_t  fCountRecoPPR;          //  Reco particle found that satisfy cuts in PPR
  Int_t  fCountDStar;            //  MC particle that are D* in acc and with D0->kpi.
  Int_t  fEvents;                //  n. of events
  Int_t  fMinITSClusters;        //  min n. of ITS clusters for RecoDecay
  Bool_t fComputeD0;             //  select analysis type: D*+ (kTRUE), D*- (kFALSE)
  Bool_t fUseMCInfo;             //  Use MC info
  Bool_t ftopologicalCut;        //  if false apply relaxed PPR cuts alse cut on the space of D0pt and softpipt  
  Bool_t fRequireNormalization;  //  normalization 

  TLorentzVector fLorentzTrack1; // lorentz 4 vector
  TLorentzVector fLorentzTrack2; // lorentz 4 vector
  TLorentzVector fLorentzTrack3; // lorentz 4 vector
  TLorentzVector fLorentzTrack4; // lorentz 4 vector

  TList *fOutput;         //! user output

  // define the histograms 
  // 2D
  TH2F *fD0ptvsSoftPtSignal;    //!
  TH2F *fD0ptvsSoftPt;          //!
 
  //1D
  TH1F *ftrigger;        //!
  TH1F *fPtPion;         //!
  TH1F *fInvMass;        //!
  TH1F *fRECOPtDStar;    //!
  TH1F *fDStar;          //!
  TH1F *fDiff;           //!
  TH1F *fDiffSideBand;   //!
  TH1F *fDStarMass;      //!
  TH1F *fPhi;            //!
  TH1F *fPhiBkg;         //!
  TH1F *fTrueDiff;       //!
  TH1F *fResZ;           //!
  TH1F *fResZBkg;        //!
  TH1F *fcharmpt;        //!
  TH1F *fdstarE;         //!
  TH1F *fEjet;           //!
  TH1F *fPhijet;         //!
  TH1F *fEtaJet;         //!
  TH1F *fdstarpt;        //!

  ClassDef(AliAnalysisTaskSEDStarJets,2); // class for HF corrections as a function of many variables
};

#endif
