#ifndef ALIANALYSISTASK_JPSIJET_H
#define ALIANALYSISTASK_JPSIJET_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
 *  Analysis Task for J/psi production in jets
 * 
 * Pre-requisities
 *  - Physics selection
 *  - PID response
 *  - CDB connect
 *  - EMCal correction
 * 
 * Dependencies
 *  - PWGDQ/dielectron
 *  - PWGJE/EMCalJet
 *  - AOD Extension + Branch Replicator
 * 
 * Supported
 *  - Event: AOD
 *  - Period: Run2 pp 13TeV
 * 
 * Outputs
 *  - AliAOD.Dielectron.root: filtered AOD event with dielectron and jets
 *  - JpsiJetResult.root: histograms for QA and physic objects
 * 
 * Author: [Yìtāo WÚ](yitao.wu@cern.ch), USTC
 * Create Date: 2019.08.17
 * More details: [Analysis Repo](https://github.com/yatowoo/JpsiJet) 
**/

#include "THistManager.h"

#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAnalysisCuts.h"
#include "AliDielectron.h"
#include "AliDielectronVarManager.h"
#include "AliJetContainer.h"
#include "AliFJWrapper.h"
class AliFJWrapper;
#include "AliEmcalJetTask.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisCuts;
class AliDielectron;
class AliDielectronVarManager;
class AliEmcalJet;
class AliJetContainer;

class AliAnalysisTaskJpsiJet : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskJpsiJet();
  AliAnalysisTaskJpsiJet(const char *taskName);
  virtual ~AliAnalysisTaskJpsiJet();
  // Interface for anslysis manager (based on TaskSE)
public:
  virtual void LocalInit();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * /*option*/);
  virtual void UserExecMix(Option_t * /*option*/) { ; }
  virtual Bool_t UserNotify() { return kTRUE; }
  virtual void NotifyRun() { ; }
  virtual void Terminate(Option_t *option="");

  // QA Histograms
public:
  void   SetTriggerQA(Bool_t enableTriggerQA = kTRUE){ fIsTriggerQA = enableTriggerQA;}
  Bool_t GetTriggerQA(){return fIsTriggerQA;}
  void   SetCellQA(Bool_t enableCellQA = kTRUE){ fIsCellQA = enableCellQA;}
  Bool_t GetCellQA(){return fIsCellQA;}
private:
  void InitHistogramsForEventQA(const char* histClass);
  void InitHistogramsForClusterQA(const char* histClass);
  void InitHistogramsForJetQA();
  void FillHistogramsForEventQA(const char* histClass);
  void FillHistogramsForClusterQA(const char* histClass);
  void FillHistogramsForJetQA(const char* jetTag);
  void FillTH2(const char* histName, const char* labelX, Double_t value, Double_t weight = 1.0);
  void InitHistogramsForRunwiseQA(const char* histClass);
  void FillHistogramsForRunwiseQA(const char* histClass);
  void FillHistogramsForRunwiseElectronQA(const char* histClass, AliAODTrack* ele);
  
  // Dielectron framework
private:
  void InitDielectron();
  void InitHistogramsForDielectron(const char* histMgrName);
  Bool_t FindDaughters(AliVTrack *trk);
  Double_t GetPseudoProperDecayTime(AliDielectronPair* pair);
  void AddTrackFromPair(AliAODTrack *trkTemplate);
  Bool_t RunDetectoreLevelAnalysis();

  // FF analysis
private:
  Bool_t FindTrackInJet(AliEmcalJet* jet, AliVParticle* p, Int_t trackID);
  void InitHistogramsForTaggedJet(const char* histClass);
  Bool_t FillHistogramsForTaggedJet(const char* histClass);

  // Jet finder task
  // -- These sub-tasks are managed by this task, not by AliAnalysisManager.
public:
  void AddTaskEmcalJet(
    const TString nTracks = "usedefault",
    const TString nClusters = "",
    const AliJetContainer::EJetAlgo_t jetAlgo = AliJetContainer::antikt_algorithm,
    const Double_t radius = 0.4,
    const AliJetContainer::EJetType_t jetType = AliJetContainer::kFullJet,
    const Double_t minTrPt = 0.15,
    const Double_t minClPt = 0.30,
    const Double_t ghostArea = 0.005,
    const AliJetContainer::ERecoScheme_t reco = AliJetContainer::pt_scheme,
    const TString tag = "Jet",
    const Double_t minJetPt = 0.0,
    const Bool_t lockTask = kTRUE,
    const Bool_t bFillGhosts = kFALSE);
  void InitJetFinders();
  Bool_t RunJetFinder(const char* jetTag);
  void EnableJetFinder(Bool_t enableJetFinder = kTRUE){fIsJetFinder = enableJetFinder;}
  Bool_t CheckJetFinder(){return fIsJetFinder;}

  // MC correction
public:
  void   SetMC(Bool_t isMC = kTRUE){fIsMC = isMC;}
  Bool_t IsMC(){return fIsMC;}
private:
  void   InitHistogramsForMC();
  void   InitHistogramsForJpsiMC(const char* histClass);
  void   InitHistogramsForJetMC();
  Bool_t ApplyEmcalCut(AliVParticle* par, Bool_t isMCTruth);
  Bool_t RunParticleLevelAnalysis();
  void   SetJpsiGeneratorType();
  void   FillHistogramsForParticle(const char* histName, AliVParticle* par);
  void   FillHistogramsForElectronPID(const TObjArray* eleArray);
  void   FillHistogramsForJpsiMC();
  Double_t GetJetMCPt(AliEmcalJet* jet);
  void   FillHistogramsForJetMC(const char* jetTag);
  Int_t  CheckDielectronDaughter(AliVParticle* par);
  AliAODMCParticle* AddParticleFromJpsi(AliAODMCParticle* jpsi);
  Bool_t FillHistogramsForTaggedJetMC(AliAODMCParticle* jpsi);

  // Event selection
public:
  void     SetTrigger(UInt_t trigger){fSelectedTrigger = trigger;}
  UInt_t GetTrigger(){return fSelectedTrigger;}
  void     SetTriggerClasses(TString trigClasses){fSelectedTriggerClasses = trigClasses;}
  TString  GetTriggerClasses(){return fSelectedTriggerClasses;}
  TString  GetTriggerTag(){return fFiredTriggerTag;}
  void     SetRejectPileup(Bool_t rejectPileup){fRejectPileup = rejectPileup;}
  Bool_t   GetRejectPileup(){return fRejectPileup;}
  void     SetEventFilter(AliAnalysisCuts *eventCuts){fEventFilter = eventCuts;}
  void     PrintEventFilter(){;}
  void     SetEMCalEnergyThreshold(Double_t eth){fEMCalEth = eth;}
  Double_t GetEMCalThreshold(){return fEMCalEth;}

/**
 *  Data members 
 **/
public:
enum EventStatus{ // for histogram event stats
  kAllInAOD,
  kPhysSelected,
  kFiltered,
  kV0ANDtrigger,
  kTRDtrigger,
  kAfterPileUp,
  kWithSinglePair,
  kWithMultiPair,
  kWithPairInJet,
  kEventStatusN
};

private:
enum PdgCode{
  PDG_PROTON = 2212,
  PDG_ELECTRON = 11,
  PDG_GAMMA = 22,
  PDG_PION = 211,
  PDG_KAON = 321,
  PDG_JPSI = 443
};

private:
  AliAODEvent       *fAOD; // Input AOD event
  AliDielectron     *fDielectron; // Dielectron framework
  TClonesArray      *fPairs; // Dielectron pair
  TClonesArray      *fDaughters; // Dielectron tracks
  TClonesArray      *fTracksWithPair; // Track array with daughters replaced by pair
  TObjArray         *fJetTasks; // Jet finder tasks
  TObjArray         *fJets; // Jet containers
  UInt_t             fSelectedTrigger; // Event offline trigger
  TString            fSelectedTriggerClasses; // Event fired trigger classes (separated by ';')
  TString            fFiredTriggerTag; // MB, EG1, EG2, DG1, DG2 
  Bool_t             fRejectPileup;
  Bool_t             fIsPileup;
  Bool_t             fIsTriggerQA; // Add QA histograms by trigger classes
  Bool_t             fIsCellQA; // Add QA histograms by Cell ID
  Bool_t             fIsJetFinder; // ONLY 1 task should enable the jet finder in LEGO train
  Bool_t             fIsMC;  // Input MC production
  Double_t           fEMCalEth; // EMCal gamma trigger energy threshold
  TClonesArray      *fMCParticles; // AOD event object in MC
  AliAODMCHeader    *fMCHeader;    // AOD event object in MC
  TString            fMCGenType;   // Generator type - Prompt/Jpsi2ee, Bdecay/B2Jpsi2ee
  AliAnalysisCuts   *fEventFilter;

// J/psi tagged jet
  AliVParticle      *fJpsiPair;
  AliEmcalJet       *fTaggedJet;
  AliVParticle      *fJpsiMC;
  AliEmcalJet       *fTaggedJetMC;

// Histograms
  THistManager      *fHistos;
  THistManager      *fHistosMC;

// Runwise QA
private:
  TString            fRunNo;
  Double_t          *fVars;

private:
  ClassDef(AliAnalysisTaskJpsiJet, 6);
};

#endif /* ALIANALYSISTASK_JPSIJET_H */
