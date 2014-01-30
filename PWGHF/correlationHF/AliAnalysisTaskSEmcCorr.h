#ifndef ALIANALYSISTASKSEMCCORR_H
#define ALIANALYSISTASKSEMCORR_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEmcCorr
// AliAnalysisTaskSE for studying HF-(hadron,electrons) and hadron-hadron correlations
//     at MC level
// Author: Andrea Rossi, andrea.rossi@cern.ch
//
//*************************************************************************

class TH1F;
class TH2F;
class TClonesArray;
class AliAODDEvent;
class AliAODMCHeader;
class AliAODRecoDecayHF2Prong;
class AliAODRecoDecayHF;
class AliAODMCParticle;
class AliAnalysisVertexingHF;
class AliRDHFCutsD0toKpi;
class AliNormalizationCounter;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSEmcCorr : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSEmcCorr();
  AliAnalysisTaskSEmcCorr(const char* name);
  virtual ~AliAnalysisTaskSEmcCorr();

 // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);  
  void SetReadMC(Bool_t readMC){fReadMC=readMC;}
  void SetKineOnly(Bool_t readKineOnly){fOnlyKine=readKineOnly;}
  void DoHadronHadron(Bool_t dohh){fDoHadronHadron=dohh;}
  void SetMinPtD(Double_t ptmin){fminDpt=ptmin;}
  void SetMaxYtrigger(Double_t maxy){fYrangeTrig=maxy;}
  void SetMaxEtaEleTrigger(Double_t maxeta){fEtarangeEleTrig=maxeta;}
  void SelectAssociatedParticles(TClonesArray *arrayMC);
  void SetEtaMaxAssPart(Double_t eta){fmaxEtaAss=eta;}
  void FillCorrelationPlots(AliAODMCParticle *part,TClonesArray *arrayMC);
  void FillCorrelationPlotsHadrons(AliAODMCParticle *part,TClonesArray *arrayMC);
  void FillSkipParticleArray(AliAODMCParticle *part,TClonesArray *arrayMC,Int_t partPos);

 private:
  AliAnalysisTaskSEmcCorr(const AliAnalysisTaskSEmcCorr&); // copy constructo not implemented yet
  AliAnalysisTaskSEmcCorr& operator=(const AliAnalysisTaskSEmcCorr&); // assignment operator not implemented yet
  Bool_t fReadMC;                            // read MC flag
  Bool_t fOnlyKine;                          // kine only
  Bool_t fDoHadronHadron;                    // flag to activate h-h correlation
  TH1F *fNentries;                           //  histo with event selection properties
  TH1F *fhNprongsD0;                         //  histo with D0 nprongs
  TH1F *fhNprongsD0chargedOnly;                         //  histo with D0 nprongs
  TH1F *fhNprongsD0chargedRef;                         //  histo with D0 nprongs
  Double_t fYrangeTrig;                       // max D y range
  Double_t fEtarangeEleTrig;                       // max eta ele
  Double_t   fminDpt;                       // min D pt
  TArrayI *fArrayDaugh;                       //! array with index of D daughter particles 
  Int_t  flastdaugh;                       // last filled position of farraydaugh 
  TArrayI *fArrayAssoc;                       //! array with index of associated particle
  Int_t  fLastAss;                       // last filled posistion of farrayassoc
  Double_t   fminPtass;                       // min pt assoc
  Double_t   fmaxEtaAss;                       // max eta assoc
  THnSparseF*  fhMCcorrel;                       //!  histo with correlations
  THnSparseF*  fhMCtrigPart;                       //! histo for counting trigger particles
  THnSparseF*  fhMChadroncorrel;                       //!  histo with correlations for hadron-hadron
  THnSparseF*  fhMChadrontrigPart;                       //! histo for counting trigger particles for hadron-hadron correl

  ClassDef(AliAnalysisTaskSEmcCorr,1); // analysis task for MC study

};

#endif
