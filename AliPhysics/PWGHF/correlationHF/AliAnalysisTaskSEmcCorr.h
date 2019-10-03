#ifndef ALIANALYSISTASKSEMCCORR_H
#define ALIANALYSISTASKSEMCCORR_H

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
class TString;

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
  void SetCheckDecay(Bool_t checkDecay){fCheckDecay=checkDecay;}
  void SetDoHFCorrelations(Bool_t doHFcor){fDoHFCorrelations=doHFcor;}
  void SetIsFastSimulation(Bool_t fastSimul){fFastSimul=fastSimul;}
  void CheckDzeroChannels(AliAODMCParticle *part,TClonesArray *mcarray);
  void CheckDplusChannels(AliAODMCParticle *part,TClonesArray *mcarray);
  //  void CheckBallChannels(AliAODMCParticle *part,TClonesArray *mcarray);
  //   void CheckLambdaChannels(AliAODMCParticle *part,TClonesArray *mcarray);
  void GetFinalStateParticles(AliAODMCParticle *part,TClonesArray *clarr);
  void SetGeneratorToBeChecked(TString strgen){fGeneratorString=strgen;}

 private:
  AliAnalysisTaskSEmcCorr(const AliAnalysisTaskSEmcCorr&); // copy constructo not implemented yet
  AliAnalysisTaskSEmcCorr& operator=(const AliAnalysisTaskSEmcCorr&); // assignment operator not implemented yet
  Bool_t fReadMC;                            // read MC flag
  Bool_t fFastSimul;                         // flag to access the MC array with standard name (kFALSE) or with the "magic number" (23) for the fast simulation
  Bool_t fCheckDecay;                       // check decay properties
  Bool_t fDoHFCorrelations;                 // activate checks for D0-h correlations
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
  TH1F* fhDzeroDecay;                           //! histo with Dzero decay channels counting
  TH1F* fhDplusDecay;                           //! histo with Dplus decay channels counting 
  TH1F* fhLambdaCDecay;                        //! histo with LambdaC decay channels counting
  TH1F* fhAllBDecay;                            //! histo with AllB decay channels counting
  Int_t fnpionPlus;                         //! internal variable
  Int_t fnpionMinus;                       //! internal variable
  Int_t fnpionZero;                        //! internal variable 
  Int_t fnkaonPlus;                        //! internal variable
  Int_t fnkaonMinus;                       //! internal variable
  Int_t fnkaonZeroShort;                        //! internal variable
  Int_t fnProton;                         //! internal variable
  Int_t fnAntiProton;                     //! internal variable                      
  Int_t fnElePlus;                        //! internal variable
  Int_t fnEleMinus;                        //! internal variable
  Int_t fnMuonPlus;                        //! internal variable
  Int_t fnMuonMinus;                        //! internal variable
  Int_t fnNeutrino;                        //! internal variable (counting both ele and muon neutrinos)
  Int_t fnJpsi;                            //! internal variable
  Int_t fnPhysPrim;                        //! internal variable to count number of daughters (phys prim + pi0+ J/Psi)
  Double_t fpxtotdaugh;                     //! internal variable for checking total momentum of counted daughter
  Double_t fpytotdaugh;                       //! internal variable for checking total momentum of counted daughter
  Double_t fpztotdaugh;                      //! internal variable for checking total momentum of counted daughter
  TString fGeneratorString;                 // generator title that for event generator selection
  ClassDef(AliAnalysisTaskSEmcCorr,2); // analysis task for MC study

};

#endif
