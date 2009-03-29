#ifndef ALIANALYSISTASKGAMMACONVERSION_H
#define ALIANALYSISTASKGAMMACONVERSION_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to do analysis on conversion pairs
//---------------------------------------------
////////////////////////////////////////////////
 
#include "AliAnalysisTaskSE.h"
#include <vector>
#include "AliV0Reader.h"
#include "TNtuple.h"

class AliGammaConversionHistograms;
class AliESDv0;
class AliKFParticle;
class AliESDInputHandler;
class AliESDEvent;
class AliAODEvent;
class TList;
class AliStack;

class AliAnalysisTaskGammaConversion : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskGammaConversion();
  AliAnalysisTaskGammaConversion(const char* name);
  virtual ~AliAnalysisTaskGammaConversion() ;// virtual destructor
 
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void ConnectInputData(Option_t *);
	
  void ProcessMCData();
  void ProcessV0sNoCut();
  void ProcessV0s();
  void ProcessGammasForNeutralMesonAnalysis();
  void SetHistograms(AliGammaConversionHistograms *histograms){fHistograms=histograms;}
  void SetDoMCTruth(Bool_t flag){fDoMCTruth=flag;}
  void SetElectronMass(Double_t electronMass){fElectronMass = electronMass;}
  void SetGammaMass(Double_t gammaMass){fGammaMass = gammaMass;}
  void SetGammaWidth(Double_t gammaWidth){fGammaWidth = gammaWidth;}
  void SetPi0Mass(Double_t pi0Mass){fPi0Mass = pi0Mass;}
  void SetPi0Width(Double_t pi0Width){fPi0Width = pi0Width;}
  void SetEtaMass(Double_t etaMass){fEtaMass = etaMass;}
  void SetEtaWidth(Double_t etaWidth){fEtaWidth = etaWidth;}
  void SetMinOpeningAngleGhostCut(Double_t ghostCut){fMinOpeningAngleGhostCut = ghostCut;}
  void SetV0Reader(AliV0Reader* reader){fV0Reader=reader;}
  void SetCalculateBackground(Bool_t bg){fCalculateBackground=bg;}
  void CalculateBackground();
  void SetWriteNtuple(Bool_t writeNtuple){fWriteNtuple = writeNtuple;}
  void FillNtuple();
  Double_t GetMCOpeningAngle(TParticle* daughter0, TParticle* daughter1) const;

 private:
  AliAnalysisTaskGammaConversion(const AliAnalysisTaskGammaConversion&); // Not implemented
  AliAnalysisTaskGammaConversion& operator=(const AliAnalysisTaskGammaConversion&); // Not implemented

  AliV0Reader* fV0Reader;

  AliStack * fStack;

  TList * fOutputContainer ; // Histogram container

  AliGammaConversionHistograms *fHistograms;

  Bool_t fDoMCTruth;
    
  vector<TParticle*> fMCAllGammas;
  vector<TParticle*> fMCPi0s;
  vector<TParticle*> fMCEtas;
  vector<TParticle*> fMCGammaChic;

  vector<AliKFParticle> fKFReconstructedGammas;
  vector<Bool_t> fIsTrueReconstructedGammas;
  vector<Int_t> electronv1;
  vector<Int_t> electronv2;

  //mass defines
  Double_t fElectronMass;
  Double_t fGammaMass;
  Double_t fPi0Mass;
  Double_t fEtaMass;

  // width defines
  Double_t fGammaWidth;
  Double_t fPi0Width;
  Double_t fEtaWidth;

  Double_t fMinOpeningAngleGhostCut;

  Bool_t fCalculateBackground;
  Bool_t fWriteNtuple;
  TNtuple *fGammaNtuple;
  TNtuple *fNeutralMesonNtuple;

  Int_t fTotalNumberOfAddedNtupleEntries;

  ClassDef(AliAnalysisTaskGammaConversion, 2); // Analysis task for gamma conversions
};
 
#endif //ALIANALYSISTASKGAMMA_H
