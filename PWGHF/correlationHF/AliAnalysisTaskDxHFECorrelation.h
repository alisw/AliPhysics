//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliAnalysisTaskDxHFECorrelation.h
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  AnalysisTask D0 - HFE correlation
///

#ifndef ALIANALYSISTASKDXHFECORRELATION_H
#define ALIANALYSISTASKDXHFECORRELATION_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"

class AliPID;
class AliPIDResponse;
class TList;
class AliDxHFEParticleSelection;
class AliDxHFEParticleSelectionD0;
class AliDxHFEParticleSelectionEl;
class AliDxHFECorrelation;
class AliAnalysisCuts;
class AliHFEpid;
class AliHFEcuts;
class AliHFAssociatedTrackCuts;

/**
 * @class AliAnalysisTaskDxHFECorrelation
 * Task for D0-HFE correleations
 */
class AliAnalysisTaskDxHFECorrelation : public AliAnalysisTaskSE {
  public:
  /// constructor
  AliAnalysisTaskDxHFECorrelation(const char* opt="");
  /// destructor
  virtual ~AliAnalysisTaskDxHFECorrelation();

  /// inherited from AliAnalysisTask: connect tree branches at input slots
  virtual void ConnectInputData(Option_t *option="") {
    return AliAnalysisTaskSE::ConnectInputData(option);
  }

  /// inherited from AliAnalysisTaskSE: create output objects
  virtual void UserCreateOutputObjects();
  /// inherited from AliAnalysisTaskSE: event processing
  virtual void UserExec(Option_t*);
  /// inherited from AliAnalysisTask: called in SlaveTerminate phase for each task
  virtual void FinishTaskOutput();
  /// inherited from AliAnalysisTask: final step
  virtual void Terminate(Option_t*);

  /// set options
  // TODO: Some of them are not in use, as the members are set by parsing arguments.
  // Keep it for now.
  void SetOption(const char* opt) { fOption = opt; }
  virtual void SetUseMC(Bool_t useMC){fUseMC=useMC;}
  virtual void SetCutsD0(AliAnalysisCuts* cuts){fCutsD0=cuts;}
  virtual void SetCutsHFE(TList* cuts){fListHFE=cuts;}

  void SetCuts(AliAnalysisCuts* cuts){fCuts=cuts;}
  void SetUseEventMixing(Bool_t useMixing) {fUseEventMixing=useMixing;}
  void SetSystem(Bool_t system){fSystem=system;}
  void SetTriggerParticle(int trigger){fTriggerParticle=trigger;} 
  void SetD0EffMap(TH1* eff, int which=kPrompt){
    if(which==kPrompt) fD0EffMapP=eff;
    if(which==kFeedDown) fD0EffMapFD=eff;
  }

  /// overloaded from TObject: get option
  virtual Option_t* GetOption() const { return fOption;}
  Bool_t GetUseMC() const {return fUseMC;}

  //TODO: Use enums in AliDxHFECorrelation
  enum {
    kPrompt=0,
    kFeedDown=1
  };

 protected:

 private:
  /// copy constructor prohibited: might change
  AliAnalysisTaskDxHFECorrelation(const AliAnalysisTaskDxHFECorrelation&);
  /// assignment operator prohibited: might change
  AliAnalysisTaskDxHFECorrelation& operator=(const AliAnalysisTaskDxHFECorrelation&);

  int ParseArguments(const char* arguments);
  int DefineSlots();

  TList* fOutput;                        //! list send on output slot 1
  TString fOption;                       //  option string
  AliDxHFECorrelation* fCorrelation;     //  correlation worker class
  AliDxHFECorrelation* fCorrelationCharm;//  correlation worker class - kine only store C
  AliDxHFECorrelation* fCorrelationBeauty; //  correlation worker class - kine only store B
  AliDxHFECorrelation* fCorrelationNonHF;//  correlation worker class - kine only store NonHF
  AliDxHFECorrelation* fCorrelationHadron;//  correlation worker class - kine only store hadrons
  AliDxHFEParticleSelection* fD0s;       //  selection of D0s
  AliDxHFEParticleSelection* fElectrons; //  selection of electrons
  AliAnalysisCuts *fCutsD0;              //  Cuts D0 
  AliAnalysisCuts *fCuts;                // Cuts which holds info for AliHFCorrelator 
  Bool_t fUseMC;                         // use MC info
  Bool_t fUseEventMixing;                // Run Event Mixing analysis
  Int_t fSystem;                         // Which system pp/PbPb
  TObjArray *fSelectedD0s;               // Array for selected D0s
  TObjArray *fSelectedElectrons;         // Array for selected Electrons
  TList* fListHFE;                       // List containing cut and pid objects for HFE
  Int_t fTriggerParticle;                // Which particle to trigger on 
  Bool_t fUseKine;                       // To run over MC or reconstructed data
  TObjArray* fMCArray;                   // Array to hold MCarray
  TString fCorrelationArguments;         // String argument for correlation
  TH1* fD0EffMapP;                       // histo containing efficiency map for D0 for prompt D0
  TH1* fD0EffMapFD;                      // histo containing efficiency map for D0 for feeddown D0 (only useful for MC)
  Bool_t fStoreSeparateOrigins;          // Whether to create correlation objects for various origins

  ClassDef(AliAnalysisTaskDxHFECorrelation, 7);
};

#endif

