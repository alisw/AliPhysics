//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliAnalysisTaskDxHFEParticleSelection.h
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  AnalysisTask electron selection for D0 - HFE correlation
///

#ifndef ALIANALYSISTASKDXHFEPARTICLESELECTION_H
#define ALIANALYSISTASKDXHFEPARTICLESELECTION_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"

class AliDxHFEParticleSelection;
class AliAnalysisCuts;
class TList;
class TObjArray;

/**
 * @class AliAnalysisTaskDxHFEParticleSelection
 * Selection task for particles uesd in the D0 - HFE correlation studies
 * Task performs the selection based on a configured AliDxHFEParticleSelection
 * instance.
 */
class AliAnalysisTaskDxHFEParticleSelection : public AliAnalysisTaskSE {
  public:
  /// constructor
  AliAnalysisTaskDxHFEParticleSelection(const char* opt="");
  /// destructor
  virtual ~AliAnalysisTaskDxHFEParticleSelection();

  enum {
    kD0=0,
    kElectron=1
  };

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

  void SetOption(const char* opt) { fOption = opt; }
  void SetFillOnlyD0D0bar(Int_t flagfill){fFillOnlyD0D0bar=flagfill;}
  void SetParticleType(int particle){fParticleType=particle;}
  virtual void SetUseMC(Bool_t useMC){fUseMC=useMC;}
  //  virtual void SetCuts(AliAnalysisCuts* cuts){fCuts=cuts;}
  virtual void SetCutList(TList* cuts){fCutList=cuts;}
  Bool_t GetUseMC() const {return fUseMC;}

 protected:

 private:
  /// copy constructor prohibited: might change
  AliAnalysisTaskDxHFEParticleSelection(const AliAnalysisTaskDxHFEParticleSelection&);
  /// assignment operator prohibited: might change
  AliAnalysisTaskDxHFEParticleSelection& operator=(const AliAnalysisTaskDxHFEParticleSelection&);

  int DefineSlots();
  int ParseArguments(const char* arguments);

  TList* fOutput;                       // list send on output slot 1
  TString fOption;                      // option string
  TList* fCutList;                         // TList containg cut objects
  AliAnalysisCuts *fCutsD0;             // Cut Object for D0 
  AliDxHFEParticleSelection* fSelector; // selector instance
  bool fUseMC;                          // use MC info
  Int_t     fFillOnlyD0D0bar;           // flag to set what to fill (0 = both, 1 = D0 only, 2 = D0bar only)
  TObjArray *fSelectedTracks;           // Array for selected Tracks
  TObjArray *fMCArray;                  // MC array
  Int_t fParticleType;                   // Holds which particle to run on
  Int_t fSystem;                        // holds collisions system (0=pp, 1=PbPb(,2=pPb))
  TString fSelectionParticleOptions;    // String to hold options for the particle selection
  Bool_t fUseKine;                      // Whether or not to run on MC stack


  ClassDef(AliAnalysisTaskDxHFEParticleSelection, 4);
};

#endif
