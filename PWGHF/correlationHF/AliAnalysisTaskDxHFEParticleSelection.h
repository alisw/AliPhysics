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
class TList;
class AliDxHFEParticleSelection;

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

 protected:

 private:
  /// copy constructor prohibited: might change
  AliAnalysisTaskDxHFEParticleSelection(const AliAnalysisTaskDxHFEParticleSelection&);
  /// assignment operator prohibited: might change
  AliAnalysisTaskDxHFEParticleSelection& operator=(const AliAnalysisTaskDxHFEParticleSelection&);

  int DefineSlots();

  TList* fOutput;                  //! list send on output slot 1
  TString fOption;                 // option string
  AliDxHFEParticleSelection* fSelector; // selector instance

  ClassDef(AliAnalysisTaskDxHFEParticleSelection, 1);
};

#endif
