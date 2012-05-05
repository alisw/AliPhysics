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
class TList;
class AliDxHFEParticleSelection;
class AliDxHFECorrelation;

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
  void SetOption(const char* opt) { fOption = opt; }
  /// overloaded from TObject: get option
  virtual Option_t* GetOption() const { return fOption;}

 protected:

 private:
  /// copy constructor prohibited: might change
  AliAnalysisTaskDxHFECorrelation(const AliAnalysisTaskDxHFECorrelation&);
  /// assignment operator prohibited: might change
  AliAnalysisTaskDxHFECorrelation& operator=(const AliAnalysisTaskDxHFECorrelation&);

  int DefineSlots();

  TList* fOutput;                        //! list send on output slot 1
  TString fOption;                       // option string
  AliDxHFECorrelation* fCorrelation;     // correlation worker class
  AliDxHFEParticleSelection* fD0s;       // selection of D0s
  AliDxHFEParticleSelection* fElectrons; // selection of electrons

  ClassDef(AliAnalysisTaskDxHFECorrelation, 1);
};

#endif
