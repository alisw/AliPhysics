#ifndef ALIANALYSISTASKESDMCLABELADDITION_H
#define ALIANALYSISTASKESDMCLABELADDITION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

#include <TString.h>
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskESDMCLabelAddition : public AliAnalysisTaskSE
{
  
public:
  AliAnalysisTaskESDMCLabelAddition();
  AliAnalysisTaskESDMCLabelAddition(const char* name);
  virtual ~AliAnalysisTaskESDMCLabelAddition() {;}
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void NotifyRun();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fDefaultStorage = ocdbPath; }
  
  
private:
  
  AliAnalysisTaskESDMCLabelAddition(const AliAnalysisTaskESDMCLabelAddition&);
  AliAnalysisTaskESDMCLabelAddition& operator=(const AliAnalysisTaskESDMCLabelAddition&);
  
  TString  fDefaultStorage; ///< location of the default OCDB storage
  Double_t fSigmaCut;       //!< sigma cut to associate clusters with TrackRefs
  Double_t fSigmaCutTrig;   //!< sigma cut to associate trigger track to triggerable track
  
  ClassDef(AliAnalysisTaskESDMCLabelAddition, 2); // Analysis task for standard ESD filtering
  
};

#endif

