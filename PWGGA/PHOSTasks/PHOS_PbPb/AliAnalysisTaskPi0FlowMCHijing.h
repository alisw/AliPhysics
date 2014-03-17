/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Authors: Henrik Qvigstad
// Date   : 20.06.2013
/* $Id$ */

#ifndef ALIANALYSISTASKPI0FLOWMCHIJING_H
#define ALIANALYSISTASKPI0FLOWMCHIJING_H

#include "AliAnalysisTaskPi0Flow.h"
#include "AliAnalysisTaskPi0FlowMC.h"

class AliAnalysisTaskPi0FlowMCHijing :  public AliAnalysisTaskPi0FlowMC
{
public:
  AliAnalysisTaskPi0FlowMCHijing(const char* name = "AliAnalysisTaskPi0Flow", Period period = AliAnalysisTaskPi0Flow::kUndefinedPeriod);
  virtual ~AliAnalysisTaskPi0FlowMCHijing();

protected: // member functions:
  AliAnalysisTaskPi0FlowMCHijing(const AliAnalysisTaskPi0FlowMC&); // not implemented
  AliAnalysisTaskPi0FlowMCHijing& operator=(const AliAnalysisTaskPi0FlowMC&); // not implemented

  virtual Double_t PrimaryWeight(Int_t primary);
  virtual Double_t PrimaryParticleWeight(TParticle * particle);

  ClassDef(AliAnalysisTaskPi0FlowMCHijing, 1); // PHOS analysis task
};

#endif // ALIANALYSISTASKPI0FLOWMCHIJING_H
