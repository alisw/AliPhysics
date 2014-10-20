/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Authors: Boris Polishchuk
// Date   : 09.07.2013

#ifndef ALIANALYSISTASKPI0FLOWMCPARAMWEIGHTS_H
#define ALIANALYSISTASKPI0FLOWMCPARAMWEIGHTS_H

#include "AliAnalysisTaskPi0FlowMC.h"

class AliAnalysisTaskPi0FlowMCParamWeights :  public AliAnalysisTaskPi0FlowMC
{
public:
  
  AliAnalysisTaskPi0FlowMCParamWeights(const char* name = "AliAnalysisTaskPi0Flow", Period period = AliAnalysisTaskPi0Flow::kUndefinedPeriod);
  virtual ~AliAnalysisTaskPi0FlowMCParamWeights();
  
  void UserCreateOutputObjects();
  // void SetParameters(Double_t p0,Double_t p1,Double_t p2,Double_t p3,Double_t p4,Double_t p5);

private:

  AliAnalysisTaskPi0FlowMCParamWeights(const AliAnalysisTaskPi0FlowMCParamWeights&); // not implemented
  AliAnalysisTaskPi0FlowMCParamWeights& operator=(const AliAnalysisTaskPi0FlowMCParamWeights&); // not implemented
  
protected:
  
  Double_t PrimaryWeight(Int_t primary);
  Double_t PrimaryParticleWeight(TParticle * particle);
  
  ClassDef(AliAnalysisTaskPi0FlowMCParamWeights, 1); // PHOS analysis task
};

#endif // ALIANALYSISTASKPI0FLOWMCPARAMWEIGHTS_H
