#ifndef ALIANALYSISTASKJETS_H
#define ALIANALYSISTASKJETS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"
class AliJetFinder;
class AliESD;
class TChain;


class AliAnalysisTaskJets : public AliAnalysisTask
{
 public: 
  AliAnalysisTaskJets(const char* name);
  virtual ~AliAnalysisTaskJets() {;}
  // Implementation of interface methods
  virtual void Init(Option_t *option);
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option); 
 private:
  AliJetFinder* fJetFinder;
  TChain*       fChain;
  AliESD*       fESD;
  
  ClassDef(AliAnalysisTaskJets, 1); // Analysis task for standard jet analysis
};
 
#endif
