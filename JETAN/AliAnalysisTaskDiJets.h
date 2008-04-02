#ifndef ALIANALYSISTASKDIJETS_H
#define ALIANALYSISTASKDIJETS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTaskSE.h"
class AliAnalysisTaskDiJets : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskDiJets();
    AliAnalysisTaskDiJets(const char* name);
    virtual ~AliAnalysisTaskDiJets() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
 private:
  AliAnalysisTaskDiJets(const AliAnalysisTaskDiJets &det);
  AliAnalysisTaskDiJets &operator=(const AliAnalysisTaskDiJets &det);
    
 private:
  TClonesArray* fDiJets;    // Array of dijets
  TClonesArray* fDiJetsIn;  // Array of dijets
  ClassDef(AliAnalysisTaskDiJets, 1); // Analysis task for standard jet analysis
};
 
#endif
