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
    AliAnalysisTaskJets();
    AliAnalysisTaskJets(const char* name);
    virtual ~AliAnalysisTaskJets() {;}
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option = "");
    virtual void Init();
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    
 private:
    Int_t         fDebug;     //  Debug flag
    AliJetFinder* fJetFinder; //  Pointer to the jet finder 
    TChain*       fChain;     //! chained files
    AliESD*       fESD;       //! ESD
    
    ClassDef(AliAnalysisTaskJets, 1); // Analysis task for standard jet analysis
};
 
#endif
