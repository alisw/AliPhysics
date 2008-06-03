#ifndef ALIANALYSISTASKJETSCDF_H
#define ALIANALYSISTASKJETSCDF_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTaskSE.h"

class AliCdfJetFinder;
class TList;

class AliAnalysisTaskJetsCDF : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskJetsCDF();
    AliAnalysisTaskJetsCDF(const char* name);
    virtual ~AliAnalysisTaskJetsCDF() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
 private:
  AliAnalysisTaskJetsCDF(const AliAnalysisTaskJetsCDF &det);
  AliAnalysisTaskJetsCDF &operator=(const AliAnalysisTaskJetsCDF &det);
    
 private:
    AliCdfJetFinder* fJetFinder;    //  Pointer to the jet finder 
    TList*           fListOfHistos; //  Output list of histograms

    ClassDef(AliAnalysisTaskJetsCDF, 1); // Analysis task for CDF jet analysis
};
 
#endif
