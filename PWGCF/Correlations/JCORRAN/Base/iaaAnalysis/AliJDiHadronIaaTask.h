#ifndef ALIJDIHADRONIAATASK_H
#define ALIJDIHADRONIAATASK_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim, J. Viinikainen, M. Vargyas
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////


#include <AliAnalysisTaskSE.h>
#include "AliJIaaAnalysis.h"

class AliJCORRANTask;
//==============================================================

using namespace std;

class AliJDiHadronIaaTask : public AliAnalysisTaskSE {

 public:
    AliJDiHadronIaaTask();
    AliJDiHadronIaaTask(const char *name,  TString inputformat);
    AliJDiHadronIaaTask(const AliJDiHadronIaaTask& ap);
    AliJDiHadronIaaTask& operator = (const AliJDiHadronIaaTask& ap);
    virtual ~AliJDiHadronIaaTask();

    // methods to fill from AliAnalysisTaskSE
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() { Init(); }
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
    virtual Bool_t UserNotify() { std::cout<<"DEBUG UserNotify"<<std::endl; return kTRUE;}

    void SetFilterTaskName(TString name){ fFilterTaskName=name; }
    void SetIaaAnalysis(AliJIaaAnalysis * ana){ fIaaAnalysis=ana; }

 private:
  
    AliJCORRANTask  * fFilterTask;
    TString           fFilterTaskName;
    AliJIaaAnalysis * fIaaAnalysis;
    TDirectory      * fOutput;

    ClassDef(AliJDiHadronIaaTask, 1);
};

#endif // AliJDiHadronIaaTask_H
