#ifndef ALIANALYSISTASKGAMMA_H
#define ALIANALYSISTASKGAMMA_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTaskSE.h"
class AliAnaGamma;
class AliESDEvent;
class AliAODEvent;
class TList;

class AliAnalysisTaskGamma : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskGamma();
    AliAnalysisTaskGamma(const char* name);
    virtual ~AliAnalysisTaskGamma() ;// virtual dtor
 
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    void SetConfigFileName(TString name ) {fConfigName = name ; }
    TString GetConfigFileName() const {return fConfigName ; }

 private:
    AliAnalysisTaskGamma(const AliAnalysisTaskGamma&); // Not implemented
    AliAnalysisTaskGamma& operator=(const AliAnalysisTaskGamma&); // Not implemented

    AliAnaGamma* fAna; //  Pointer to the jet finder 
    TList * fOutputContainer ; // Histogram container
    TString fConfigName ; //Configuration file name

    ClassDef(AliAnalysisTaskGamma, 2); // Analysis task for standard gamma correlation analysis
};
 
#endif //ALIANALYSISTASKGAMMA_H
