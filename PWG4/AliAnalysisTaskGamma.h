#ifndef ALIANALYSISTASKGAMMA_H
#define ALIANALYSISTASKGAMMA_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"
class AliAnaGamma;
class AliESDEvent;
class AliAODEvent;
class TChain;
class TList;

class AliAnalysisTaskGamma : public AliAnalysisTask
{
 public:
    AliAnalysisTaskGamma();
    AliAnalysisTaskGamma(const char* name);
    virtual ~AliAnalysisTaskGamma() ;// virtual dtor
 
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option = "");
    virtual void CreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t *option);

    void SetConfigFileName(TString name ) {fConfigName = name ; }
    TString GetConfigFileName() const {return fConfigName ; }

 private:

    AliAnaGamma* fAna; //  Pointer to the jet finder 
    TChain*       fChain;     //! chained files
    AliESDEvent*       fESD;       //! ESD
    AliAODEvent*       fAOD;       //! AOD
    TTree*        fTreeG;     //  tree of prompt gamma, does nothing for the moment 
    TList * fOutputContainer ; // Histogram container
    TString fConfigName ; //Configuration file name

    ClassDef(AliAnalysisTaskGamma, 1); // Analysis task for standard gamma correlation analysis
};
 
#endif //ALIANALYSISTASKGAMMA_H
