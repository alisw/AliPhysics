#ifndef AliAnalysisTaskParticleCorrelation_H
#define AliAnalysisTaskParticleCorrelation_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTaskSE.h"
class AliAnaPartCorrMaker;
class AliESDEvent;
class AliAODEvent;
class TList;

class AliAnalysisTaskParticleCorrelation : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskParticleCorrelation();
    AliAnalysisTaskParticleCorrelation(const char* name);
    virtual ~AliAnalysisTaskParticleCorrelation() ;// virtual dtor
 
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    void SetConfigFileName(TString name ) {fConfigName = name ; }
    TString GetConfigFileName() const {return fConfigName ; }

 private:
    AliAnalysisTaskParticleCorrelation(const AliAnalysisTaskParticleCorrelation&); // Not implemented
    AliAnalysisTaskParticleCorrelation& operator=(const AliAnalysisTaskParticleCorrelation&); // Not implemented

    AliAnaPartCorrMaker* fAna; //  Pointer to the jet finder 
    TList * fOutputContainer ; //! Histogram container
    TClonesArray * fAODBranch; //! AOD branch
    TString fConfigName ; //Configuration file name

    ClassDef(AliAnalysisTaskParticleCorrelation, 1); // Analysis task for standard gamma correlation analysis
};
 
#endif //AliAnalysisTaskParticleCorrelation_H
