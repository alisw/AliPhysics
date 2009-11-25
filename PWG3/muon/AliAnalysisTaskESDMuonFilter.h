#ifndef ALIANALYSISTASKESDMUONFILTER_H
#define ALIANALYSISTASKESDMUONFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TList.h> 
#include "AliAnalysisTaskSE.h"

class AliAnalysisFilter;
class AliStack;

class AliAnalysisTaskESDMuonFilter : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskESDMuonFilter();
    AliAnalysisTaskESDMuonFilter(const char* name);
    virtual ~AliAnalysisTaskESDMuonFilter() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    virtual void ConvertESDtoAOD();

    // Setters
    virtual void SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
    void SetWriteMuonAOD(Bool_t enableMuonAOD){fEnableMuonAOD = enableMuonAOD;}

 private:
    AliAnalysisTaskESDMuonFilter(const AliAnalysisTaskESDMuonFilter&);
    AliAnalysisTaskESDMuonFilter& operator=(const AliAnalysisTaskESDMuonFilter&);
    void PrintMCInfo(AliStack *pStack,Int_t label); // for debugging
    AliAnalysisFilter* fTrackFilter; //  Track Filter
    Bool_t fEnableMuonAOD; // flag for enabling Muon AOD production
    ClassDef(AliAnalysisTaskESDMuonFilter, 1); // Analysis task for standard ESD filtering

};
 
#endif
