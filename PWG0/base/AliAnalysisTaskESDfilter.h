#ifndef ALIANALYSISTASKESDFILTER_H
#define ALIANALYSISTASKESDFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TList.h> 
#include "AliAnalysisTask.h"

class AliESDEvent;
class TChain;
class AliAODEvent;
class AliAnalysisFilter;

class AliAnalysisTaskESDfilter : public AliAnalysisTask
{
 public:
    AliAnalysisTaskESDfilter();
    AliAnalysisTaskESDfilter(const char* name);
    virtual ~AliAnalysisTaskESDfilter() {;}
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option = "");
    virtual void CreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t *option);
    // Setters
    virtual void SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
    virtual void SetKinkFilter (AliAnalysisFilter*  KinkF) {fKinkFilter  =  KinkF;}
    virtual void SetV0Filter   (AliAnalysisFilter*    V0F) {fV0Filter    =    V0F;}
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    
 private:
    Int_t              fDebug;       //  Debug flag
    TTree*             fTree;        //! chained files
    AliESDEvent*       fESD;         //! ESD
    AliAODEvent*       fAOD;         //! AOD event 
    TTree*             fTreeA;       //! AOD tree
    AliAnalysisFilter* fTrackFilter; //  Track Filter
    AliAnalysisFilter* fKinkFilter;  //  Kink  Filter
    AliAnalysisFilter* fV0Filter;    //  V0    Filter    
    ClassDef(AliAnalysisTaskESDfilter, 1); // Analysis task for standard ESD filtering
};
 
#endif
