#ifndef ALIANALYSISTASKSE_H
#define ALIANALYSISTASKSE_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliAnalysisTask.h"
class AliVEvent;
class AliAODEvent;
class AliMCEvent;
class TTree;



class AliAnalysisTaskSE : public AliAnalysisTask
{
 public:
    AliAnalysisTaskSE();
    AliAnalysisTaskSE(const char* name);
    virtual ~AliAnalysisTaskSE() {;}
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option = "");
    virtual void CreateOutputObjects();
    virtual void Exec(Option_t* option);
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    virtual void Init() {;}
    // To be implemented by user
    virtual void UserCreateOutputObjects()  {;}
    virtual void UserExec(Option_t* option) {;}
    
    // Getters
    virtual AliVEvent*   InputEvent()  {return fInputEvent;}
    virtual AliAODEvent* AODEvent()    {return fOutputAOD;}
    virtual TTree*       OutputTree()  {return fTreeA;}
    virtual AliMCEvent*  MCEvent()     {return fMCEvent;}
 protected:
    Int_t         fDebug;        //  Debug flag
    AliVEvent*    fInputEvent;   //! VEvent Input
    AliAODEvent*  fOutputAOD;    //! AOD out 
    AliMCEvent*   fMCEvent;      //! MC
    TTree*        fTreeA;        //  AOD output Tree
    ClassDef(AliAnalysisTaskSE, 1); // Analysis task for standard jet analysis
};
 
#endif
