#ifndef ALIANALYSISTASKME_H
#define ALIANALYSISTASKME_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliAnalysisTask.h"
class AliVEvent;
class AliAODEvent;
class AliInputEventHandler;
class TTree;
class AliMultiEventInputHandler;



class AliAnalysisTaskME : public AliAnalysisTask
{
 public:
    AliAnalysisTaskME();
    AliAnalysisTaskME(const char* name);
    AliAnalysisTaskME(const AliAnalysisTaskME& obj);
    AliAnalysisTaskME& operator=(const AliAnalysisTaskME& other);
    virtual ~AliAnalysisTaskME() {;}
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option = "");
    virtual void CreateOutputObjects();
    virtual void Exec(Option_t* option);
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    virtual void Init() {;}
    virtual void RequireFreshBuffer() {fFreshBufferOnly = kTRUE;}
    // To be implemented by user
    virtual void UserCreateOutputObjects()  {;}
    virtual void UserExec(Option_t* /*option*/) {;}
    // Helpers for adding branches to the AOD
   virtual void AddAODBranch(const char* cname, void* addobj);
// Getters
    virtual Int_t          DebugLevel()              {return fDebug;     }
    virtual AliVEvent*     GetEvent(Int_t iev);
    virtual AliAODEvent*   AODEvent()                {return fOutputAOD; }
    virtual TTree*         OutputTree()              {return fTreeA;     }
    virtual Long64_t       Entry()                   {return fEntry;     }
    virtual const char*    CurrentFileName();
  protected:
    Int_t                      fDebug;           //  Debug flag
    Int_t                      fEntry;           //  Current entry in the chain
    Bool_t                     fFreshBufferOnly; //  Flag for Exec call for fresh buffer only
    AliMultiEventInputHandler* fInputHandler;    //! Input Handler
    AliAODEvent*               fOutputAOD;       //! AOD out 
    TTree*                     fTreeA;           //  AOD output Tree
    ClassDef(AliAnalysisTaskME, 1); // Analysis task for standard jet analysis
};
 
#endif
