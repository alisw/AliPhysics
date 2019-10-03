/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

#ifndef ALIPHYSICSSELECTIONTASK_H
#define ALIPHYSICSSELECTIONTASK_H

#include "AliAnalysisTaskSE.h"

class AliPhysicsSelection;

class AliPhysicsSelectionTask : public AliAnalysisTaskSE {
  public:
    AliPhysicsSelectionTask();
    AliPhysicsSelectionTask(const char* opt);

    virtual ~AliPhysicsSelectionTask();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t*);
    virtual void   FinishTaskOutput();
    virtual void   Terminate(Option_t*);

    void SetOption(const char* opt) { fOption = opt; }
    
    void SetPhysicsSelection(AliPhysicsSelection* physicsSelection) { fPhysicsSelection = physicsSelection; }
    AliPhysicsSelection* GetPhysicsSelection() const { return fPhysicsSelection; }
    void   SetUseSpecialOutput(Bool_t v=kTRUE)         {fUseSpecialOutput = v;}
    Bool_t GetUseSpecialOutput()               const { return fUseSpecialOutput; }
    static AliPhysicsSelectionTask* AddTaskPhysicsSelection ( const Bool_t mCAnalysisFlag = kFALSE, const Bool_t applyPileupCuts = kFALSE, const UInt_t deprecatedFlag2 = 0, const Bool_t useSpecialOutput=kFALSE);

 protected:
    TList* fOutput;                  //! list send on output slot 1
    TString fOption;                 // option string  
    Bool_t  fUseSpecialOutput;       // do we use special output instead of merging?
    AliPhysicsSelection* fPhysicsSelection; // event selection class
 private:
    AliPhysicsSelectionTask(const AliPhysicsSelectionTask&);
    AliPhysicsSelectionTask& operator=(const AliPhysicsSelectionTask&);

  ClassDef(AliPhysicsSelectionTask, 2);
};

#endif
