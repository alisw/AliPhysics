/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

#ifndef ALICOLLISIONNORMALIZATIONTASK_H
#define ALICOLLISIONNORMALIZATIONTASK_H

#include "AliAnalysisTaskSE.h"

class AliCollisionNormalization;

class AliCollisionNormalizationTask : public AliAnalysisTaskSE {
  public:
    AliCollisionNormalizationTask();
    AliCollisionNormalizationTask(const char* name);

    virtual ~AliCollisionNormalizationTask();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t*);
    virtual void   Terminate(Option_t*);

  Bool_t IsEventInBinZero(); // returns true if the event has to be put in the bin0.
  void SetMC(Bool_t flag = kTRUE) { fIsMC = flag; }
  
  //    void SetOption(const char* opt) { fOption = opt; }
    
  //    void SetCollisionNormalization(AliCollisionNormalization* physicsSelection) { fCollisionNormalization = physicsSelection; }
  AliCollisionNormalization* GetCollisionNormalization() const { return fCollisionNormalization; }

 protected:
  TList* fOutput;                  //! list send on output slot 1
  //    TString fOption;                 // option string  
  Bool_t fIsMC;
  
  AliCollisionNormalization* fCollisionNormalization; // collision normalization class

 private:
    AliCollisionNormalizationTask(const AliCollisionNormalizationTask&);
    AliCollisionNormalizationTask& operator=(const AliCollisionNormalizationTask&);

  ClassDef(AliCollisionNormalizationTask, 1);
};

#endif
