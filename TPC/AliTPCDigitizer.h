#ifndef ALITPCDIGITIZER_H
#define ALITPCDIGITIZER_H
/* Copyright(c) 1998-2001, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliDigitizer.h"

class AliRunDigitizer;

class AliTPCDigitizer : public AliDigitizer {
 public:    
    AliTPCDigitizer();
    AliTPCDigitizer(AliRunDigitizer * manager);
    virtual ~AliTPCDigitizer();
    // Initialize merging and digitization
    virtual Bool_t Init();
    // Do the main work
    virtual void Exec(Option_t* option=0);    
    Int_t GetDebug() const {return fDebug;}       // get debug level
    void SetDebug(Int_t level){fDebug = level;}   // set debug level        
 private:    
    Int_t fDebug;
 private:
    ClassDef(AliTPCDigitizer,1)  // MUON merging/digitization
};    
#endif

