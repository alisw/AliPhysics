#ifndef ALITPCDIGITIZER_H
#define ALITPCDIGITIZER_H
/* Copyright(c) 1998-2001, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigitizer.h"

class AliDigitizationInput;

class AliTPCDigitizer : public AliDigitizer {
 public:    
    AliTPCDigitizer();
    AliTPCDigitizer(AliDigitizationInput * digInput);
    virtual ~AliTPCDigitizer();
    // Initialize merging and digitization
    virtual Bool_t Init();
    // Do the main work
    virtual void Digitize(Option_t* option=0);    
    Int_t GetDebug() const {return fDebug;}       // get debug level
    void SetDebug(Int_t level){fDebug = level;}   // set debug level        
 private: 
    void DigitizeFast(Option_t* option=0); //digitize - using row pointers
    void DigitizeSave(Option_t* option=0); // digitize using controlled arrays   
    Int_t fDebug;
 private:
    ClassDef(AliTPCDigitizer,2)  // MUON merging/digitization
};    
#endif

