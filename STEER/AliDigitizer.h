#ifndef ALIDIGITIZER_H
#define ALIDIGITIZER_H
/* Copyright(c) 1998-2001, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$*/

////////////////////////////////////////////////////////////////////////
//
//  Abstract Base Class for Detector specific Merging/Digitization   
//                  
//  Author: Jiri Chudoba (CERN)
//
////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
class AliRunDigitizer;

class AliDigitizer: public TNamed {

 public:
    AliDigitizer();                               // default ctor - dummy
    AliDigitizer(AliRunDigitizer *manager);       // ctor to be used          
    virtual ~AliDigitizer();
    virtual Bool_t Init() = 0;
    virtual void Digitize() = 0;
    Int_t GetDebug() {return fDebug;}             // get debug level
    void SetDebug(Int_t level){fDebug = level;}   // set debug level
    
private:

    Int_t fDebug;                                 // debug level

    ClassDef(AliDigitizer,1)
};

#endif // ALIDIGITIZER_H
