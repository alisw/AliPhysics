#ifndef ALIDIGITIZER_H
#define ALIDIGITIZER_H
/* Copyright(c) 1998-2001, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$*/

////////////////////////////////////////////////////////////////////////
//
//  Base Class for Detector specific Merging/Digitization   
//                  
//  Author: Jiri Chudoba (CERN)
//
////////////////////////////////////////////////////////////////////////

#include "TTask.h"

class AliRunDigitizer;

class AliDigitizer: public TTask {

 public:
    AliDigitizer();                               // default ctor - dummy
    AliDigitizer(AliRunDigitizer *manager);       // ctor to be used          
    virtual ~AliDigitizer();
    virtual Bool_t Init() {return kTRUE;}
//    virtual void Digitize() = 0;

 protected:
    AliRunDigitizer *fManager;
    
    ClassDef(AliDigitizer,1)
};

#endif // ALIDIGITIZER_H
