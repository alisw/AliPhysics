#ifndef ALIITSFDIGITZER_H
#define ALIITSFDIGITZER_H
/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

class TObjArray;
class TTree;

#include <TClonesArray.h> // function of this class used in inline functions.

class AliRunDigitizer;

#include "AliDigitizer.h" // Base class from which this one is derived
#include "AliITS.h"   // ITS class functions used in inline functions.
class AliITSmodule;

class AliITSFDigitizer : public AliDigitizer{
 public:
    AliITSFDigitizer();
    AliITSFDigitizer(AliRunDigitizer *manager);
    virtual ~AliITSFDigitizer();
    // Standard routines.
    virtual Bool_t Init();
    // Perform SDigits to Digits, with or without merging, depending on the
    // number of files.
    virtual void Exec(Option_t* opt=0);
 private:
    // Routines used internaly
    // Returns a pointer to the TObjecArray of Modules.
    TObjArray* GetModules(){return fITS->GetModules();}
    // Returns a pointer to a  specific module.
    AliITSmodule* GetModule(Int_t i){return fITS->GetModule(i);}
 private:
    AliITS *fITS;      //! local pointer to ITS
    Bool_t  fInit;     //! flag to indicate Initilization when well.


    ClassDef(AliITSFDigitizer,1) // Task to Digitize ITS from summable hits.
};
#endif
