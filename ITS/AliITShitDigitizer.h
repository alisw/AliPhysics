#ifndef ALIITSDIGITZER_H
#define ALIITSDIGITZER_H
/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

class TObjArray;

#include <TClonesArray.h> // function of this class used in inline functions.

class AliRunDigitizer;

#include "AliDigitizer.h" // Base class from which this one is derived
#include "AliITS.h"   // ITS class functions used in inline functions.
class AliITShit;
class AliITSmodule;



class AliITSDigitizer : public AliDigitizer{
 public:
    AliITSDigitizer();
    AliITSDigitizer(AliRunDigitizer *manager);
    virtual ~AliITSDigitizer();
    // Standard routines.
    virtual Bool_t Init();
    virtual void Exec(Option_t* opt=0);
 private:
    // Routines used internaly
    TClonesArray* GetHits(){return fITS->Hits();}
    AliITShit* GetHit(Int_t h){return (AliITShit*)(GetHits()->UncheckedAt(h));}
    TObjArray* GetModules(){return fITS->GetModules();}
    AliITSmodule* GetModule(Int_t i){return fITS->GetModule(i);}
    AliRunDigitizer* GetManager(){return fManager;}
 private:
    AliITS *fITS;  //! local pointer to ITS

    ClassDef(AliITSDigitizer,1) // Task to Digitize ITS from summable hits.
};
#endif
