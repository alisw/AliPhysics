#ifndef ALIDIGITIZER_H
#define ALIDIGITIZER_H
/* Copyright(c) 1998-2001, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$*/

////////////////////////////////////////////////////////////////////////
//
//  Base Class for Detector specific Merging/Digitization   
//  Detector specific digitization classes derive from this
//  Author: Jiri Chudoba (CERN)
//
////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "AliDigitizationInput.h"

class AliDigitizer: public TNamed {

 public:
// ctor with name and title
    AliDigitizer(const Text_t* name="AliDigitizer",
                const Text_t* title="AliDigitizer");
// ctor to be used with name and title
    AliDigitizer(AliDigitizationInput *manager,
                 const Text_t* name="AliDigitizer",
                 const Text_t* title="AliDigitizer");
// Copy ctor needed because there is a pointer
    AliDigitizer(const AliDigitizer &dig);
    AliDigitizer& operator=(const AliDigitizer &dig)
      {dig.Copy(*this);return *this;}
      
    virtual ~AliDigitizer();
    virtual Bool_t Init() {return kTRUE;}
    virtual void Digitize(Option_t* option) = 0;
    Bool_t GetRegionOfInterest() const {return fDigInput ? fDigInput->GetRegionOfInterest() : kFALSE;}

protected:
    Int_t GetNInputStreams() const;
    void Copy(TObject &dig) const;

    AliDigitizationInput *fDigInput;   //! Pointer to the Digitizer input
    
    ClassDef(AliDigitizer,3) // Base class for detector digitizers
};

#endif // ALIDIGITIZER_H

