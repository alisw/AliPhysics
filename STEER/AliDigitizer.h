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

#include "TTask.h"

class AliRunDigitizer;

class AliDigitizer: public TTask {

 public:
// ctor with name and title
    AliDigitizer(const Text_t* name="AliDigitizer",
                const Text_t* title="AliDigitizer");
// ctor to be used with name and title
    AliDigitizer(AliRunDigitizer *manager,
                 const Text_t* name="AliDigitizer",
                 const Text_t* title="AliDigitizer");
// Copy ctor needed because there is a pointer
    AliDigitizer(const AliDigitizer &dig);
    AliDigitizer& operator=(const AliDigitizer &dig)
      {dig.Copy(*this);return *this;}
      
    virtual ~AliDigitizer();
    virtual Bool_t Init() {return kTRUE;}
    void SetRegionOfInterest(Bool_t flag) {fRegionOfInterest = flag;};
//    virtual void Digitize() = 0;

protected:
    Int_t GetNInputStreams() const;
    void Copy(TObject &dig) const;

    AliRunDigitizer *fManager;   //! Pointer to the Digitizer manager
    Bool_t fRegionOfInterest;    // Flag for digitization only in region of interest
    
    ClassDef(AliDigitizer,2) // Base class for detector digitizers
};

#endif // ALIDIGITIZER_H

