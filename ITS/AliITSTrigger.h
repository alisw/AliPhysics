#ifndef ALIITSTRIGGER_H
#define ALIITSTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////
//
//  ITS SPD Trigger Detector Class
//
//
//
/////////////////////////////////////////////////

#include "AliTriggerDetector.h"

class AliITSTrigger : public AliTriggerDetector
{
 public:
                   AliITSTrigger();   // constructor
        virtual   ~AliITSTrigger(){}  // destructor
   virtual void    CreateInputs();
   virtual void    Trigger();

private:
         Bool_t    RequireZ10cm(Int_t iFOperChipinStave[][40][2], Int_t stave1, Int_t stave2);

  ClassDef( AliITSTrigger, 1 )  // ITS SPD Trigger Detector class
};

#endif
