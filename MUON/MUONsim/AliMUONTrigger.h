/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALIMUONTRIGGER_H
#define ALIMUONTRIGGER_H

// $Id$

/// \ingroup sim
/// \class AliMUONTrigger
/// \brief MUON trigger detector class
//  Author: E. Lopez Torres


#include "AliTriggerDetector.h"

class AliMUONVTriggerStore;

class AliMUONTrigger : public AliTriggerDetector
{
 public:
   AliMUONTrigger();  // constructor
   virtual ~AliMUONTrigger();  // destructor
   virtual void    CreateInputs();
   virtual void    Trigger();

private:
   /// Not implemented
   AliMUONTrigger(const AliMUONTrigger&);
   /// Not implemented
   AliMUONTrigger& operator=(const AliMUONTrigger&);
   
   AliMUONVTriggerStore* fTriggerStore; //!<! trigger store
   
  ClassDef(AliMUONTrigger,1)  // MUON Trigger Detector class
};
#endif








