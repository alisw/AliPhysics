/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALIMUONTRIGGER_H
#define ALIMUONTRIGGER_H

/// \ingroup sim
/// \class AliMUONTrigger
/// \brief MUON trigger detector class

#include "AliTriggerDetector.h"

class AliMUONTrigger : public AliTriggerDetector
{
 public:
   AliMUONTrigger();  // constructor
  virtual ~AliMUONTrigger(){}  // destructor
   virtual void    CreateInputs();
   virtual void    Trigger();

  ClassDef(AliMUONTrigger,1)  // MUON Trigger Detector class
};
#endif








