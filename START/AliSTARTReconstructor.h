#ifndef ALISTARTRECONSTRUCTOR_H
#define ALISTARTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"

class AliSTARTReconstructor: public AliReconstructor {
public:
  AliSTARTReconstructor(): AliReconstructor() {};
  virtual ~AliSTARTReconstructor() {};
  virtual void         Reconstruct(AliRunLoader* /*runLoader*/) const;
   virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;
 
  ClassDef(AliSTARTReconstructor, 0)   // class for the START reconstruction
};

#endif
