#ifndef ALIVZERORECONSTRUCTOR_H
#define ALIVZERORECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///                                                                          //
/// class for VZERO reconstruction                                           //
///                                                                          //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"


class AliVZEROReconstructor: public AliReconstructor {
public:
  AliVZEROReconstructor(): AliReconstructor() {};
  virtual ~AliVZEROReconstructor() {};

  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;

  ClassDef(AliVZEROReconstructor, 0)   // class for the VZERO reconstruction
};

#endif
