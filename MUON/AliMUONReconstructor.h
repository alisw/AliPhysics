#ifndef ALIMUONRECONSTRUCTOR_H
#define ALIMUONRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

#include "AliReconstructor.h"

class AliMUONReconstructor: public AliReconstructor 
{
  public:
    AliMUONReconstructor();
    virtual ~AliMUONReconstructor();

    virtual void         Reconstruct(AliRunLoader* runLoader) const;
    virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;
 
  ClassDef(AliMUONReconstructor, 0)   // class for the MUON reconstruction
};

#endif
