#ifndef ALIRICHRECONSTRUCTOR_H
#define ALIRICHRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"

class AliRICH;


class AliRICHReconstructor: public AliReconstructor {
public:
  AliRICHReconstructor(): AliReconstructor() {};
  virtual ~AliRICHReconstructor() {};

  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;

private:
  AliRICH*             GetRICH(AliRunLoader* runLoader) const;

  ClassDef(AliRICHReconstructor, 0)   // class for the RICH reconstruction
};

#endif
