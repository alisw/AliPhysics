#ifndef ALITPCRECONSTRUCTOR_H
#define ALITPCRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"

class AliTPCParam;


class AliTPCReconstructor: public AliReconstructor {
public:
  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual AliTracker*  CreateTracker(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;

private:
  AliTPCParam*         GetTPCParam(AliRunLoader* runLoader) const;

  ClassDef(AliTPCReconstructor, 0)   // class for the TPC reconstruction
};

#endif
