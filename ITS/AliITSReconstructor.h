#ifndef ALIITSRECONSTRUCTOR_H
#define ALIITSRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"

class AliITSgeom;


class AliITSReconstructor: public AliReconstructor {
public:
  AliITSReconstructor(): AliReconstructor() {};
  virtual ~AliITSReconstructor() {};

  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual AliTracker*  CreateTracker(AliRunLoader* runLoader) const;
  virtual AliVertexer* CreateVertexer(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;

private:
  AliITSgeom*          GetITSgeom(AliRunLoader* runLoader) const;

  ClassDef(AliITSReconstructor, 0)   // class for the ITS reconstruction
};

#endif
