#ifndef ALIVERTEXGENERATOR_H
#define ALIVERTEXGENERATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <TVector3.h>


class AliVertexGenerator: public TObject {
 public:
  virtual TVector3 GetVertex() = 0;

  ClassDef(AliVertexGenerator, 1)    // Base class for vertex generators
};

#endif














