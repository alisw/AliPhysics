#ifndef ALIPOSITIONOBJ_H
#define ALIPOSITIONOBJ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"

#include "AliPosition.h"
 
class AliPositionObj : public TObject,public AliPosition
{
 public:
  AliPositionObj();              // Default constructor
  AliPositionObj(Ali3Vector& q); // Constructor
  ~AliPositionObj();             // Destructor
  void Load(Ali3Vector& q);      // Load all attributes of input AliPosition

 ClassDef(AliPositionObj,1) // Handling of positions in various reference frames.
};
#endif
