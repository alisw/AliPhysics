#ifndef ALI3VECTOROBJ_H
#define ALI3VECTOROBJ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"

#include "Ali3Vector.h"
 
class Ali3VectorObj : public TObject,public Ali3Vector
{
 public:
  Ali3VectorObj();               // Default constructor
  Ali3VectorObj(Ali3Vector& q);  // Constructor
  ~Ali3VectorObj();              // Destructor
  void Load(Ali3Vector& q);      // Load all attributes of input Ali3Vector

 ClassDef(Ali3VectorObj,1) // Handling of 3-vectors in various reference frames.
};
#endif
