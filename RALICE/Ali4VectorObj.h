#ifndef ALI4VECTOROBJ_H
#define ALI4VECTOROBJ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include "TObject.h"

#include "Ali4Vector.h"
 
class Ali4VectorObj : public TObject,public Ali4Vector
{
 public:
  Ali4VectorObj();                       // Default constructor
  Ali4VectorObj(Ali4Vector& q);          // Constructor
  virtual ~Ali4VectorObj();              // Destructor
  Ali4VectorObj(Ali4VectorObj& q);       // Copy constructor
  void Load(Ali4Vector& q);              // Load all attributes of input Ali4Vector

 ClassDef(Ali4VectorObj,2) // Handling of Lorentz 4-vectors in various reference frames.
};
#endif
