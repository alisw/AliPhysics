#ifndef ALIATTRIBOBJ_H
#define ALIATTRIBOBJ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include "TObject.h"

#include "AliAttrib.h"

class AliAttribObj : public TObject,public AliAttrib
{
 public:
  AliAttribObj();                                    // Default constructor
  AliAttribObj(AliAttrib& a);                        // Constructor
  virtual ~AliAttribObj();                           // Destructor
  AliAttribObj(const AliAttribObj& a);               // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer

 ClassDef(AliAttribObj,4) // Generic handling of detector signal (calibration) attributes.
};
#endif
