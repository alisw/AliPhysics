#ifndef IceTDOM_h
#define IceTDOM_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "IceDOM.h"

class IceTDOM : public IceDOM
{
 public:
  IceTDOM();                                         // Default constructor
  virtual ~IceTDOM();                                // Default destructor
  IceTDOM(const IceTDOM& m);                         // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer

 ClassDef(IceTDOM,1) // Signal (Hit) handling of an IceTop Digital Optical Module (TDOM).
};
#endif
