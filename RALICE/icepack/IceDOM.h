#ifndef IceDOM_h
#define IceDOM_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "IceGOM.h"

class IceDOM : public IceGOM
{
 public:
  IceDOM();                                          // Default constructor
  virtual ~IceDOM();                                 // Default destructor
  IceDOM(const IceDOM& m);                           // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer

 ClassDef(IceDOM,1) // Signal (Hit) handling of a generic IceCube Digital Optical Module (DOM).
};
#endif
