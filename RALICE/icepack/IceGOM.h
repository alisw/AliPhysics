#ifndef IceGOM_h
#define IceGOM_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "AliDevice.h"

class IceGOM : public AliDevice
{
 public:
  IceGOM();                                          // Default constructor
  virtual ~IceGOM();                                 // Default destructor
  IceGOM(const IceGOM& m);                           // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer

 ClassDef(IceGOM,1) // Signal (Hit) handling of a generic IceCube Optical Module (GOM).
};
#endif
