#ifndef IceEvent_h
#define IceEvent_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "AliEvent.h"

class IceEvent : public AliEvent
{
 public:
  IceEvent();                                        // Default constructor
  virtual ~IceEvent();                               // Default destructor
  IceEvent(const IceEvent& evt);                     // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer

 ClassDef(IceEvent,1) // Handling of IceCube event data.
};
#endif
