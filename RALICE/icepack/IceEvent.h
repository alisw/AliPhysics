#ifndef IceEvent_h
#define IceEvent_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "AliEvent.h"
#include "IceGOM.h"

class IceEvent : public AliEvent
{
 public:
  IceEvent();                                        // Default constructor
  virtual ~IceEvent();                               // Default destructor
  IceEvent(const IceEvent& evt);                     // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer
  Int_t GetNstrings(TString classname);              // Number of fired strings for good mods. of the specified class for this event 
  Int_t GetNstrings(AliTrack& t,TString classname);  // Number of fired strings for good mods. of the specified class associated to this track
  Int_t GetNstrings(AliJet& j,TString classname);    // Number of fired strings for good mods. of the specified class associated to this jet

 protected:
  TArrayI* fStrings; //! Temp. array to hold the string ids of fired modules

 ClassDef(IceEvent,2) // Handling of IceCube event data.
};
#endif
