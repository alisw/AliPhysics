#ifndef EvtAna_h
#define EvtAna_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"

#include "AliJob.h"
#include "IceEvent.h"
#include "IceAOM.h"

class EvtAna : public TTask
{
 public :
  EvtAna(const char* name="",const char* title=""); // Constructor
  virtual ~EvtAna();                                // Destructor
  virtual void Exec(Option_t* opt);                 // The actions to be performed 

 ClassDef(EvtAna,1) // Just a class to test the IceTask concept
};
#endif
