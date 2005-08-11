#ifndef Analyse_h
#define Analyse_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"

#include "AliJob.h"
#include "IceAOM.h"

class Analyse : public TTask
{
 public :
  Analyse(const char* name="",const char* title=""); // Constructor
  virtual ~Analyse();                                // Destructor
  virtual void Exec(Option_t* opt);                  // The actions to be performed 

 ClassDef(Analyse,1) // Just a class to demonstrate the TTask concept
};
#endif
