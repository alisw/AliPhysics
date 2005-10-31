#ifndef IceDwalk_h
#define IceDwalk_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"

#include "AliJob.h"
#include "AliSample.h"
#include "IceEvent.h"
#include "IceGOM.h"

class IceDwalk : public TTask
{
 public :
  IceDwalk(const char* name="",const char* title=""); // Constructor
  virtual ~IceDwalk();                                // Destructor
  virtual void Exec(Option_t* opt);                   // Direct walk reconstruction
  void SetDmin(Float_t d);      // Set minimum hit distance to form a track element
  void SetDtmarg(Int_t dt);     // Set maximum hit time difference margin for track elements
  void SetTangsep(Float_t ang); // Set angular separation within which track candidates are clustered in a jet
  void SetJangsep(Float_t ang); // Set angular separation within which jets are merged into 1 single track

 protected :
  Float_t fDmin;    // Minimum hit distance (in m) to form a track element 
  Int_t fDtmarg;    // Maximum hit time difference margin (in ns) for track elements
  Float_t fTangsep; // Angular separation (in deg) within which track candidates are clustered in a jet
  Float_t fJangsep; // Angular separation (in deg) within which jets are merged into 1 single track

 ClassDef(IceDwalk,1) // TTask derived class to perform direct walk reconstruction
};
#endif
