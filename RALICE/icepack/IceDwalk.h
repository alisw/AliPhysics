#ifndef IceDwalk_h
#define IceDwalk_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"
#include "TArrayI.h"

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
  void SetDmin(Float_t d);       // Set minimum hit distance to form a track element
  void SetDtmarg(Int_t dt);      // Set maximum hit time difference margin for track elements
  void SetTangmax(Float_t ang);  // Set max. angular separation for track candidate clustering into jets
  void SetRtangmax(Float_t ang); // Set max. ang. separation for relative r0 direction for track candidate clustering
  void SetRtdmax(Float_t d);     // Set maximum r0 distance for track candidate clustering
  void SetJangmax(Float_t ang);  // Set max. angular separation for jet merging into 1 single track
  void SetRjangmax(Float_t ang); // Set max. angular separation for relative r0 direction for jet merging
  void SetRjdmax(Float_t d);     // Set maximum r0 distance for jet merging
  void SetMaxModA(Int_t nmax);   // Set max. number of good fired Amanda modules for events to be processed
  void SetMinModA(Int_t nmin);   // Set min. number of good fired Amanda modules for events to be processed
  void SetMaxHitsA(Int_t nmax);  // Set max. number of good hits per Amanda module to be processed
  void SetTrackName(TString s);  // Set (alternative) name for the produced first guess tracks
  void SetCharge(Float_t charge);// Set user defined charge for the produced first guess tracks

 protected :
  Float_t fDmin;     // Minimum hit distance (in m) to form a track element 
  Int_t fDtmarg;     // Maximum hit time difference margin (in ns) for track elements
  Float_t fTangmax;  // Angular separation (in deg) within which track candidates are clustered in a jet
  Float_t fRtangmax; // Relative r0 angular separation (in deg) for track candidate clustering
  Float_t fRtdmax;   // Maximum r0 distance (in m) for track candidate clustering
  Float_t fJangmax;  // Angular separation (in deg) within which jets are merged into 1 single track
  Float_t fRjangmax; // Relative r0 angular separation (in deg) for jet clustering
  Float_t fRjdmax;   // Maximum r0 distance (in m) for jet merging
  Int_t fMaxmodA;    // The max. number of good fired Amanda modules for events to get processed
  Int_t fMinmodA;    // The min. number of good fired Amanda modules for events to get processed
  Int_t fMaxhitsA;   // The maximum number of good hits per Amanda module to be processed
  TString fTrackname;// The name identifier for the produced first guess tracks
  Float_t fCharge;   // User defined charge of the produced first guess tracks

 ClassDef(IceDwalk,6) // TTask derived class to perform direct walk reconstruction
};
#endif
