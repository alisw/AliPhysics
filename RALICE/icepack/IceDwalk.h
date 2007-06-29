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
  IceDwalk(const char* name="IceDwalk",const char* title="Direct walk reconstruction"); // Constructor
  virtual ~IceDwalk();                                // Destructor
  virtual void Exec(Option_t* opt);                   // Direct walk reconstruction
  void SetDmin(Float_t d);       // Set minimum hit distance to form a track element
  void SetDtmarg(Int_t dt);      // Set maximum hit time difference margin for track elements
  void SetMaxDhit(Float_t d);    // Set maximum distance (in scat. length) for hit association 
  void SetTangmax(Float_t ang);  // Set max. angular separation for track candidate clustering into jets
  void SetTdistmax(Float_t d,Int_t invol=1);// Set maximum track distance for track candidate clustering
  void SetJangmax(Float_t ang,Int_t iter=1);// Set max. angular separation for jet merging into 1 single track
  void SetJdistmax(Float_t d,Int_t invol=1);// Set maximum jet distance for jet merging
  void SetMaxModA(Int_t nmax);   // Set max. number of good fired Amanda modules for events to be processed
  void SetMinModA(Int_t nmin);   // Set min. number of good fired Amanda modules for events to be processed
  void SetMaxHitsA(Int_t nmax);  // Set max. number of good hits per Amanda module to be processed
  void SetVgroupUsage(Int_t flag); // (De)activate usage of distinct phase and group velocities
  void SetTrackName(TString s);  // Set (alternative) name for the produced first guess tracks
  void SetCharge(Float_t charge);// Set user defined charge for the produced first guess tracks

 protected :
  IceEvent* fEvt;    // Pointer to the event structure
  Float_t fDmin;     // Minimum hit distance (in m) to form a track element 
  Int_t fDtmarg;     // Maximum hit time difference margin (in ns) for track elements
  Float_t fMaxdhit;  // Maximum distance (in scat. length) for hit association
  Float_t fTangmax;  // Angular separation (in deg) within which track candidates are clustered in a jet
  Float_t fTdistmax; // Maximum track distance (in m) for track candidate clustering
  Int_t fTinvol;     // Flag to denote maximum track distance testing inside/outside detector volume
  Float_t fJangmax;  // Angular separation (in deg) within which jets are merged into 1 single track
  Int_t fJiterate;   // Flag to indicate iteration in the jet merging process
  Float_t fJdistmax; // Maximum jet distance (in m) for jet merging
  Int_t fJinvol;     // Flag to denote maximum jet distance testing inside/outside detector volume
  Int_t fMaxmodA;    // The max. number of good fired Amanda modules for events to get processed
  Int_t fMinmodA;    // The min. number of good fired Amanda modules for events to get processed
  Int_t fMaxhitsA;   // The maximum number of good hits per Amanda module to be processed
  Int_t fVgroup;     // Flag to indicate usage of distinct phase and group velocities
  TString fTrackname;// The name identifier for the produced first guess tracks
  Float_t fCharge;   // User defined charge of the produced first guess tracks

  virtual void AssociateHits(TObjArray& tes,TObjArray& hits,Float_t& qmax,Int_t& nahmax); // Hit association
  virtual void SelectQvalue(TObjArray& tes,Float_t qmax);                  // TC selection via Q-value
  virtual void ClusterTracks(TObjArray& tes,TObjArray& jets,Float_t qmax); // Track clustering  
  virtual void MergeJets(TObjArray& jets);                                 // Jet Merging
  virtual void StoreTracks(TObjArray& jets);                               // Final track storage

 ClassDef(IceDwalk,8) // TTask derived class to perform (improved) direct walk reconstruction
};
#endif
