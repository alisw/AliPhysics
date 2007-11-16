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
  virtual ~IceDwalk();                                    // Destructor
  virtual void Exec(Option_t* opt);                       // Direct walk reconstruction
  void SetDmin(Float_t d,TString s="A");                  // Set minimum hit distance to form a track element
  void SetDtmarg(Int_t dt,TString s="A");                 // Set maximum hit time difference margin for track elements
  void SetMaxDhit(Float_t d,TString s="A");               // Set maximum distance (in scat. length) for hit association 
  void SetTangmax(Float_t ang,TString s="A");             // Set max. angular separation for track candidate clustering into jets
  void SetTdistmax(Float_t d,TString s="A",Int_t invol=1);// Set maximum track distance for track candidate clustering
  void SetJangmax(Float_t ang,TString s="A",Int_t iter=1);// Set max. angular separation for jet merging into 1 single track
  void SetJdistmax(Float_t d,TString s="A",Int_t invol=1);// Set maximum jet distance for jet merging
  void SetMaxMod(Int_t nmax,TString s="A");               // Set max. number of good fired (D)OMs for events to be processed
  void SetMinMod(Int_t nmin,TString s="A");               // Set min. number of good fired (D)OMs for events to be processed
  void SetMaxHits(Int_t nmax,TString s="A");              // Set max. number of good hits per (D)OM to be processed
  void SetVgroupUsage(Int_t flag,TString s="A");          // (De)activate usage of distinct phase and group velocities
  void SetTrackName(TString s);                           // Set (alternative) name for the produced first guess tracks
  void SetCharge(Float_t charge);                         // Set user defined charge for the produced first guess tracks

 protected :
  IceEvent* fEvt;     // Pointer to the event structure
  Float_t fDminA;     // Minimum Amanda OM hit distance (in m) to form a track element 
  Float_t fDminI;     // Minimum InIce DOM hit distance (in m) to form a track element 
  Int_t fDtmargA;     // Maximum Amanda OM hit time difference margin (in ns) for track elements
  Int_t fDtmargI;     // Maximum InIce DOM hit time difference margin (in ns) for track elements
  Float_t fMaxdhitA;  // Maximum Amanda OM hit distance (in scat. length) for hit association
  Float_t fMaxdhitI;  // Maximum InIce DOM hit distance (in scat. length) for hit association
  Float_t fTangmaxA;  // Amanda angular separation (in deg) within which track candidates are clustered in a jet
  Float_t fTangmaxI;  // InIce angular separation (in deg) within which track candidates are clustered in a jet
  Float_t fTdistmaxA; // Maximum Amanda track distance (in m) for track candidate clustering
  Float_t fTdistmaxI; // Maximum InIce track distance (in m) for track candidate clustering
  Int_t fTinvolA;     // Amanda flag to denote maximum track distance testing inside/outside detector volume
  Int_t fTinvolI;     // InIce flag to denote maximum track distance testing inside/outside detector volume
  Float_t fJangmaxA;  // Amanda angular separation (in deg) within which jets are merged into 1 single track
  Float_t fJangmaxI;  // InIce angular separation (in deg) within which jets are merged into 1 single track
  Int_t fJiterateA;   // Amanda flag to indicate iteration in the jet merging process
  Int_t fJiterateI;   // InIce flag to indicate iteration in the jet merging process
  Float_t fJdistmaxA; // Amanda maximum jet distance (in m) for jet merging
  Float_t fJdistmaxI; // InIce maximum jet distance (in m) for jet merging
  Int_t fJinvolA;     // Amanda flag to denote maximum jet distance testing inside/outside detector volume
  Int_t fJinvolI;     // InIce flag to denote maximum jet distance testing inside/outside detector volume
  Int_t fMaxmodA;     // The max. number of good fired Amanda OMs for events to get processed
  Int_t fMaxmodI;     // The max. number of good fired InIce DOMs for events to get processed
  Int_t fMinmodA;     // The min. number of good fired Amanda OMs for events to get processed
  Int_t fMinmodI;     // The min. number of good fired InIce DOMs for events to get processed
  Int_t fMaxhitsA;    // The maximum number of good hits per Amanda OM to be processed
  Int_t fMaxhitsI;    // The maximum number of good hits per InIce DOM to be processed
  Int_t fVgroupA;     // Amanda flag to indicate usage of distinct phase and group velocities
  Int_t fVgroupI;     // InIce flag to indicate usage of distinct phase and group velocities
  TString fTrackname; // The name identifier for the produced first guess tracks
  Float_t fCharge;    // User defined charge of the produced first guess tracks

  virtual void Amanda(); // Direct walk reconstruction for Amanda OM signals
  virtual void InIce();  // Direct walk reconstruction for InIce DOM signals
  virtual void AssociateHits(TObjArray& tes,TObjArray& hits,Int_t vgroup,Float_t maxdhit,Float_t& qmax,Int_t& nahmax);// Hit association
  virtual void SelectQvalue(TObjArray& tes,Float_t qmax);                  // TC selection via Q-value
  virtual void ClusterTracks(TObjArray& tes,TObjArray& jets,Float_t tangmax,Int_t tinvol,Float_t tdistmax,Float_t qmax);// Track clustering  
  virtual void MergeJets(TObjArray& jets,Float_t jangmax,Float_t jdistmax,Int_t jinvol,Int_t jiterate);// Jet Merging
  virtual void StoreTracks(TObjArray& jets,Float_t jangmax,TString name,TString title); // Final track storage

 ClassDef(IceDwalk,9) // TTask derived class to perform (improved) direct walk reconstruction
};
#endif
