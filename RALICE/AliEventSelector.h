#ifndef ALIEVENTSELECTOR_H
#define ALIEVENTSELECTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include "AliJob.h"
#include "AliAstrolab.h"
#include "AliEvent.h"

class AliEventSelector : public AliAstrolab
{
 public :
  AliEventSelector(const char* name="AliEventSelector",const char* title="Event selection"); // Constructor
  virtual ~AliEventSelector();                                                               // Destructor
  AliEventSelector(const AliEventSelector& q);               // Copy constructor
  virtual TObject* Clone(const char* name="") const;         // Make a deep copy and provide its pointer
  virtual void Exec(Option_t* opt);                          // Event selection
  void SetSelector(TString type,Int_t flag=1);               // Specify selection types to be used
  void SetLogic(TString type);                               // Set type of decision logic
  void UseTracks(TString name,Int_t n=-1);                   // Specify track names to be used
  void SetAstroMatch(Double_t da,Double_t dt,TString dir);   // Set parameters for Astro matching
  void SetRange(TString type,TString obs,Double_t low,Double_t up); // Set range of various observables
  void SetRange(TString type,TString obs,TString name,Int_t nlow,Int_t nup); // Set range of various observables

 protected :
  Int_t fFirst;                 // Flag to indicate first invokation
  AliEvent* fEvt;               // Pointer to the current event structure
  AliDevice* fParams;           // The device containing all parameter settings and final select flag
  Int_t fTrackflag;             // Flag to indicate usage of individual track selection criteria
  Int_t fEventflag;             // Flag to indicate usage of total event selection criteria
  Int_t fAstroflag;             // Flag to indicate usage of Astrolab selection criteria
  Int_t fLogic;                 // Decision logic (0=unknown  1=and  2=or)
  TObjArray* fUseNames;         // The track names to be used 
  TArrayI* fUseNtk;             // The max. numbers of the various track names to be used
  Int_t fSelect;                // Event selection flag (-1=reject  0=unknown  1=accept)
  Double_t fAstroDa;            // Maximum angular distance (deg.) w.r.t. the reference object
  Double_t fAstroDt;            // Maximum absolute time difference (sec.) w.r.t. the reference signal
  Int_t fAstroDir;              // Direction flag for pointing to external objects
  Float_t fTrackCharges[2];     // Track charge range selections
  Float_t fTrackMasses[2];      // Track mass range selections
  Double_t fTrackMomenta[6];    // Track momentum range selections
  Double_t fTrackEnergies[6];   // Track energy range selections 
  Double_t fTrackRapidities[4]; // Track rapidity range selections
  Int_t fTrackDevices[2];       // Range of number of track associated devices
  TString fTrackDevClass;       // (Derived) class name of track associated devices
  Float_t fEventCharges[2];     // Event charge range selections
  Float_t fEventMasses[2];      // Event mass range selections
  Double_t fEventMomenta[6];    // Event momentum range selections
  Double_t fEventEnergies[6];   // Event energy range selections 
  Int_t fEventDevices[2];       // Range of number of event associated devices
  TString fEventDevClass;       // (Derived) class name of event associated devices
  Int_t fEventTracks[10];       // Range of number of various track types
  TString fEventTrkName;        // Name of the tracks for total number of tracks selection
  void Track(Int_t mode);       // Check criteria for individual track observables
  void Event();                 // Check criteria for total event observables
  void Astro();                 // Check for matches with external (astrophysical) objects

 ClassDef(AliEventSelector,1) // TTask derived class to perform generic event selection
};
#endif
