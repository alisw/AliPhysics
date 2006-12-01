#ifndef ALICOLLIDER2_H
#define ALICOLLIDER2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id: AliCollider.h,v 1.9 2004/05/04 15:33:04 nick Exp $

#include "TPythia6.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TArrayI.h" 

#include "AliEvent.h"
#include "AliRandom.h"
 
class AliCollider : public TPythia6
{
 public:
  AliCollider();                                        // Default constructor
  virtual ~AliCollider();                               // Default destructor
  void SetOutputFile(TString name);                     // Initialise the ROOT output data file
  void SetVertexMode(Int_t mode);                       // Select mode for (sec.) vertex structure creation
  Int_t GetVertexMode() const;                          // Provide vertex structure creation mode
  void SetResolution(Double_t res);                     // Set resolution (in meter) for resolving (sec.) vertices
  Double_t GetResolution() const;                       // Provide (sec.) vertex resolving resolution (in meter)
  void SetRunNumber(Int_t run);                         // Set user defined run number
  Int_t GetRunNumber() const;                           // Provide the user defined run number
  void SetPrintFreq(Int_t n);                           // Set print frequency for every 'n' events
  Int_t GetPrintFreq() const;                           // Provide the print frequency
  void Init(char* frame,char* beam,char* target,Float_t win); // Standard Pythia initialisation
  void Init(char* frame,Int_t zp,Int_t ap,Int_t zt,Int_t at,Float_t win); // Nucl-Nucl initialisation
  void SetStable(Int_t id,Int_t mode=1);                // Declare a certain particle as stable or not
  void SelectEvent(Int_t id);                           // Select only events containing specified particles
  Int_t GetSelectionFlag() const;                       // Return the selection flag for this event
  void MakeEvent(Int_t npt=0,Int_t mlist=-1,Int_t medit=1);// Generate a single event with npt participant nucleons
  void EndRun();                                        // Properly close all buffers and output file
  AliEvent* GetEvent(Int_t select=0) const;             // Provide pointer to the generated event structure
  void SetSpectatorPmin(Float_t pmin);                  // Set minimal momentum for spectator track to be stored
  Float_t GetSpectatorPmin() const;                     // Provide the minimal momentum for spectator tracks
  void SetUserControl(Int_t flag);                      // Selection of full user control w.r.t. MC parameters. 
  Int_t GetUserControl() const;                         // Provide the value of the user control flag.
  void SetElastic(Int_t flag);                          // Selection flag for elastic and diffractive processes.
  Int_t GetElastic() const;                             // Provide the value of the elastic selection flag.

 protected:
  Int_t fVertexmode;    // The vertex structure creation mode
  Double_t fResolution; // The resolution (in meter) for resolving (sec.) vertices 
  Int_t fRunnum;        // The user defined run number
  Int_t fEventnum;      // The automatically updated event number
  Int_t fPrintfreq;     // The user selected print frequency
  char* fFrame;         // The Pythia frame indicator
  Float_t fWin;         // The Pythia energy indicator
  Int_t fNucl;          // Flag to denote nucleus-nucleus (1) or standard Pythia (0) running
  Int_t fZproj;         // Z of the projectile particle
  Int_t fAproj;         // A of the projectile particle
  Int_t fZtarg;         // Z of the target particle
  Int_t fAtarg;         // A of the target particle
  Float_t fFracpp;      // Fraction of p+p collisions
  Float_t fFracnp;      // Fraction of n+p collisions
  Float_t fFracpn;      // Fraction of p+n collisions
  Float_t fFracnn;      // Fraction of n+n collisions
  AliRandom fRan;       // Random number generator
  AliEvent* fEvent;     // The produced event structure
  Float_t fSpecpmin;    // The minimal momentum for spectator tracks to be stored
  Int_t fUserctrl;      // Flag to denote the user control selection w.r.t. MC parameters
  Int_t fElastic;       // Flag to denote inclusion of elastic and difractive processes.

  TFile* fOutFile;      // The user defined output data file 
  TTree* fOutTree;      // The standard ROOT output tree

  TArrayI* fSelections; // The particle KC codes for event selection
  Int_t fSelect;        // Flag to indicate whether the total event is selected (1) or not (0)

  Int_t IsSelected();   // Check whether (sub)event passed the selection criteria
  void GetFractions(Float_t zp,Float_t ap,Float_t zt,Float_t at); // Determine various N-N collision fractions
  TString GetPyname(Int_t kf); // Provide the correctly truncated Pythia particle name for PDG code kf  

 ClassDef(AliCollider,8) // Pythia based universal physics event generator
};
#endif
