#ifndef ALIJET_H
#define ALIJET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
 
#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
#include "TObjArray.h"

#include "Ali4Vector.h"
#include "AliTrack.h"
 
class AliJet : public TObject,public Ali4Vector
{
 public:
  AliJet();                          // Default constructor
  AliJet(Int_t n);                   // Create a Jet to hold initially n Tracks
  ~AliJet();                         // Default destructor
  void Reset();                      // Reset all values
  virtual void Add(AliTrack& t);     // Add a track to the jet
  virtual void Add(AliTrack* t) { Add(*t); }
  void Info(TString f);              // Print jet information in coordinate frame f 
  void List(TString f="car");        // Print jet prim. track information for coord. frame f
  void ListAll(TString f="car");     // Print jet prim. and decay track information for coord. frame f
  Double_t GetEnergy();              // Provide the total jet energy
  Double_t GetMomentum();            // Provide the value of the total jet 3-momentum
  Ali3Vector Get3Momentum();         // Provide the total jet 3-momentum
  Double_t GetInvmass();             // Provide the invariant mass  
  Float_t GetCharge();               // Provide the total charge of the jet
  Int_t GetNtracks();                // Return the number of tracks in the jet
  AliTrack* GetTrack(Int_t i);       // Provide i-th track of the jet (1=first track)

 protected:
  void SetNtinit(Int_t n=2); // Set the initial max. number of tracks for this Jet
  Int_t fNtinit;             // The initial max. number of tracks for this jet
  Int_t fNtmax;              // The maximum number of tracks for this Jet
  Float_t fQ;                // The total charge of the jet 
  Int_t fNtrk;               // The number of tracks in the jet
  TObjArray* fTracks;        // Array to hold the pointers to the tracks of the jet
 
 ClassDef(AliJet,1) // Creation and investigation of a jet of particle tracks.
};
#endif
