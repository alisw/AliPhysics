#ifndef ALIJET_H
#define ALIJET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
 
#include <math.h>
 
#include "TNamed.h"
#include "TObjArray.h"

#include "Ali4Vector.h"
#include "AliTrack.h"
 
class AliJet : public TNamed,public Ali4Vector
{
 public:
  AliJet();                                // Default constructor
  AliJet(Int_t n);                         // Create a Jet to hold initially n Tracks
  virtual ~AliJet();                       // Default destructor
  AliJet(const AliJet& j);                 // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer
  virtual void SetOwner(Bool_t own=kTRUE); // Set ownership of all added objects
  virtual void Reset();                    // Reset all values
  void AddTrack(AliTrack& t);              // Add a track to the jet
  void AddTrack(AliTrack* t) { AddTrack(*t); }
  virtual void Data(TString f="car");      // Print jet information in coordinate frame f 
  virtual void List(TString f="car");      // Print jet prim. track information for coord. frame f
  virtual void ListAll(TString f="car");   // Print jet prim. and decay track information for coord. frame f
  Double_t GetEnergy();                    // Provide the total jet energy
  Double_t GetMomentum();                  // Provide the value of the total jet 3-momentum
  Ali3Vector Get3Momentum() const;         // Provide the total jet 3-momentum
  Double_t GetInvmass();                   // Provide the invariant mass  
  Float_t GetCharge() const;               // Provide the total charge of the jet
  Int_t GetNtracks() const;                // Return the number of tracks in the jet
  AliTrack* GetTrack(Int_t i) const;       // Provide i-th track of the jet (1=first track)
  AliTrack* GetIdTrack(Int_t id) const;    // Provide the track with user identifier "id"
  Double_t GetPt();                        // Provide trans. momentum w.r.t. z-axis
  Double_t GetPl();                        // Provide long. momentum w.r.t. z-axis
  Double_t GetEt();                        // Provide trans. energy w.r.t. z-axis
  Double_t GetEl();                        // Provide long. energy w.r.t. z-axis
  Double_t GetMt();                        // Provide trans. mass w.r.t. z-axis
  Double_t GetRapidity();                  // Provide rapidity value w.r.t. z-axis
  void SetTrackCopy(Int_t j);              // (De)activate creation of private copies in fTracks
  Int_t GetTrackCopy() const;              // Provide TrackCopy flag value      
  void SetId(Int_t id);                    // Set the user defined identifier
  Int_t GetId() const;                     // Provide the user defined identifier

 protected:
  void Init();                           // Initialisation of pointers etc...
  void SetNtinit(Int_t n=2);             // Set the initial max. number of tracks for this Jet
  void AddTrack(AliTrack& t,Int_t copy); // Internal memberfunction to add a track to the jet
  void AddTrack(AliTrack* t,Int_t copy) { AddTrack(*t,copy); }
  Int_t fNtinit;                         // The initial max. number of tracks for this jet
  Int_t fNtmax;                          // The maximum number of tracks for this Jet
  Float_t fQ;                            // The total charge of the jet 
  Int_t fNtrk;                           // The number of tracks in the jet
  TObjArray* fTracks;                    // Array to hold the pointers to the tracks of the jet
  Int_t fTrackCopy;                      // Flag to denote creation of private copies in fTracks
  Int_t fUserId;                         // The user defined identifier
 
 ClassDef(AliJet,9) // Creation and investigation of a jet of particle tracks.
};
#endif
