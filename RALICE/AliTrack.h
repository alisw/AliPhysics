#ifndef ALITRACK_H
#define ALITRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "TObject.h"
#include "TObjArray.h"
 
#include "AliSignal.h"
#include "AliBoost.h"
#include "AliPosition.h"
 
class AliTrack : public TObject,public Ali4Vector
{
 public:
  AliTrack();                       // Default constructor
  ~AliTrack();                      // Destructor
  void Reset();                     // Reset all values to 0
  void Set4Momentum(Ali4Vector& p); // Set track 4-momentum
  void Set3Momentum(Ali3Vector& p); // Set track 3-momentum
  void SetMass(Double_t m,Double_t dm=0); // Set particle mass and error
  void SetCharge(Float_t q);        // Set particle charge
  void Info(TString f="car");       // Print track information for coord. frame f
  void List(TString f="car");       // Print track and decay level 1 information for coord. frame f
  void ListAll(TString f="car");    // Print track and all decay level information for coord. frame f
  Ali3Vector Get3Momentum();        // Provide track 3-momentum
  Double_t GetMomentum();           // Provide value of track 3-momentum
  Double_t GetMass();               // Provide particle mass
  Float_t GetCharge();              // Provide particle charge
  Double_t GetEnergy();             // Provide particle total energy
  void Decay(Double_t m1,Double_t m2,Double_t thcms,Double_t phicms); // Perform 2-body decay
  Int_t GetNdecay();                // Provide number of decay products
  AliTrack* GetDecayTrack(Int_t j); // Access to decay produced track number j
  void AddSignal(AliSignal& s);     // Relate an AliSignal to this track
  void RemoveSignal(AliSignal& s);  // Remove related AliSignal from this track
  Int_t GetNsignals();              // Provide number of related AliSignals
  AliSignal* GetSignal(Int_t j);    // Access to the related AliSignal number j
  void SetBeginPoint(AliPosition p);// Set the track begin-point
  AliPosition GetBeginPoint();      // Provide the track begin-point
  void SetEndPoint(AliPosition p);  // Set the track end-point
  AliPosition GetEndPoint();        // Provide the track end-point
 
 protected:
  Float_t fQ;          // The charge of the particle
  Int_t fNdec;         // The number of decay products
  TObjArray* fDecays;  // The array of decay produced tracks
  Int_t fNsig;         // The number of related AliSignals
  TObjArray* fSignals; // The array of related AliSignals
  AliPosition fBegin;  // The begin-point of the track 
  AliPosition fEnd;    // The end-point of the track 

 private:
  void Dump(AliTrack* t,Int_t n,TString f); // Recursively print all decay levels
 
 ClassDef(AliTrack,1) // Handling of the attributes of a reconstructed particle track.
};
#endif
