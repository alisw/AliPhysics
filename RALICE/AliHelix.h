#ifndef ALIHELIX_H
#define ALIHELIX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include "THelix.h"
#include "TObjArray.h"

#include "Ali3Vector.h"
#include "AliTrack.h"
#include "AliEvent.h"
 
class AliHelix : public THelix
{
 public:
  AliHelix();                  // Default constructor
  virtual ~AliHelix();         // Destructor
  AliHelix(const AliHelix& h); // Copy constructor
  void SetB(Ali3Vector& b);    // Set the magnetic field vector in Tesla
  Ali3Vector& GetB();          // Provide the magnetic field vector in Tesla
  void SetTofmax(Float_t tof); // Set maximum time of flight
  Float_t GetTofmax() const;   // Provide the maximum time of flight
  void Display(AliTrack* t,Double_t* range=0,Int_t iaxis=3,Double_t scale=-1);// Show curve for this track
  void Display(AliEvent* e,Double_t* range=0,Int_t iaxis=3,Double_t scale=-1);// Show curves for this event
  void Refresh(Int_t mode=0);  // Refresh the view before drawing the next one
  AliPosition* Extrapolate(AliTrack* t,Double_t* pars=0,Double_t scale=-1); // Extrapolate this track
  void MakeCurve(AliTrack* t,Double_t* range=0,Int_t iaxis=3,Double_t scale=-1); // Helix curve for this track

 protected:
  Ali3Vector fB;                               // The magnetic field vector in Tesla
  Float_t fTofmax;                             // The maximum time of flight
  Int_t fRefresh;                              // Auto-refresh flag for drawings
  TObjArray* fCurves;                          //! Temp. storage for the curves on the drawing
  AliPosition* fExt;                           //! The extrapolation result
 
 ClassDef(AliHelix,1) // Representation and extrapolation of AliTracks in a magnetic field.
};
#endif
