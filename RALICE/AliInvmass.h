#ifndef ALIINVMASS_H
#define ALIINVMASS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
#include "TObjArray.h"

#include "AliRandom.h"
#include "AliTrack.h"

class AliInvmass : public TObject
{
 public:
  AliInvmass();                                    // Default constructor
  ~AliInvmass();                                   // Destructor
  void SetStorageMode(Int_t m);                    // Set storage mode (1=single, 2=multiple)
  void SetThetaSwitch(Int_t i=1);                  // Enable (1/0) new theta for comb. bkg. reco.
  void SetPhiSwitch(Int_t i=1);                    // Enable (1/0) new phi for comb. bkg. reco.
  Int_t GetStorageMode();                          // Provide storage mode
  Int_t GetThetaSwitch();                          // Provide theta switch flag
  Int_t GetPhiSwitch();                            // Provide phi switch flag
  TObjArray* Invmass(TObjArray* a1,TObjArray* a2); // Two-particle inv. mass reco.
  TObjArray* CombBkg(TObjArray* a1,TObjArray* a2); // Two-particle comb. background reco.

 protected:
  Double_t fPi;     // Value of pi
  Int_t fMode;      // Storage mode for signal and bkg. results (2=separate arrays)
  Int_t fBkg;       // Flag to denote comb. background processing
  AliRandom fRndm;  // The random number generator for the comb. bkg. reconstruction
  Int_t fNewtheta;  // Flag to denote enabling of switching theta for comb. bkg. reco.
  Int_t fNewphi;    // Flag to denote enabling of switching phi for comb. bkg. reco.
  TObjArray* fMinv; // Array with reconstructed invariant mass 'tracks'
  TObjArray* fMbkg; // Array with reconstructed comb. background 'tracks'

 private:
  void Combine(TObjArray* a1,TObjArray* a2); // Make two-particle combinations

 ClassDef(AliInvmass,1) // Construction of invariant mass and combinatorial background.
};
#endif
