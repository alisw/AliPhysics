#ifndef ALIINVMASS_H
#define ALIINVMASS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
// Class AliInvmass
// Construction of invariant mass and combinatorial background.
//
// Example :
// ---------
//
// TObjArray* photons=new TObjArray(); // Array with photon tracks for pi0 rec.
//
// // Code to create some photon tracks from pi0 decays
// Int_t ntracks=200;
// for (Int_t i=0; i<ntracks; i++)
// {
//  photons->Add(new Alitrack);
//  ...
//  ...
//  ...
// }  
//
// // Perform the invariant mass and comb. bkg. reconstruction
//
// TObjArray* allm=q.Invmass(photons,photons); // All reconstructed invariant masses
//
// TH1F* hall=new TH1F("hall","hall",200,0,2); // Histo with M_inv of all combinations
//
// Int_t nall=0;
// if (allm) nall=allm->GetEntries();
//
// AliTrack* t;
// Float_t minv;
// for (Int_t j=0; j<nall; j++)
// {
//  t=(AliTrack*)allm->At(j);
//  if (t)
//  {
//   minv=t->GetMass();
//   hall->Fill(minv);
//  }
// }
//
// TObjArray* bkgm=q.CombBkg(photons,photons); // Reconstructed comb. background
//
// TH1F* hbkg=new TH1F("hbkg","hbkg",200,0,2); // Histo with M_inv. of comb. background
//
// Int_t nbkg=0;
// if (bkgm) nbkg=bkgm->GetEntries();
//
// for (Int_t j=0; j<nbkg; j++)
// {
//  t=(AliTrack*)bkgm->At(j);
//  if (t)
//  {
//   minv=t->GetMass();
//   hbkg->Fill(minv);
//  }
// }
//
// TH1F* hsig=new TH1F("sig","sig",200,0,2);   // Histo with the bkg. subtracted signal
// hsig->Sumw2();
// hsig->Add(hall,hbkg,1,-1);
//
//
// Note : By default the storage of the reconstructed information is performed
//        in separate TObjArrays for the signal and comb. background resp.
//        In order to limit the memory usage, AliInvmass::SetStorageMode(1) may be
//        used to activate only a single TObjArray to store the reconstructed information.
//        Consequently, the following statements 
//
//         TObjArray* allm=q.Invmass(photons,photons);
//         TObjArray* bkgm=q.CombBkg(photons,photons);
//
//        will result in the fact that after he invokation of CombBkg
//        the information of "allm" is lost due to the fact that the storage is
//        is re-used for "bkgm" in case the "single storage" option has been selected.
//        Usage of the, in that case invalid, pointer "allm" may cause your
//        program to crash.
//
//        * Thus : In case of single storage usage, all invokations of the returned
//                 array pointer have to be completed before invoking any memberfunction
//                 of the same AliInvmass object again.
//        
//        
//
//--- NvE 12-apr-1999 UU-SAP Utrecht
////////////////////////////////////////////////////////////////////////////////

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

 ClassDef(AliInvmass,1) // Class definition to enable ROOT I/O
};
#endif
