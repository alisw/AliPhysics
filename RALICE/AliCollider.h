#ifndef ALICOLLIDER_H
#define ALICOLLIDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id: AliCollider.h,v 1.3 2002/12/11 14:45:12 nick Exp $

#include "TPythia6.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h" 

#include "AliEvent.h"
#include "AliRandom.h"
 
class AliCollider : public TPythia6
{
 public:
  AliCollider();                                        // Default constructor
  virtual ~AliCollider();                               // Default destructor
  void SetOutputFile(TString name);                     // Initialise the ROOT output data file
  void SetVertexMode(Int_t mode);                       // Select mode for (sec.) vertex structure creation
  Int_t GetVertexMode();                                // Provide vertex structure creation mode
  void SetResolution(Double_t res);                     // Set resolution (in cm) for resolving (sec.) vertices
  Double_t GetResolution();                             // Provide (sec.) vertex resolving resolution (in cm)
  void SetRunNumber(Int_t run);                         // Set user defined run number
  Int_t GetRunNumber();                                 // Provide the user defined run number
  void SetPrintFreq(Int_t n);                           // Set print frequency for every 'n' events
  Int_t GetPrintFreq();                                 // Provide the print frequency
  void Init(char* frame,char* beam,char* target,Float_t win); // Standard Pythia initialisation
  void Init(char* frame,Int_t zp,Int_t ap,Int_t zt,Int_t at,Float_t win); // Nucl-Nucl initialisation
  void MakeEvent(Int_t npt,Int_t mlist=-1,Int_t medit=1);// Generate a single event with npt participant nucleons
  void EndRun();                                        // Properly close all buffers and output file
  AliEvent* GetEvent();                                 // Provide pointer to the generated event structure

 protected:
  Int_t fVertexmode;    // The vertex structure creation mode
  Double_t fResolution; // The resolution (in cm) for resolving (sec.) vertices 
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

  TFile* fOutFile;      // The user defined output data file 
  TTree* fOutTree;      // The standard ROOT output tree

  void GetFractions(Float_t zp,Float_t ap,Float_t zt,Float_t at); // Determine various N-N collision fractions

 ClassDef(AliCollider,3) // Pythia based universal physics event generator
};
#endif
