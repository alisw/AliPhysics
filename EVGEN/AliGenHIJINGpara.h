#ifndef ALIGENHIJINGPARA_H
#define ALIGENHIJINGPARA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Parameterisation of pi and K, eta and pt distributions
// used for the ALICE TDRs.
// Author: andreas.morsch@cern.ch

#include "AliGenerator.h"
class TF1;
class AliDecayerPythia;

class AliGenHIJINGpara : public AliGenerator
{
 public:

  AliGenHIJINGpara();
  AliGenHIJINGpara(Int_t npart);
  AliGenHIJINGpara(const AliGenHIJINGpara &HIJINGpara);
     
  virtual ~AliGenHIJINGpara();
  virtual void SetCutVertexZ(Float_t cut=999999.) {fCutVertexZ = cut;}
  virtual void Generate();
  virtual void Init();
  AliGenHIJINGpara & operator=(const AliGenHIJINGpara & rhs);
  virtual void SetPtRange(Float_t ptmin = 0., Float_t ptmax=15.);
  virtual void SetPi0Decays(Bool_t flag = kFALSE) {fPi0Decays = flag;}
 private:
  void DecayPi0(Float_t* orig, Float_t * p);
 protected:
  Int_t   fNt;         // CurrentTrack;
  Float_t fCutVertexZ; // Vertex truncation
  Bool_t  fPi0Decays;  // Flag for pi0 decays
  TF1* fPtpi;          // Parametrised pt distribution for pi
  TF1* fPtka;          // Parametrised pt distribution for ka
  TF1* fETApic;        // Parametrised eta distribution for pi
  TF1* fETAkac;        // Parametrised eta distribution fro ka
  AliDecayerPythia  *fDecayer; // ! Pointer to pythia object for decays
  ClassDef(AliGenHIJINGpara,3) // Hijing parametrisation generator
};
#endif










