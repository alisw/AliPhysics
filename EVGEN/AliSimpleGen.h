#ifndef ALISIMPLEGEN_H
#define ALISIMPLEGEN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////
//                                                       //
//  Class to generate the particles for the MC           //
//  The base class is empty                              //
//                                                       //
///////////////////////////////////////////////////////////

#include "AliGenerator.h"
class TF1;



class AliGenHIJINGpara : public AliGenerator
{
 public:

  AliGenHIJINGpara();
  AliGenHIJINGpara(Int_t npart);
  virtual ~AliGenHIJINGpara();
  virtual void Generate();
  virtual void Init();

 protected:

  TF1* fPtpi; // Parametrised pt distribution for pi
  TF1* fPtka; // Parametrised pt distribution for ka
  TF1* fETApic; // Parametrised eta distribution for pi
  TF1* fETAkac; // Parametrised eta distribution fro ka

  ClassDef(AliGenHIJINGpara,1) // Hijing parametrisation generator
};

class AliGenFixed : public AliGenerator
{
 public:
  AliGenFixed();
  AliGenFixed(Int_t npart);
  virtual ~AliGenFixed() {}
  virtual void Generate();
  virtual void Init() {}
  virtual void SetSigma(Float_t sx, Float_t sy, Float_t sz);
  virtual void SetMomentum(Float_t pmom) {fPMin=pmom; fPMax=pmom;}
  virtual void SetPhi(Float_t phi) {fPhiMin=phi*TMath::Pi()/180; fPhiMax=phi*TMath::Pi()/180;}
  virtual void SetTheta(Float_t theta) {fThetaMin=theta*TMath::Pi()/180; fThetaMax=theta*TMath::Pi()/180;}
  virtual void SetPart(Int_t part) {fIpart=part;}
 
protected:

  Int_t fIpart; // Particle type

  ClassDef(AliGenFixed,1) // Single particle generator
};


class AliGenBox : public AliGenerator
{
 public:

  AliGenBox();
  AliGenBox(Int_t npart);
  virtual ~AliGenBox() {}
  virtual void Generate();
  virtual void Init();
  virtual void SetPart(Int_t part) {fIpart=part;}
protected:

  Int_t fIpart; // Particle type

  ClassDef(AliGenBox,1) // Square box random generator
};

#endif
