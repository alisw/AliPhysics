#ifndef ALISIMPLEGEN_H
#define ALISIMPLEGEN_H
///////////////////////////////////////////////////////////
//                                                       //
//  Class to generate the particles for the MC           //
//  The base class is empty                              //
//                                                       //
///////////////////////////////////////////////////////////

#include "AliGenerator.h"
#include "AliPDG.h"
#include "TF1.h"


class AliGenHIJINGpara : public AliGenerator
{
 
protected:

  TF1* fPtpi; // Parametrised pt distribution for pi
  TF1* fPtka; // Parametrised pt distribution for ka
  TF1* fETApic; // Parametrised eta distribution for pi
  TF1* fETAkac; // Parametrised eta distribution fro ka

public:

  AliGenHIJINGpara();
  AliGenHIJINGpara(Int_t npart);
  virtual ~AliGenHIJINGpara();
  virtual void Generate();
  virtual void Init();

  ClassDef(AliGenHIJINGpara,1) // Hijing parametrisation generator
};

class AliGenFixed : public AliGenerator
{
 
protected:

  Int_t fIpart; // Particle type

public:

  AliGenFixed();
  AliGenFixed(Int_t npart);
  virtual ~AliGenFixed() {}
  virtual void Generate();
  virtual void Init() {}
  virtual void SetSigma(Float_t, Float_t, Float_t);
  //
  virtual void SetMomentum(Float_t pmom) {fPMin=pmom; fPMax=pmom;}
  virtual void SetPhi(Float_t phi) {fPhiMin=phi*TMath::Pi()/180; fPhiMax=phi*TMath::Pi()/180;}
  virtual void SetTheta(Float_t theta) {fThetaMin=theta*TMath::Pi()/180; fThetaMax=theta*TMath::Pi()/180;}
  virtual void SetPart(Int_t part) {fIpart=part;}

  ClassDef(AliGenFixed,1) // Single particle generator
};


class AliGenBox : public AliGenerator
{
 
protected:

  Int_t fIpart; // Particle type

public:

  AliGenBox();
  AliGenBox(Int_t npart);
  virtual ~AliGenBox() {}
  virtual void Generate();
  virtual void Init() {}
  virtual void SetPart(Int_t part) {fIpart=part;}

  ClassDef(AliGenBox,1) // Square box random generator
};

#endif
