#ifndef ALIGEVSIMPARTICLE_H
#define ALIGEVSIMPARTICLE_H

#include "TObject.h"


class AliGeVSimParticle : public TObject {


  Int_t fPDG;            // Particle type code

  Float_t fN;            // Multiplicity (subject to scalling)
  Float_t fT;            // Slope Parameter (subject to scalling)
  Float_t fSigmaY;       // Rapidity Width
  Float_t fExpansion;    // Expansion Velocity in c units (subject to scalling)

  Float_t fV1;           // Direct Flow coefficient (subject to scalling)
  Float_t fV2;           // Elliptical flow coefficient (subject to scalling)

 public:

  ////////////////////////////////////////////////////////////////////////////

  AliGeVSimParticle() {}
  AliGeVSimParticle(Int_t pdg); 
  AliGeVSimParticle(Int_t pdg, Int_t n, 
		    Float_t T, Float_t dY = 1., Float_t exp=0.);

  ~AliGeVSimParticle() {}

  ////////////////////////////////////////////////////////////////////////////

  Int_t GetPdgCode() {return fPDG;}


  Float_t GetMultiplicity() {return fN;}
  Float_t GetTemperature() {return fT;}
  Float_t GetSigmaY() {return fSigmaY;}
  Float_t GetExpansionVelocity() {return fExpansion;}
  
  void SetMultiplicity(Float_t n) {fN = n;}
  void SetExpansionVelocity(Float_t vel) {fExpansion = vel;}

  // Flow

  void SetDirectFlow(Float_t v1) {fV1 = v1;}
  void SetEllipticalFlow(Float_t v2) {fV2 = v2;}

  Float_t GetDirectFlow() {return fV1;}
  Float_t GetEllipticalFlow() {return fV2;}


  ////////////////////////////////////////////////////////////////////////////

  ClassDef(AliGeVSimParticle, 1)

};


#endif
