#ifndef ALIGEVSIMPARTICLE_H
#define ALIGEVSIMPARTICLE_H

////////////////////////////////////////////////////////////////////////////////
//
// AliGeVSimParticle is a helper class for GeVSim (AliGenGeVSim) event generator.
// An object of this class represents one particle type and contain 
// information about particle type thermal parameters.
//
// For examples, parameters and testing macros refer to:
// http:/home.cern.ch/radomski
//
// Author:
// Sylwester Radomski,
// GSI, March 2002
// 
// S.Radomski@gsi.de
//
////////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliGeVSimParticle : public TObject {

 public:
  
  ////////////////////////////////////////////////////////////////////////////
  
  AliGeVSimParticle() {}
  AliGeVSimParticle(Int_t pdg); 
  AliGeVSimParticle(Int_t pdg, Int_t n, 
		    Float_t T, Float_t dY = 1., Float_t exp=0.);
  
  ~AliGeVSimParticle() {}
  
  ////////////////////////////////////////////////////////////////////////////
  
  Int_t GetPdgCode() const {return fPDG;}
  
  Float_t GetMultiplicity() const {return fN;}
  Float_t GetTemperature() const {return fT;}
  Float_t GetSigmaY() const {return fSigmaY;}
  Float_t GetExpansionVelocity() const {return fExpansion;}
  
  void SetMultiplicity(Float_t n) {fN = n;}
  void SetExpansionVelocity(Float_t vel) {fExpansion = vel;}
  
  // Flow
  
  void SetDirectedFlow(Float_t v11, Float_t v12=0, Float_t v13=1, Float_t v14=0);
  void SetEllipticFlow(Float_t v21, Float_t v22=0, Float_t v23=0);
  
  Float_t GetDirectedFlow(Float_t pt, Float_t y);
  Float_t GetEllipticFlow(Float_t pt, Float_t y);
  
  Float_t GetDirectedFlow();
  Float_t GetEllipticFlow();

  
  ////////////////////////////////////////////////////////////////////////////
  
 private:
  
  Int_t fPDG;            // Particle type code
  
  Float_t fN;            // Multiplicity (subject to scalling)
  Float_t fT;            // Slope Parameter (subject to scalling)
  Float_t fSigmaY;       // Rapidity Width
  Float_t fExpansion;    // Expansion Velocity in c units (subject to scalling)
  
  Float_t fV1[4];        // Direct Flow coefficient parameters (subject to scalling)
  Float_t fV2[3];        // Elliptical flow coefficient parameters (subject to scalling)
  
 public:
  
  ClassDef(AliGeVSimParticle, 1)
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
