////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "AliGeVSimParticle.h"

ClassImp(AliGeVSimParticle);

////////////////////////////////////////////////////////////////////////////////////////////////////

AliGeVSimParticle::AliGeVSimParticle(Int_t pdg, Int_t n, Float_t T, Float_t dY, Float_t exp) {
  //
  //  pdg - Particle type code in PDG standard (see: http://pdg.lbl.gov)
  //  n   - Multiplicity of particle type
  //  T   - Inverse slope parameter ("temperature")
  //  dY  - Raridity Width (only for model 1)
  //  exp - expansion velocity (only for model 4) 
  
  fPDG = pdg;
  fN = n;
  fT = T;
  fSigmaY = dY;
  fExpansion = exp;

  fV1[0] = fV1[1] = fV1[2] = fV1[3] = 0.;
  fV2[0] = fV2[1] = fV2[2] = 0.;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

AliGeVSimParticle::AliGeVSimParticle(Int_t pdg) {
  //
  //  pdg - Particle type code in PDG standard (see: http://pdg.lbl.gov)
  // 

  fPDG = pdg;
  fN = 0;
  fT = 0.;
  fSigmaY = 0.;
  fExpansion = 0.;
  
  fV1[0] = fV1[1] = fV1[2] = fV1[3] = 0.;
  fV2[0] = fV2[1] = fV2[2] = 0.; 
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGeVSimParticle::SetDirectedFlow(Float_t v11, Float_t v12, Float_t v13, Float_t v14) {
  //
  // Set Directed Flow parameters.
  // Acctual flow coefficient is calculated as follows
  //
  // V1(Pt,Y) = (V11 + V12*Pt) * sign(Y) * (V13 + V14 * Y^3)
  //
  // where sign = 1 for Y > 0 and -1 for Y < 0
  // 
  // Defaults values
  // v12 = v14 = 0
  // v13 = 1
  // 
  // Note 1: In many cases it is sufficient to set v11 only.
  // Note 2: Be carefull with parameter v14
  // 


  fV1[0] = v11;
  fV1[1] = v12;
  fV1[2] = v13;
  fV1[3] = v14;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGeVSimParticle::SetEllipticFlow(Float_t v21, Float_t v22, Float_t v23) {
  //
  // Set Elliptic Flow parameters.
  // Acctual flow coefficient is calculated as follows
  //
  // V2 = (V21 + V22 * Pt^2) * exp( -V23 * Y^2)
  // 
  // Default values:
  // v22 = v23 = 0
  //
  // Note: In many cases it is sufficient to set v21 only
  //
  
  fV2[0] = v21;
  fV2[1] = v22;
  fV2[2] = v23;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t AliGeVSimParticle::GetDirectedFlow(Float_t pt, Float_t y) {
  //
  // Return coefficient of a directed flow for a given pt and y.
  // For coefficient calculation method refer to SetDirectedFlow()
  // 
  
  Float_t v;
  
  v = (fV1[0] + fV1[1]* pt) * TMath::Sign(1.,y) *
    (fV1[2] + fV1[3] * TMath::Abs(y*y*y) );

  return v;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t AliGeVSimParticle::GetEllipticFlow(Float_t pt, Float_t y) {
  //
  // Return coefficient of a elliptic flow for a given pt and y.
  // For coefficient calculation method refer to SetEllipticFlow()
  // 
    
  return  (fV2[0] + fV2[2] * pt * pt) * TMath::Exp( -fV2[3] * y*y );
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t AliGeVSimParticle::GetDirectedFlow() {
  //
  // Simplified version of GetDirectedFlow(pt,y) for backward compatibility
  // Return fV1[0]
  //
  
  return fV1[0];
}

////////////////////////////////////////////////////////////////////////////////////////////////////
Float_t AliGeVSimParticle::GetEllipticFlow() {
  //
  // Simplified version of GetEllipticFlow(pt,y) for backward compatibility
  // Return fV2[0]
  //
  
  return fV2[0];
}

////////////////////////////////////////////////////////////////////////////////////////////////////














