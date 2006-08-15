/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//////////////////////////////////////////////////////////////////////////////
//
// AliGeVSimParticle is a helper class for GeVSim (AliGenGeVSim) event generator.
// An object of this class represents one particle type and contain 
// information about particle type thermal parameters.
//
//////////////////////////////////////////////////////////////////////////////
//
// For examples, parameters and testing macros refer to:
// http:/home.cern.ch/radomski
// 
// for more detailed description refer to ALICE NOTE
// "GeVSim Monte-Carlo Event Generator"
// S.Radosmki, P. Foka.
//  
// Author:
// Sylwester Radomski,
// GSI, March 2002
//  
// S.Radomski@gsi.de
//
////////////////////////////////////////////////////////////////////////////////
//
// Updated and revised: September 2002, S. Radomski, GSI
//
////////////////////////////////////////////////////////////////////////////////


#include "TMath.h"
#include "AliGeVSimParticle.h"

ClassImp(AliGeVSimParticle)


////////////////////////////////////////////////////////////////////////////////////////////////////
AliGeVSimParticle::AliGeVSimParticle():
    fPDG(0),
    fModel(0),
    fN(0),
    fMultTotal(kTRUE),
    fIsSetMult(kFALSE),
    fT(0.),
    fSigmaY(0.),
    fExpansion(0.),
    fIsDirectedSimple(kTRUE),
    fIsEllipticSimple(kTRUE),
    fIsEllipticOld(kFALSE)
{
    // Default constructor
}

AliGeVSimParticle::AliGeVSimParticle(Int_t pdg, Int_t model, Float_t multiplicity,
				     Float_t T, Float_t dY, Float_t exp):
    fPDG(pdg),
    fModel(model),
    fN(multiplicity),
    fMultTotal(kTRUE),
    fIsSetMult(kFALSE),
    fT(T),
    fSigmaY(dY),
    fExpansion(exp),
    fIsDirectedSimple(kTRUE),
    fIsEllipticSimple(kTRUE),
    fIsEllipticOld(kFALSE)
{
  //
  //  pdg          - Particle type code in PDG standard (see: http://pdg.lbl.gov)
  //  model        - momentum distribution model (1 - 7)
  //  multiplicity - multiplicity of particle type
  //  T            - Inverse slope parameter ("temperature")
  //  dY           - Raridity Width (only for model 1)
  //  exp          - expansion velocity (only for model 4) 
  fV1[0] = fV1[1] = fV1[2] = fV1[3] = 0.;
  fV2[0] = fV2[1] = fV2[2] = 0.;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

AliGeVSimParticle::AliGeVSimParticle(Int_t pdg, Int_t model, Float_t multiplicity):
    fPDG(pdg),
    fModel(model),
    fN(multiplicity),
    fMultTotal(kTRUE),
    fIsSetMult(kFALSE),
    fT(0.),
    fSigmaY(0.),
    fExpansion(0.),
    fIsDirectedSimple(kTRUE),
    fIsEllipticSimple(kTRUE),
    fIsEllipticOld(kFALSE)
 {
  //
  // pdg - Particle type code in PDG standard (see: http://pdg.lbl.gov)
  //  
  // Note that multiplicity can be interpreted by GeVSim 
  // either as Total multiplicity in the acceptance or dN/dY
  // 
  fV1[0] = fV1[1] = fV1[2] = fV1[3] = 0.;
  fV2[0] = fV2[1] = fV2[2] = 0.; 
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void  AliGeVSimParticle::SetModel(Int_t model) {
  //
  // Set Model (1-7) 
  // For details about standard and custom models refer to ALICE NOTE
  //

  if (model < 1 || model > 7) 
    Error("SetModel","Model Id ( %d ) out of range [1..7]", model);

  fModel = model;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void  AliGeVSimParticle::SetMultiplicity(Float_t mult) {
  //
  // Set multiplicity. The value is interpreted either as a total multiplciity
  // in the acceptance or as a multiplicity density - dN/dY at midrapidity
  //  

  fN = mult;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGeVSimParticle::SetMultTotal(Bool_t isTotal) {
  //
  // Switch between total multiplicity (kTRUE) and 
  // multiplciity density (kFALSE)
  //
  // If this method is used its overrides mode in AliGenGeVSim 
  //
  
  fMultTotal = isTotal;
  fIsSetMult = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
 
void AliGeVSimParticle::SetDirectedSimple(Float_t v1) {
  //
  // Set directed flow coefficient to a value independent
  // of transverse momentum and rapidity
  //

  fV1[0] = v1;
  fIsDirectedSimple = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGeVSimParticle::SetEllipticSimple(Float_t v2) {
  //
  // Set elliptic flow coefficient to a value independent
  // of transverse momentum and rapidity
  //

  fV2[0] = v2;
  fIsEllipticSimple = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t AliGeVSimParticle::IsFlowSimple() const
{
  //
  // Function used by AliGenGeVSim 
  //
  // Returns true if both Elliptic and Directed flow has a simple model.
  // If at least one is parametrised returns false. 
  // 

  return (fIsDirectedSimple && fIsEllipticSimple);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGeVSimParticle::SetDirectedParam(Float_t v11, Float_t v12, Float_t v13, Float_t v14) {
  //
  // Set parameters for Directed Flow 
  // Actual flow coefficient is calculated as follows
  //
  // V1(Pt,Y) = (V11 + V12*Pt) * sign(Y) * (V13 + V14 * Y^3)
  //
  // where sign = 1 for Y > 0 and -1 for Y < 0
  // 
  // Defaults values
  // v12 = v14 = 0
  // v13 = 1
  // 

  fV1[0] = v11;
  fV1[1] = v12;
  fV1[2] = v13;
  fV1[3] = v14;
  
  fIsDirectedSimple = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGeVSimParticle::SetEllipticParam(Float_t v21, Float_t pTmax, Float_t v22) {
  //
  // Set parameters for Elliptic Flow
  // Actual flow coefficient is calculated as follows
  //
  // pTmax is in GeV 
  // v21 - flow value at saturation
  //
  //
  // V2 = v21 * (pT/pTMax ) * exp (-v22 * y^2)    where pT <= pTmax  
  //      v21 * exp (-v22 * y^2)                   where pT > pTmax  
  //
  // Default values:
  // v22 = 0
  //
  // The parametrisation is suitable for relativistic particles
  // eg. Pions (at RHIC energies)
  //


  fV2[0] = v21;
  fV2[1] = pTmax;
  fV2[2] = v22;

  fIsEllipticSimple = kFALSE;
  fIsEllipticOld = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGeVSimParticle::SetEllipticOld(Float_t v21, Float_t v22, Float_t v23) {
  //
  // Set parameters for Elliptic Flow
  // Actual flow coefficient is calculated as follows
  //
  // V2 = (V21 + V22 pT^2) * exp (-v22 * y^2)
  //
  // The parameterisation is suitable for heavy particles: proton, kaon
  //

  fV2[0] = v21;
  fV2[1] = v22;
  fV2[2] = v23;

  fIsEllipticSimple = kFALSE;
  fIsEllipticOld = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t AliGeVSimParticle::GetDirectedFlow(Float_t pt, Float_t y) {
  //
  // Return coefficient of a directed flow for a given pt and y.
  // For coefficient calculation method refer to SetDirectedParam()
  // 
  
  if (fIsDirectedSimple) return fV1[0];

  Float_t v;
  
  v = (fV1[0] + fV1[1]* pt) * TMath::Sign((Float_t)1.,y) *
    (fV1[2] + fV1[3] * TMath::Abs(y*y*y) );

  return v;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t AliGeVSimParticle::GetEllipticFlow(Float_t pt, Float_t y) {
  //
  // Return coefficient of a elliptic flow for a given pt and y.
  // For coefficient calculation method refer to SetEllipticParam()
  // 

  if (fIsEllipticSimple) return fV2[0];

  if (fIsEllipticOld) {
    
    // old parametrisation
    return (fV2[0]+fV2[1]*pt*pt) * TMath::Exp(-fV2[2]*y*y);

  } else {

    // new "pionic" parameterisation
    if (pt < fV2[1]) return ( (pt / fV2[1]) * fV2[0] * TMath::Exp(-fV2[2]*y*y) ); 
    else  return ( fV2[0] * TMath::Exp(-fV2[2]*y*y) );
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////















