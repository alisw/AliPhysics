
////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// AliGenAfterBurnerFlow is a After Burner event generator applying flow.
// The generator changes Phi coordinate of the particle momentum.
// Flow (directed and elliptical) can be defined on particle type level
//
// For examples, parameters and testing macros refer to:
// http:/home.cern.ch/radomski
//
// Author:
// Sylwester Radomski,
// GSI, April 2002
// 
// S.Radomski@gsi.de
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include "TParticle.h"
#include "TLorentzVector.h"
#include "AliStack.h"
#include "AliGenAfterBurnerFlow.h"
#include "AliGenCocktailAfterBurner.h"

ClassImp(AliGenAfterBurnerFlow)

////////////////////////////////////////////////////////////////////////////////////////////////////

AliGenAfterBurnerFlow::AliGenAfterBurnerFlow() {
  // Deafult Construction
  
  fReactionPlane = 0;
  fCounter = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

AliGenAfterBurnerFlow::AliGenAfterBurnerFlow(Float_t reactionPlane) {
  // Standard Construction
  // 
  // reactionPlane - Reaction Plane Angle in Deg

  if (reactionPlane < 0 || reactionPlane > 360)
    Error("AliGenAfterBurnerFlow", 
	  "Reaction Plane Angle - %d - ot of bounds [0-360]", reactionPlane);

  fReactionPlane = reactionPlane;
  fCounter = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

AliGenAfterBurnerFlow::~AliGenAfterBurnerFlow() {
  // Standard Destructor

}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetDirected(Int_t pdg, Float_t v11, Float_t v12, Float_t v13, Float_t v14) {
  //
  // Set Directed Flow parameters for a given particle type.
  // Actual flow coefficient depends on Pt and Y and is caculated by
  //  
  // V1(Pt,Y) = (V11 + V12*Pt) * sign(Y) * (V13 + V14 * Y^3)
  //
  // where sign = 1 for Y > 0 and -1 for Y < 0
  // 
  // Defaults values
  // v12 = v14 = 0
  // v13 = 1
  // 
  // In many cases it is sufficient to set v11 only.
  // Note: be carefull with parameter v14
  // 

  SetFlowParameters(pdg, 1, v11, v12, v13, v14);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetElliptic(Int_t pdg, Float_t v21, Float_t v22, Float_t v23) {
  //
  // Set Elliptic Flow parameters for a given particle type.
  // Actual flow coefficient depends on Pt and Y and is caculated by
  //
  // V2 = (V21 + V22 * Pt^2) * exp( -V23 * Y^2)
  // 
  // Default values:
  // v22 = v23 = 0
  //
  // In many cases it is sufficient to set v21 only
  //

  SetFlowParameters(pdg, 2, v21, v22, v23, 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetDefDirected(Float_t v11, Float_t v12, Float_t v13, Float_t v14) {
  // 
  // Set Directed Flow parameters for all particles.
  // These parameters can be overriden for a specific type by calling
  // SetDirected() function.
  //
  // For explanation of parameters refer to SetDirected()
  //

  SetFlowParameters(0, 1, v11, v12, v13, v14);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetDefElliptic(Float_t v21, Float_t v22, Float_t v23) {
  // 
  // Set Elliptic Flow parameters for all particles.
  // These parameters can be overriden for a specific type by calling
  // SetElliptic() function.
  //
  // For explanation of parameters refer to SetElliptic()
  //

  SetFlowParameters(0, 2, v21, v22, v23, 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetFlowParameters
(Int_t pdg, Int_t order, Float_t v1, Float_t v2,Float_t v3,Float_t v4) {
  // 
  // private function
  // 
  
  Int_t index = 0;
  Int_t newEntry = 1;

  // Defaults

  if (pdg == 0) {
    index = kN - order;
    newEntry = 0;
  }

  // try to find existing entry
  for (Int_t i=0; i<fCounter; i++) {
    if (pdg == (Int_t)fParams[i][0] && 
	order == (Int_t)fParams[i][1]) {
      
      index = i;
      newEntry = 0;
    }
  }
  
  // check fCounter

  if (newEntry && (fCounter > kN-3)) {
    Error("AliAfterBurnerFlow","Overflow");
    return;
  }
  
  if (newEntry) {
    index = fCounter;
    fCounter++;
  }
  
  // Set new particle type

  fParams[index][0] = pdg;
  fParams[index][1] = order;
  fParams[index][2] = v1;
  fParams[index][3] = v2;
  fParams[index][4] = v3;
  fParams[index][5] = v4;  
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::Init() {
  // 
  // Standard AliGenerator Initializer
  //

}

////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t AliGenAfterBurnerFlow::GetCoeff
(Int_t pdg, Int_t n, Float_t Pt, Float_t Y) {
  //
  // private function
  // Return Flow Coefficient for a given particle type flow order
  // and particle momentum (Pt, Y)
  //

  Int_t index = kN - n;  // default index 
  Float_t v = 0;

  // try to find specific parametrs

  for (Int_t i=0; i<fCounter; i++) {
    
    if ((Int_t)fParams[i][0] == pdg &&
	(Int_t)fParams[i][1] == n) {
      
      index = i;
      break;
    }
  } 
  
  // calculate v
  
  if ((Int_t)fParams[index][1] == 1) { // Directed
    
    v = (fParams[index][2] + fParams[index][3] * Pt) * TMath::Sign((Float_t)1.,Y) *
      (fParams[index][4] + fParams[index][5] * TMath::Abs(Y*Y*Y) );

  } else {  // Elliptic

    v = (fParams[index][2] + fParams[index][3] * Pt * Pt) *
      TMath::Exp( - fParams[index][4] * Y * Y);
  }
  
  return v;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::Generate() {
  // 
  // AliGenerator generate function doing actual job.
  // Algorythm:
  //
  // 1. loop over particles on the stack
  // 2. find direct and elliptical flow coefficients for 
  //    a particle type ore use defaults
  // 3. calculate delta phi
  // 4. change phi in orginal particle
  // 
  // Algorythm based on:
  // A.M. Poskanzer, S.A. Voloshin
  // "Methods of analysisng anisotropic flow in relativistic nuclear collisions"
  // PRC 58, 1671 (September 1998)
  //
  
  AliGenCocktailAfterBurner *gen;
  AliStack *stack;
  TParticle *particle;
  TLorentzVector momentum;

  Int_t pdg;
  Float_t phi, dPhi;
  Float_t pt, y;

  // Get Stack of the first Generator
  gen = (AliGenCocktailAfterBurner *)gAlice->Generator();
  stack = gen->GetStack(0);

  // Loop over particles

  for (Int_t i=0; i<stack->GetNtrack(); i++) {

    particle = stack->Particle(i);

    particle->Momentum(momentum);
    pdg = particle->GetPdgCode();
    phi = particle->Phi();

    // get Pt, Y

    pt = momentum.Pt();
    y = momentum.Rapidity();

    // Calculate Delta Phi for Directed and Elliptic Flow
    
    dPhi = -2 * GetCoeff(pdg, 1, pt, y) * TMath::Sin( phi - fReactionPlane );
    dPhi -= GetCoeff(pdg, 2, pt, y) * TMath::Sin( 2 * (phi - fReactionPlane));
    
    
    // cout << i << "\t" << pt << "\t" << y << "\t" << (GetCoeff(pdg, 1, pt, y)) << "\t"
    //	 << (GetCoeff(pdg, 2, pt, y)) << "\t" << dPhi << endl;

    // Set new phi 
    
    phi += dPhi;
    momentum.SetPhi(phi);
    particle->SetMomentum(momentum);
  }

  cout << "Flow After Burner: DONE" << endl;
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////

