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

///////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include "TParticle.h"
#include "TLorentzVector.h"
#include "AliStack.h"
#include "AliGenAfterBurnerFlow.h"
#include "AliMC.h"
#include "AliGenCocktailAfterBurner.h"

// emanuele ---------------------------------------------------------------(
#include <TList.h>
#include "AliCollisionGeometry.h"
#include "AliGenCocktailEntry.h"
#include "TRandom.h"
// emanuele ---------------------------------------------------------------)

ClassImp(AliGenAfterBurnerFlow)

////////////////////////////////////////////////////////////////////////////////////////////////////

    AliGenAfterBurnerFlow::AliGenAfterBurnerFlow(): 
	fReactionPlane(0),
	fCounter(0)
{
    //
    // Default Construction
    //
}

////////////////////////////////////////////////////////////////////////////////////////////////////

AliGenAfterBurnerFlow::AliGenAfterBurnerFlow(Float_t reactionPlane):
	fReactionPlane(reactionPlane),
	fCounter(0) 
{
  //
  // Standard Construction
  // 
  // reactionPlane - Reaction Plane Angle given in Deg [0-360]
  // but stored and applied in radiants (standard for TParticle & AliCollisionGeometry) 

// emanuele ---------------------------------------------------------------(

  if(reactionPlane == 0)     { Info("AliGenAfterBurnerFlow", "Using a random R.P. Angle event by event ( ! not the same used by Hijing ! ) ") ; }
  else if(reactionPlane < 0) { Info("AliGenAfterBurnerFlow", "Using the Hijing R.P. Angle event by event ") ; }
  else if(reactionPlane > 0) { Info("AliGenAfterBurnerFlow", "Using a fixed R.P. Angle ( psi = %d deg.) for every event ", reactionPlane) ; }
  
  // it was // if(reactionPlane < 0 || reactionPlane > 360)   Error("AliGenAfterBurnerFlow", "Reaction Plane Angle - %d - out of bounds [0-360]", reactionPlane); //

// emanuele ---------------------------------------------------------------(

  fReactionPlane = 2 * TMath::Pi() * (reactionPlane/360) ;  // r.p. given in degrees (Radomski's way) but stored and applied in radiants (TParticle's way) 
  fCounter = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

AliGenAfterBurnerFlow::~AliGenAfterBurnerFlow() {
  // Standard Destructor

}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetDirectedSimple(Int_t pdg, Float_t v1) {
  //
  // Set Directed Flow 
  // The same directed flow is applied to all specified particles 
  // independently on transverse momentum or rapidity
  //
  // PDG - particle type to apply directed flow
  //       if (PDG == 0) use as default  
  //

  SetFlowParameters(pdg, 1, 0, v1, 0, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetDirectedParam
(Int_t pdg, Float_t v11, Float_t v12, Float_t v13, Float_t v14) {
  //
  // Set Directed Flow 
  // Directed flow is parameterised as follows
  //
  // V1(Pt,Y) = (V11 + V12*Pt) * sign(Y) * (V13 + V14 * Y^3)
  //
  // where sign = 1 for Y > 0 and -1 for Y < 0
  // 
  // Defaults values
  // v12 = v14 = 0
  // v13 = 1
  //  
  // PDG - particle type to apply directed flow
  //       if (PDG == 0) use as default  
  //
  
  SetFlowParameters(pdg, 1, 1, v11, v12, v13, v14);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetEllipticSimple(Int_t pdg, Float_t v2) {
  //
  // Set Elliptic Flow
  // The same Elliptic flow is applied to all specified particles
  // independently on transverse momentum or rapidity
  //
  // PDG - particle type to apply directed flow
  //       if (PDG == 0) use as default  
  //
  // V2 - flow coefficient
  //      
  // NOTE: for starting playing with FLOW 
  //       start with this function and values 0.05 - 0.1
  //

  SetFlowParameters(pdg, 2, 0, v2, 0, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetEllipticParamPion
(Int_t pdg, Float_t v21, Float_t pTmax, Float_t v22) {
  //
  // Set Elliptic Flow
  //
  // Elliptic flow is parametrised to reproduce 
  // V2 of Pions at RHIC energies and is given by:
  // 
  // V2 = v21 * (pT/pTMax ) * exp (-v22 * y^2)    where pT <= pTmax  
  //      v21 * exp (-v22 * y^2)                   where pT > pTmax  
  //
  // v21   - value at saturation
  // pTmax - saturation transverse momentum
  // v22   - rapidity decrising
  //

  SetFlowParameters(pdg, 2, 1, v21, pTmax, v22, 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetEllipticParamOld
(Int_t pdg, Float_t v21, Float_t v22, Float_t v23) {
  //
  // Set Elliptic Flow
  //
  // Elliptic flow is parameterised using 
  // old MevSim parameterisation
  // 
  // V2 = (V21 + V22 pT^2) * exp (-v22 * y^2)
  //

  SetFlowParameters(pdg, 2, 2, v21, v22, v23, 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetFlowParameters
(Int_t pdg, Int_t order, Int_t type, Float_t v1, Float_t v2,Float_t v3,Float_t v4) {
  // 
  // private function
  // 
  
  Int_t index = 0;
  Bool_t newEntry = kTRUE;

  // Defaults

  if (pdg == 0) {
    index = fgkN - order;
    newEntry = kFALSE;
  }

  // try to find existing entry
  for (Int_t i=0; i<fCounter; i++) {
    if (pdg == (Int_t)fParams[i][0] && 
	order == (Int_t)fParams[i][1]) {
      
      index = i;
      newEntry = kFALSE;
    }
  }
  
  // check fCounter

  if (newEntry && (fCounter > fgkN-3)) {
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
  fParams[index][2] = type;
  fParams[index][3] = v1;
  fParams[index][4] = v2;
  fParams[index][5] = v3;
  fParams[index][6] = v4;  
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::Init() {
  // 
  // Standard AliGenerator Initializer
  //

}

////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t AliGenAfterBurnerFlow::GetCoefficient
(Int_t pdg, Int_t n, Float_t Pt, Float_t Y) {
  //
  // private function
  // Return Flow Coefficient for a given particle type flow order
  // and particle momentum (Pt, Y)
  //

  Int_t index = fgkN - n;  // default index 
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
  
  Int_t type = (Int_t)fParams[index][2];

  if ((Int_t)fParams[index][1] == 1) { // Directed

    if (type == 0 )
      v = fParams[index][3];
    else 
      v = (fParams[index][3] + fParams[index][4] * Pt) * TMath::Sign((Float_t)1.,Y) *
	(fParams[index][5] + fParams[index][6] * TMath::Abs(Y*Y*Y) );

  } else {  // Elliptic

    if (type == 0) v = fParams[index][3];

    // Pion parameterisation 

    if (type == 1) { 
      if (Pt < fParams[index][4]) 
	v = fParams[index][3] * (Pt / fParams[index][4]) ;
      else 
	v = fParams[index][3];
      
      v *= TMath::Exp( - fParams[index][5] * Y * Y);
    }

    // Old parameterisation
    
    if (type == 2) 
      v = (fParams[index][3] + fParams[index][4] * Pt * Pt) *
	TMath::Exp( - fParams[index][5] * Y * Y);
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
  gen = (AliGenCocktailAfterBurner *)gAlice->GetMCApp()->Generator();

// emanuele ---------------------------------------------------------------(

  AliGenerator* genHijing = 0 ;
  AliCollisionGeometry* geom = 0 ;
  AliGenCocktailEntry* entry = 0 ;
  TList* fEntries = 0 ;

  TRandom* rand = new TRandom(0) ;
  Float_t fHow = fReactionPlane ;   // this is a temp. solution not to add a new data member in the .h

  for(Int_t ns=0;ns<gen->GetNumberOfEvents();ns++) 
  {
   gen->SetActiveEventNumber(ns) ;
   stack = gen->GetStack(ns);                                  // it was 0.  
   fEntries = gen->Entries() ;

   TIter next(fEntries) ;
   while((entry = (AliGenCocktailEntry*)next())) 
   {
    if(fHow == 0) // hijing R.P.
    {
     Info("Generate (e)","Using R.P. from HIJING ... ");       
     genHijing = entry->Generator() ;        //  cout <<" * GENERATOR IS "<< genHijing << "  : " << genHijing->GetName() << endl;
     if(genHijing->ProvidesCollisionGeometry()) 
     { 
      geom = gen->GetCollisionGeometry(ns) ; //  cout << " * GEOMETRY YES * " << endl ;
      fReactionPlane = geom->ReactionPlaneAngle() ; 
     }
     else 
     {
      Error("Generate (e)", "NO CollisionGeometry !!!  -  using fixed R.P. angle = 0. ") ; 
      fReactionPlane = 0. ; 
     }
    }
    else if(fHow < 0) // random R.P.
    {
     Info("Generate (e)","Using random R.P.s ... ");       
     fReactionPlane = 2 * TMath::Pi() * rand->Rndm() ;
    }
    else  // if constant R.P. -> do nothing (fReactionPlane already setted)
    {
     Info("Generate (e)","Using a fixed R.P. psi = %d rad.",fReactionPlane);       
    }    
    cout << " * Reaction Plane Angle (event " << ns << ") = " << fReactionPlane << " rad. ( = " << (360*fReactionPlane/(2*TMath::Pi())) << " deg.) * " << endl ;
   }

// emanuele ---------------------------------------------------------------)

   // Loop over particles
   
   for (Int_t i=0; i<stack->GetNtrack(); i++) 
   {
     particle = stack->Particle(i);

     particle->Momentum(momentum);
     pdg = particle->GetPdgCode();
     phi = particle->Phi();

     // get Pt, Y
     
     pt = momentum.Pt() ; 
     //y = momentum.Rapidity() ;

// emanuele ---------------------------------------------------------------(

    if(TMath::Abs(momentum.Z()) != TMath::Abs(momentum.T())) { y = momentum.Rapidity() ; }
    else { y = 0. ; }
    // cout << " * Lorentz Vector (momentum) = " << momentum.X() << " , "  << momentum.Y() << " , " << momentum.Z() << " , " << momentum.T() << " . * " << endl ;
    // cout << " *                        pt = " << momentum.Pt() << " . * " << endl ;
    // cout << " *                         Y = " << y << " . * " << endl ;

// emanuele ---------------------------------------------------------------)

     // Calculate Delta Phi for Directed and Elliptic Flow
     
     dPhi = -2 * GetCoefficient(pdg, 1, pt, y) * TMath::Sin( phi - fReactionPlane );
     dPhi -= GetCoefficient(pdg, 2, pt, y) * TMath::Sin( 2 * (phi - fReactionPlane));
     
     // Set new phi      
     
     phi += dPhi;
     momentum.SetPhi(phi);
     particle->SetMomentum(momentum);
   }

// emanuele ---------------------------------------------------------------(
  }
// emanuele ---------------------------------------------------------------)

  Info("Generate","Flow After Burner: DONE");  
}

////////////////////////////////////////////////////////////////////////////////////////////////////

