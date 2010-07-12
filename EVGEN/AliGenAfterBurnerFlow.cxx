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
// Author:
// Sylwester Radomski, 2002
// Martin Poghosyan, 2008
// Constantin Loizides, 2010
//////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TParticle.h>
#include <TLorentzVector.h>
#include <TList.h>
#include <TRandom.h>
#include "AliStack.h"
#include "AliGenAfterBurnerFlow.h"
#include "AliGenCocktailAfterBurner.h"
#include "AliMC.h"
#include "AliCollisionGeometry.h"
#include "AliGenCocktailEntry.h"


ClassImp(AliGenAfterBurnerFlow)


AliGenAfterBurnerFlow::AliGenAfterBurnerFlow():AliGenerator(),
  fReactionPlane(0),
  fHow(0),
  fCounter(0),
  fStack(0)
{
  //
  // Default Construction
  InitPrimaries();
  SetNpParams();
}

AliGenAfterBurnerFlow::AliGenAfterBurnerFlow(Float_t reactionPlane):AliGenerator(),
  fReactionPlane(TMath::Pi()*reactionPlane/180.),
  fHow(1),
  fCounter(0),
  fStack(0)
{
  // reactionPlane - Reaction Plane Angle given in Deg [0-360]
  // but stored and applied in radiants (standard for TParticle & AliCollisionGeometry) 

  InitPrimaries();
  SetNpParams();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

AliGenAfterBurnerFlow::~AliGenAfterBurnerFlow() 
{
  // def. dest.
}

void AliGenAfterBurnerFlow::SetDirectedSimple(Int_t pdg, Float_t v1) 
{
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

void AliGenAfterBurnerFlow::SetDirectedParam(Int_t pdg, Float_t v11, Float_t v12, 
                                             Float_t v13, Float_t v14) 
{
  //
  // Set Directed Flow 
  // Directed flow is parameterised as follows
  //
  // V1(Pt,Y) = (V11 + V12*Pt) * sign(Y) * (V13 + V14 * abs(Y)^3)
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

void AliGenAfterBurnerFlow::SetEllipticSimple(Int_t pdg, Float_t v2) 
{
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

void AliGenAfterBurnerFlow::SetEllipticParam(Int_t pdg, 
                                             Float_t v00, Float_t v10, Float_t v11,
                                             Float_t v22) 
{
  //
  // Set Elliptic Flow
  //
  // Elliptic flow is parametrised to reproduce 
  // V2 of Pions at RHIC energies and is given by:
  // 
  // V2 = (v00 + v10*pt + v11*pt^2) * exp (-v22 * y^2) and zero if V2<0.
  //

  SetFlowParameters(pdg, 2, 3, v00, v10, v11, v22);
}

void AliGenAfterBurnerFlow::SetEllipticParamPion(Int_t pdg, Float_t v21, 
                                                 Float_t pTmax, Float_t v22) 
{
  //
  // Set Elliptic Flow
  //
  // Elliptic flow is parametrised to reproduce 
  // V2 of Pions at RHIC energies and is given by:
  // 
  // V2 = v21 * (pT/pTMax ) * exp (-v22 * y^2)    where pT <= pTmax  
  //      v21 * exp (-v22 * y^2)                  where pT > pTmax  
  //
  // v21   - value at saturation
  // pTmax - saturation transverse momentum
  // v22   - rapidity decreasing
  //

  SetFlowParameters(pdg, 2, 1, v21, pTmax, v22, 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetEllipticParamOld(Int_t pdg, Float_t v21, Float_t v22, Float_t v23) 
{
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

void AliGenAfterBurnerFlow::SetNpParams(Int_t order, Float_t p0, Float_t p1, Float_t p2, Float_t p3)
{
  //
  // Set npart parameterization.
  //

  fNpParams[0] = order;
  fNpParams[1] = p0;
  fNpParams[2] = p1;
  fNpParams[3] = p2;
  fNpParams[4] = p3;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::SetFlowParameters(Int_t pdg, Int_t order, Int_t type, 
                                              Float_t v1, Float_t v2,Float_t v3,Float_t v4) 
{
  // 
  // private function
  // 

  if(TMath::Abs(pdg)>fgkPDG){
    Error("AliAfterBurnerFlow","Overflow");
    return;
  }
  fIsPrim[TMath::Abs(pdg)]=kTRUE;
  
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

void AliGenAfterBurnerFlow::Init() 
{
  // 
  // Standard AliGenerator Initializer
  //

  if(fHow == 0)     { Info("AliGenAfterBurnerFlow", "Using the Hijing R.P. Angle event by event "); }
  else if(fHow == 1){ Info("AliGenAfterBurnerFlow", "Using a fixed R.P. Angle for every event ") ; }
  else { Info("AliGenAfterBurnerFlow", 
              "Using a random R.P. Angle event by event ( ! not the same used by Hijing ! ) "); }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t AliGenAfterBurnerFlow::GetCoefficient(Int_t pdg, Int_t n, Float_t Pt, Float_t Y) const
{
  //
  // private function
  // Return Flow Coefficient for a given particle type flow order
  // and particle momentum (Pt, Y)
  //

  Int_t index = fgkN - n;  // default index (for all pdg)
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

    // New v2 parameterisation
    if (type == 3) {
      v = (fParams[index][3] + fParams[index][4] *Pt + fParams[index][5] *Pt*Pt) *
	TMath::Exp( - fParams[index][6] * Y*Y);
      if (v<0)
        v = 0;
    }
  }
  
  return v;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t AliGenAfterBurnerFlow::GetNpNorm(Int_t npart)
{
  //
  // Calculate npart norm.
  //

  if (npart<0)
    return 1;

  Int_t order = (Int_t)fNpParams[0];
  if (order<0)
    return 1;

  Float_t ret = 0;
  Int_t npp = 1;
  for (Int_t i=0; i<=order; i++) {
    ret += npp*fNpParams[i+1];
    npp *= npart;
  } 
  return ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t AliGenAfterBurnerFlow::IsPrimary(Int_t pdg)
{
  if(pdg>fgkPDG) return kFALSE;
  return fIsPrim[pdg];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t CalcAngle(Double_t phi, Double_t phi0, Double_t phiRP, Double_t v2, Double_t v1=0.)
{
  Double_t phi1 = phi-(phi+2*v1*TMath::Sin(phi-phiRP)+v2*TMath::Sin(2*(phi-phiRP))-phi0)/
    (1.+2*v1*TMath::Cos(phi-phiRP)+ 2*v2*TMath::Cos(2*(phi-phiRP)));
  if(TMath::Abs(phi/phi1-1.)<0.00001) return phi1;
  return CalcAngle(phi1, phi0, phiRP, v2, v1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::InitPrimaries()
{
  for(Int_t i=0; i<fgkPDG; i++) fIsPrim[i]=kFALSE;

  //mesons
  fIsPrim[211]=kTRUE;
  fIsPrim[311]=kTRUE;
  fIsPrim[321]=kTRUE;
  fIsPrim[411]=kTRUE;
  fIsPrim[421]=kTRUE;
  fIsPrim[431]=kTRUE;
  fIsPrim[511]=kTRUE;
  fIsPrim[521]=kTRUE;
  fIsPrim[531]=kTRUE;
  fIsPrim[541]=kTRUE;
  fIsPrim[111]=kTRUE;
  fIsPrim[221]=kTRUE;
  fIsPrim[331]=kTRUE;
  fIsPrim[441]=kTRUE;
  fIsPrim[551]=kTRUE;
  fIsPrim[130]=kTRUE;
  fIsPrim[310]=kTRUE;
  fIsPrim[213]=kTRUE;
  fIsPrim[313]=kTRUE;
  fIsPrim[323]=kTRUE;
  fIsPrim[413]=kTRUE;
  fIsPrim[423]=kTRUE;
  fIsPrim[433]=kTRUE;
  fIsPrim[513]=kTRUE;
  fIsPrim[523]=kTRUE;
  fIsPrim[533]=kTRUE;
  fIsPrim[543]=kTRUE;
  fIsPrim[113]=kTRUE;
  fIsPrim[223]=kTRUE;
  fIsPrim[333]=kTRUE;
  fIsPrim[443]=kTRUE;
  fIsPrim[553]=kTRUE;

  //baryons
  fIsPrim[2112]=kTRUE;
  fIsPrim[2212]=kTRUE;
  fIsPrim[3112]=kTRUE;
  fIsPrim[3122]=kTRUE;
  fIsPrim[3212]=kTRUE;
  fIsPrim[3222]=kTRUE;
  fIsPrim[3312]=kTRUE;
  fIsPrim[3322]=kTRUE;
  fIsPrim[4112]=kTRUE;
  fIsPrim[4122]=kTRUE;
  fIsPrim[4212]=kTRUE;
  fIsPrim[4222]=kTRUE;
  fIsPrim[4132]=kTRUE;
  fIsPrim[4312]=kTRUE;
  fIsPrim[4232]=kTRUE;
  fIsPrim[4322]=kTRUE;
  fIsPrim[4332]=kTRUE;
  fIsPrim[5112]=kTRUE;
  fIsPrim[5122]=kTRUE;
  fIsPrim[5212]=kTRUE;
  fIsPrim[5222]=kTRUE;
  fIsPrim[1114]=kTRUE;
  fIsPrim[2114]=kTRUE;
  fIsPrim[2214]=kTRUE;
  fIsPrim[2224]=kTRUE;
  fIsPrim[3114]=kTRUE;
  fIsPrim[3214]=kTRUE;
  fIsPrim[3224]=kTRUE;
  fIsPrim[3314]=kTRUE;
  fIsPrim[3324]=kTRUE;
  fIsPrim[3334]=kTRUE;
  fIsPrim[4114]=kTRUE;
  fIsPrim[4214]=kTRUE;
  fIsPrim[4224]=kTRUE;
  fIsPrim[4314]=kTRUE;
  fIsPrim[4324]=kTRUE;
  fIsPrim[4334]=kTRUE;
  fIsPrim[5114]=kTRUE;
  fIsPrim[5214]=kTRUE;
  fIsPrim[5224]=kTRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::Generate() 
{
  // 
  // AliGenerator generate function doing actual job.
  // Algorythm:
  //
  // 1. loop over particles on the stack and choose primaries
  // 2. calculate delta phi
  // 3. change phi of primary particle and if it is non-stable 
  //    then its daughters' phi and vertex also
  // 
  // For more details see :
  // M.G. Poghosyan 
  // PWG2 meeting on 06.05.2008 and 03.06.2008


  if (0) 
    for(Int_t ii=0; ii<fCounter;ii++)
    {
      printf("%d  %f  %f  %f  %f\n",ii,fParams[ii][0],fParams[ii][1],fParams[ii][2],fParams[ii][3]);
    } 

  AliGenCocktailAfterBurner *gen;
  
  TParticle *particle;
  TParticle *particleM;
  TLorentzVector momentum;
  TLorentzVector vertex;

  Int_t pdg;
  Float_t phi;
  Float_t pt, y;

  // Get Stack of the first Generator
  //  gen = (AliGenCocktailAfterBurner *)gAlice->Generator();
  gen = (AliGenCocktailAfterBurner *)gAlice->GetMCApp()->Generator();


  AliGenerator* genHijing = 0 ;
  AliCollisionGeometry* geom = 0 ;
  AliGenCocktailEntry* entry = 0 ;
  TList* fEntries = 0 ;

  TRandom* rand = new TRandom(0) ;
  for(Int_t ns=0;ns<gen->GetNumberOfEvents();ns++) 
  {
    gen->SetActiveEventNumber(ns) ;
   
    fStack = gen->GetStack(ns);
    fEntries = gen->Entries() ;

    TIter next(fEntries) ;
    Int_t npart = -1;

    if(fHow == 0) // hijing R.P.
    {
      while((entry = (AliGenCocktailEntry*)next())) 
      {
        Info("Generate (e)","Using R.P. from HIJING ... ");       
        genHijing = entry->Generator() ;        
        if(genHijing->ProvidesCollisionGeometry()) 
        { 
          geom = gen->GetCollisionGeometry(ns) ;
          fReactionPlane = geom->ReactionPlaneAngle() ; 
          npart =  geom->ProjectileParticipants() + geom->TargetParticipants();
          break;
        }
        else 
        {
          Error("Generate (e)", "NO CollisionGeometry !!!  -  using fixed R.P. angle = 0. ") ; 
          fReactionPlane = 0. ; 
        }
      }
    }
    else if(fHow ==1 ) //  fixed R.P. 
    {
      Info("Generate (e)","Using fixed R.P. ...");       
    }
    else
    {
      Info("Generate (e)","Using random R.P.s ... ");       
      fReactionPlane = 2 * TMath::Pi() * rand->Rndm() ;
    }    
    
    cout << " * Reaction Plane Angle (event " << ns << ") = " << fReactionPlane << 
      " rad. ( = " << (360*fReactionPlane/(2*TMath::Pi())) << " deg.) Npart = " << npart << "* " << endl ;

    Int_t nParticles = fStack->GetNprimary();
    for (Int_t i=0; i<nParticles; i++) 
    {
      particle = fStack->Particle(i);
 
      Int_t iM=particle->GetMother(0);
      pdg = particle->GetPdgCode();

      //exclude incoming protons in PYTHIA     
      if(particle->GetPdgCode()==21) continue;

      if(TMath::Abs(pdg)>fgkPDG) continue; 
      // is particle primary?
      if(!fIsPrim[TMath::Abs(pdg)]) continue; 

      if(iM>0)
      {
	particleM = fStack->Particle(iM);
	Int_t pdgM = TMath::Abs(particleM->GetPdgCode());
        // is mother primary?
        if((TMath::Abs(pdgM)<fgkPDG)&&fIsPrim[TMath::Abs(pdgM)]) continue;
      }
   
      particle->Momentum(momentum);
      phi = particle->Phi();

      // get Pt, y    
      pt = momentum.Pt() ; 
      y = 10000.;

      if(TMath::Abs(momentum.Z()) != TMath::Abs(momentum.T())) 
        y = momentum.Rapidity() ; 
      
      Double_t v1 = GetCoefficient(pdg, 1, pt, y);
      Double_t v2 = GetCoefficient(pdg, 2, pt, y);
      Double_t npartnorm = GetNpNorm(npart);
      v2 *= npartnorm;

      //printf("ntup %d %f  %f %f %f %f\n ",npart, v1, v2, pt, y, npartnorm);
     
      Double_t phi1 = CalcAngle(phi, phi, fReactionPlane,v2,v1);
     
      Rotate(i, phi1-phi);
    }
  }

  Info("Generate","Flow After Burner: DONE");  
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void AliGenAfterBurnerFlow::Rotate(Int_t i, Double_t phi, Bool_t IsPrim)
{
  TParticle*  particle = fStack->Particle(i);
  
  TLorentzVector momentum;
  particle->Momentum(momentum);
  momentum.RotateZ(phi);
  particle->SetMomentum(momentum);
 
  if(!IsPrim)
  {
    TLorentzVector vertex;
    particle->ProductionVertex(vertex);
    vertex.RotateZ(phi);
    particle->SetProductionVertex(vertex);
  }
  
  if(particle->GetFirstDaughter()<0) return;     
  for(Int_t iD=particle->GetFirstDaughter(); iD<=particle->GetLastDaughter(); iD++) Rotate(iD, phi, kFALSE);

  return;
}


