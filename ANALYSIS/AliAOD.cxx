#include "AliAOD.h"
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

/////////////////////////////////////////////////////////////
//
// base class for AOD containers
//
/////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TParticle.h>
#include <TClass.h>
#include <TString.h>
#include "AliAODParticle.h"
#include "AliTrackPoints.h"

ClassImp(AliAOD)

AliAOD::AliAOD():
 fParticles(0x0),
 fIsRandomized(kFALSE),
 fPrimaryVertexX(0.0),
 fPrimaryVertexY(0.0),
 fPrimaryVertexZ(0.0),
 fParticleClass(0x0)
{
 //ctor
// Info("AliAOD()","Entered");
// SetOwner(kTRUE);
// Info("AliAOD()","Exited");
}
/**************************************************************************/

AliAOD::~AliAOD()
{
  //Destructor
  //fParticleClass does not belong to AliAOD -> Do not delete it
  delete fParticles;
  
}
/**************************************************************************/

void AliAOD::SetParticleClassName(const char* classname)
{
//Sets type of particle that is going to be stored 
  if (gROOT == 0x0) Fatal("SetParticleClassName","ROOT System not initialized");
  TClass* pclass = gROOT->GetClass(classname);
  if ( pclass == 0x0 )
   {
     Error("SetParticleClass","Can not get TClass for class named %s",classname);
     return;
   }
  SetParticleClass(pclass);
}
/**************************************************************************/

void AliAOD::SetParticleClass(TClass* pclass)
{
//Sets type of particle that is going to be stored 

  if ( pclass == 0x0 )
   {
     Error("SetParticleClass","Parameter is NULL.");
     return;
   }
   
  if ( pclass->InheritsFrom("AliVAODParticle") == kFALSE )
   {
     Error("SetParticleClass","Class named %s does not inherit from AliVAODParticle",pclass->GetName());
     return;
   }
  if (pclass != fParticleClass)
   {
     fParticleClass = pclass;
     if (fParticleClass) delete fParticles;
     fParticles = new TClonesArray(fParticleClass);
   }
}

/**************************************************************************/

void  AliAOD::AddParticle(TParticle* part, Int_t idx)
{
  //Adds TParticle to event
  if (part == 0x0) 
   {
     Error("AddParticle(TParticle*,Int_t)","pointer to particle is NULL");
     return;
   }

  if (fParticles == 0x0) SetParticleClassName("AliAODParticle");
  AddParticle( new AliAODParticle(*part,idx) );
}
/**************************************************************************/

void  AliAOD::AddParticle(AliVAODParticle* particle)
{
 //add particle to AOD
 //MAKES ITS OWN COPY OF THE PARTICLE!!! (AOD is not going to keep and delete input pointer)
 
  if (fParticles == 0x0) SetParticleClassName("AliAODParticle");

  Int_t idx = fParticles->GetLast() + 1;
  TClonesArray& arr = *fParticles;
  
  AliVAODParticle* pp = (AliVAODParticle*)fParticleClass->New(arr[idx]);
  pp->operator=(*particle);
  
}
/**************************************************************************/

void  AliAOD::AddParticle(Int_t pdg, Int_t idx,
                          Double_t px, Double_t py, Double_t pz, Double_t etot,
                          Double_t vx, Double_t vy, Double_t vz, Double_t time)
{
  //adds particle to event (standard AOD class)

  if (fParticles == 0x0) SetParticleClassName("AliAODParticle");

  Int_t newpartidx = fParticles->GetLast() + 1;
  TClonesArray& arr = *fParticles;
  
  AliVAODParticle* p =  (AliVAODParticle*)fParticleClass->New(arr[newpartidx]);
  
  p->SetPdgCode(pdg);
  p->SetUID(idx);
  p->SetMomentum(px,py,pz,etot);
  p->SetProductionVertex(vx,vy,vz,time);
  
}
/**************************************************************************/

void AliAOD::SwapParticles(Int_t i, Int_t j)
{
//swaps particles positions; used by AliHBTEvent::Blend
  if ( (i<0) || (i>=GetNumberOfParticles()) ) return;
  if ( (j<0) || (j>=GetNumberOfParticles()) ) return;
  
  AliVAODParticle* tmpobj = (AliVAODParticle*)fParticleClass->New();
  AliVAODParticle& tmp = *tmpobj;
  AliVAODParticle& first = *(GetParticle(i));
  AliVAODParticle& second = *(GetParticle(j));
  
  tmp = first;
  first = second;
  second = tmp;
  
}
/**************************************************************************/

void  AliAOD::Reset()
{
  //deletes all particles from the event
   if (fParticles) fParticles->Clear("C");
   
   fIsRandomized = kFALSE;
}
/**************************************************************************/

void AliAOD::GetPrimaryVertex(Double_t&x, Double_t&y, Double_t&z)
{
//returns positions of the primary vertex
  x = fPrimaryVertexX;
  y = fPrimaryVertexY;
  z = fPrimaryVertexZ;
}
/**************************************************************************/

void AliAOD::SetPrimaryVertex(Double_t x, Double_t y, Double_t z)
{
//Sets positions of the primary vertex 
  fPrimaryVertexX = x;
  fPrimaryVertexY = y;
  fPrimaryVertexZ = z;
}
/**************************************************************************/

Int_t AliAOD::GetNumberOfCharged(Double_t etamin, Double_t etamax) const
{
  //reurns number of charged particles within given pseudorapidity range
  Int_t n = 0;
  Int_t npart = GetNumberOfParticles();
  for (Int_t i = 0; i < npart; i++)
   {
     AliVAODParticle* p = GetParticle(i);
     Double_t eta = p->Eta();
     if ( (eta < etamin) || (eta > etamax) ) continue;
     if (p->Charge() != 0.0) n++;
   }
  return n;
}
/**************************************************************************/

void AliAOD::Move(Double_t x, Double_t y, Double_t z)
{
 //moves all spacial coordinates about this vector
 // vertex
 // track points
 // and whatever will be added to AOD and AOD particles that is a space coordinate

  fPrimaryVertexX += x;
  fPrimaryVertexY += y;
  fPrimaryVertexZ += z;

  Int_t npart = GetNumberOfParticles();
  for (Int_t i = 0; i < npart; i++)
   {
     AliVAODParticle* p = GetParticle(i);
     AliTrackPoints* tp  = p->GetTPCTrackPoints();
     if (tp) tp->Move(x,y,z);
     tp  = p->GetITSTrackPoints();
     if (tp) tp->Move(x,y,z);
   }
}

void AliAOD::Print(Option_t* /*option*/)
{
  //Prints AOD
  TString ts;
  TString msg("\n");
  msg+="Particle Class: ";
  if (fParticleClass)
   {
     msg+=fParticleClass->GetName();
   }
  else
   {
     msg+="Not specified yet";
   } 
  msg += "\n";
  msg += "Vertex position X: ";
  msg += fPrimaryVertexX;
  msg += " Y:" ;
  msg += fPrimaryVertexY;
  msg += " Z:";
  msg += fPrimaryVertexZ;
  msg += "\n";
  
  msg += "Randomized: ";
  msg += fIsRandomized;
  msg += "\n";
  
  Info("Print","%s",msg.Data());
  
  Int_t npart = GetNumberOfParticles();
  Info("Print","Npart: %d",npart);
  for (Int_t i = 0; i < npart; i++)
   {
     Info("Print","Getting particle %d",i);
     AliVAODParticle* p = GetParticle(i);
     Info("Print","Printing particle %d, address %#x",i,p);
     p->Dump();
     p->Print();
     Info("Print","particle %d printed",i);
   }
}

void AliAOD::SetOwner(Bool_t /*owner*/)
{
//Sets the ownership of particles: if particles should be also deleted if AOD is deleted/reseted
//Since fParticles is Clones and not Object Array, it is always the owner and this method does not have sense
 
 MayNotUse("SetOwner");
 //if fParticles->SetOwner(owner);
 
}
