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

#include <TParticle.h>
#include "AliAOD.h"
#include "AliAODParticle.h"
#include "AliTrackPoints.h"

ClassImp(AliAOD)

AliAOD::AliAOD():
 fParticles(10),
 fIsRandomized(kFALSE),
 fPrimaryVertexX(0.0),
 fPrimaryVertexY(0.0),
 fPrimaryVertexZ(0.0)
{
 //ctor
 SetOwner(kTRUE);
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
  AddParticle( new AliAODParticle(*part,idx) );
}
/**************************************************************************/

void  AliAOD::AddParticle(Int_t pdg, Int_t idx,
                          Double_t px, Double_t py, Double_t pz, Double_t etot,
                          Double_t vx, Double_t vy, Double_t vz, Double_t time)
{
  //adds particle to event
  AddParticle(new  AliAODParticle(pdg,idx,px,py,pz,etot,vx,vy,vz,time));
}
/**************************************************************************/

void AliAOD::SwapParticles(Int_t i, Int_t j)
{
//swaps particles positions; used by AliHBTEvent::Blend
  if ( (i<0) || (i>=GetNumberOfParticles()) ) return;
  if ( (j<0) || (j>=GetNumberOfParticles()) ) return;

  AliVAODParticle* tmp = (AliVAODParticle*)fParticles.At(i);
  fParticles.AddAt(fParticles.At(j),i);
  fParticles.AddAt(tmp,j);
}
/**************************************************************************/

void  AliAOD::Reset()
{
  //deletes all particles from the event
   for(Int_t i =0; i<GetNumberOfParticles(); i++)
    {
      for (Int_t j = i+1; j<GetNumberOfParticles(); j++)
        if ( fParticles.At(j) == fParticles.At(i) ) fParticles.RemoveAt(j);
      delete fParticles.RemoveAt(i);
    }
//   fRandomized = kFALSE;
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
  Int_t npart = fParticles.GetEntries();
  for (Int_t i = 0; i < npart; i++)
   {
     AliVAODParticle* p = (AliVAODParticle*)fParticles.At(i);
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

  Int_t npart = fParticles.GetEntries();
  for (Int_t i = 0; i < npart; i++)
   {
     AliVAODParticle* p = (AliVAODParticle*)fParticles.At(i);
     AliTrackPoints* tp  = p->GetTPCTrackPoints();
     if (tp) tp->Move(x,y,z);
     tp  = p->GetITSTrackPoints();
     if (tp) tp->Move(x,y,z);
   }
}
