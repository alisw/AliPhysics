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

#include "AliAODStdParticle.h"

ClassImp(AliAOD)

/**************************************************************************/

void  AliAOD::AddParticle(TParticle* part, Int_t idx)
{
  //Adds TParticle to event
  if (part == 0x0) 
   {
     Error("AddParticle(TParticle*,Int_t)","pointer to particle is NULL");
     return;
   }
  AddParticle( new AliAODStdParticle(*part,idx) );
}
/**************************************************************************/

void  AliAOD::AddParticle(Int_t pdg, Int_t idx,
                          Double_t px, Double_t py, Double_t pz, Double_t etot,
                          Double_t vx, Double_t vy, Double_t vz, Double_t time)
{
  //adds particle to event
  AddParticle(new  AliAODStdParticle(pdg,idx,px,py,pz,etot,vx,vy,vz,time));
}
/**************************************************************************/

void AliAOD::SwapParticles(Int_t i, Int_t j)
{
//swaps particles positions; used by AliHBTEvent::Blend
  if ( (i<0) || (i>=GetNumberOfParticles()) ) return;
  if ( (j<0) || (j>=GetNumberOfParticles()) ) return;

  AliAODParticle* tmp = (AliAODParticle*)fParticles.At(i);
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
