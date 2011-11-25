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

//_________________________________________________________________________
//  Hits class for PHOS     
//  A hit in PHOS is the sum of all hits in a single crystal
//*--
//*-- Author: Maxime Volkov (RRC KI) & Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSHit.h"

ClassImp(AliPHOSHit)
  
  //____________________________________________________________________________
AliPHOSHit::AliPHOSHit(const AliPHOSHit & hit) :
  AliHit(hit),fId(hit.fId),fELOS(hit.fELOS),fTime(hit.fTime)
{
  // copy ctor
} 

//____________________________________________________________________________
AliPHOSHit::AliPHOSHit(Int_t shunt, Int_t track, Int_t id, Float_t *hits) : 
  AliHit(shunt, track),fId(0),fELOS(0),fTime(0)
{
  //
  // Create a CPV hit object
  //
  
  fX          = hits[0] ;
  fY          = hits[1] ;
  fZ          = hits[2] ;
  fTime       = hits[3] ;
  fId         = id ;
  fELOS       = hits[4] ;
}
//____________________________________________________________________________
Bool_t AliPHOSHit::operator==(AliPHOSHit const &rValue) const
{ 
  // Two hits are identical if they have the same Id and originate from the same primary 

  Bool_t rv = kFALSE ; 

  if ( (fId == rValue.GetId()) && ( fTrack == rValue.GetPrimary() ) )
    rv = kTRUE;
  
  return rv;
}

//____________________________________________________________________________
AliPHOSHit & AliPHOSHit::operator = (const AliPHOSHit &)
{
  Fatal("operator =", "not implemented");
  return *this;
}
//____________________________________________________________________________
AliPHOSHit AliPHOSHit::operator+(const AliPHOSHit &rValue)
{
  // Add the energy of the hit
  
  fELOS += rValue.GetEnergy() ;

  if(rValue.GetTime() < fTime)
    fTime = rValue.GetTime() ;
    
   return *this;

}

