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
//               
//*-- Author: Maxime Volkov (RRC KI) & Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- Standard library ---
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strstream.h>

// --- AliRoot header files ---
#include "AliPHOSHit.h"
#include "AliRun.h"
#include "AliConst.h"


ClassImp(AliPHOSHit)

//____________________________________________________________________________
AliPHOSHit::AliPHOSHit(const AliPHOSHit & hit) 
{
   // copy ctor
   
  fId      = hit.fId ; 
  fELOS    = hit.fELOS ;
  fPrimary = hit.fPrimary ; 
  fTrack   = hit.fTrack ; 
  fMomentum= hit.fMomentum ;
  fX       = hit.fX ; 
  fY       = hit.fY ; 
  fZ       = hit.fZ ; 
  fPid     = hit.fPid ;
 
} 

//____________________________________________________________________________
AliPHOSHit::AliPHOSHit(Int_t shunt, Int_t primary, Int_t track, Int_t id, Float_t *hits, Int_t pid, TLorentzVector p, Float_t *xy): AliHit(shunt, track)
{
  //
  // Create a CPV hit object
  //

  fId         = id ;
  fELOS       = hits[3] ;
  fPrimary    = primary ;
  fPid        = pid ; 
  fMomentum  = p;
  fX         = xy[0];  //position of particle first entering cristall/pad
  fY         = xy[1];
  fZ         = xy[2];  
}

//____________________________________________________________________________
Bool_t AliPHOSHit::operator==(AliPHOSHit const &rValue) const
{ 
  // Two hits are identical if they have the same Id and originate from the same primary 

  Bool_t rv = kFALSE ; 

  if ( (fId == rValue.GetId()) && ( fPrimary == rValue.GetPrimary() ) && (fPid*rValue.fPid == 0) )
    rv = kTRUE;
  
  return rv;
}

//____________________________________________________________________________
AliPHOSHit AliPHOSHit::operator+(const AliPHOSHit &rValue)
{
  // Add the energy of the hit
  
  fELOS += rValue.GetEnergy() ;
    
  if((fPid == 0) && (rValue.fPid != 0)){
    fPid = rValue.fPid ;
    fX    = rValue.fX ;
    fY    = rValue.fY ;
    fZ    = rValue.fZ ;
   }     

   return *this;

}

//____________________________________________________________________________
ostream& operator << (ostream& out, const AliPHOSHit& hit) 
{
  // Print out Id and energy 
  
  out << "AliPHOSHit = " << hit.GetId() << " " << hit.GetEnergy() << endl ;
  return out ;
}



