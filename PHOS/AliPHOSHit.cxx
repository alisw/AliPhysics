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

//_________________________________________________________________________
// Hit classes for PHOS
//*-- Author : Maxim Volkov, RRC KI
//////////////////////////////////////////////////////////////////////////////

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
AliPHOSHit::AliPHOSHit(Int_t shunt, Int_t track, Int_t id, Float_t *hits):
AliHit(shunt, track)
{

   fId      = id ;
   fX       = hits[0];
   fY       = hits[1];
   fZ       = hits[2];
   fELOS    = hits[3];
}

//____________________________________________________________________________
Bool_t AliPHOSHit::operator==(AliPHOSHit const &rValue) const
{ 
  if ( fId != rValue.GetId() ) return kFALSE;
  
  return kTRUE;
}

//____________________________________________________________________________
AliPHOSHit AliPHOSHit::operator+(const AliPHOSHit &rValue) const
{
  
  AliPHOSHit added(*this);

   added.fX    = rValue.fX  ;
   added.fY    = rValue.fY ;
   added.fZ    = rValue.fZ ;

  added.fELOS += rValue.GetEnergy() ;

  return added;

}

//____________________________________________________________________________
ostream& operator << (ostream& out, const AliPHOSHit& hit) 
{
  out << "AliPHOSHit = " << hit.GetId() << " " << hit.GetEnergy() << endl ;
  return out ;
}



