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

/*
$Log$

*/

#include "AliCRThit.h"

ClassImp(AliCRThit)

//____________________________________________________________________________
AliCRThit::AliCRThit(const AliCRThit & hit)
{
   //
   // copy ctor for AliCRThit object
   //

  fnmou = hit.fnmou;
  fId   = hit.fId;
  fX    = hit.fX;
  fY    = hit.fY;
  fZ    = hit.fZ;
  fpxug = hit.fpxug;
  fpyug = hit.fpyug;
  fpzug = hit.fpzug;
  flay  = hit.flay;
  fxver = hit.fxver;
  fyver = hit.fyver;
  fzver = hit.fzver;
}
 
//______________________________________________________________________________
AliCRThit::AliCRThit(Int_t shunt, Int_t track, Int_t *vol,
                     Float_t *hits) :AliHit(shunt, track)
{
//
// Constructor of hit object
//

  fnmou = hits[0];
  fId   = hits[1];
  fX    = hits[2];
  fY    = hits[3];
  fZ    = hits[4];
  fpxug = hits[5];
  fpyug = hits[6];
  fpzug = hits[7];
  flay  = hits[8];
  fxver = hits[9];
  fyver = hits[10];
  fzver = hits[11];
}

