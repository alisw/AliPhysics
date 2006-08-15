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
//                                                                           //
//  Hit object for the TRD                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDhit.h"

ClassImp(AliTRDhit)

//_____________________________________________________________________________
AliTRDhit::AliTRDhit()
  :AliHit()
  ,fDetector(0)
  ,fQ(0)
{
  //
  // AliTRDhit default constructor
  //

}

//_____________________________________________________________________________
AliTRDhit::AliTRDhit(Int_t shunt, Int_t track, Int_t det
                   , Float_t *hits, Int_t q)
  :AliHit(shunt,track)
  ,fDetector((UShort_t) det)
  ,fQ((Short_t) q)
{
  //
  // Create a TRD hit
  //

  // Store position 
  fX = hits[0];
  fY = hits[1];
  fZ = hits[2];

}

//_____________________________________________________________________________
AliTRDhit::~AliTRDhit()
{
  //
  // AliTRDhit destructor
  //

}
