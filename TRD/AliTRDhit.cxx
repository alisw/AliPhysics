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
Revision 1.3  2000/06/07 16:25:37  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.2.1  2000/05/08 14:48:31  cblume
AliTRDhit class now in separate files

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Hit object for the TRD                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDhit.h"

ClassImp(AliTRDhit)

//_____________________________________________________________________________
AliTRDhit::AliTRDhit():AliHit()
{
  //
  // AliTRDhit default constructor
  //

}

//_____________________________________________________________________________
AliTRDhit::AliTRDhit(Int_t shunt, Int_t track, Int_t *det, Float_t *hits)
          :AliHit(shunt, track)
{
  //
  // Create a TRD hit
  //

  // Store volume hierarchy
  fDetector = det[0];

  // Store position and charge
  fX        = hits[0];
  fY        = hits[1];
  fZ        = hits[2];
  fQ        = hits[3];

}

//_____________________________________________________________________________
AliTRDhit::~AliTRDhit()
{
  //
  // AliTRDhit destructor
  //

}
