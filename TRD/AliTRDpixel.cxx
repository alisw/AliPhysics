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
Revision 1.3  2000/02/28 19:10:26  cblume
Include the new TRD classes

Revision 1.2.4.1  2000/02/28 17:59:27  cblume
Initialize fTrack with -1

Revision 1.2  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Contains the information for one TRD pixel                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDpixel.h"

ClassImp(AliTRDpixel)

//_____________________________________________________________________________
AliTRDpixel::AliTRDpixel():TObject()
{
  //
  // Create a TRD pixel
  // 

  fSignal   =  0;
  fTrack[0] = -1;
  fTrack[1] = -1;
  fTrack[2] = -1;

}

//_____________________________________________________________________________
AliTRDpixel::~AliTRDpixel()
{
  //
  // AliTRDpixel destructor
  //

}

//_____________________________________________________________________________
void AliTRDpixel::Copy(AliTRDpixel &p)
{
  //
  // Copy function
  //

  p.fSignal = fSignal;
  for (Int_t iTrackPixel = 0; iTrackPixel < kTrackPixel; iTrackPixel++) {
    p.fTrack[iTrackPixel] = fTrack[iTrackPixel];
  }

}
