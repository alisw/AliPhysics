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
//  Contains the information for one TRD pixel                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDpixel.h"

ClassImp(AliTRDpixel)

//_____________________________________________________________________________

  // Maximal number of stored tracks
  const Int_t AliTRDpixel::fgkNTrackPixel = kNTrackPixel;

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
void AliTRDpixel::Copy(TObject &p) const
{
  //
  // Copy function
  //

  ((AliTRDpixel &) p).fSignal = fSignal;
  for (Int_t iTrackPixel = 0; iTrackPixel < kNTrackPixel; iTrackPixel++) {
    ((AliTRDpixel &) p).fTrack[iTrackPixel] = fTrack[iTrackPixel];
  }

}
