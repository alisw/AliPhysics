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

#include "AliCRThit.h"

ClassImp(AliCRThit)

//____________________________________________________________________________
AliCRThit::AliCRThit()
{
   //
   // default ctor for AliCRThit object
   //

  fId     = 0.;
  fX      = 0.;
  fY      = 0.;
  fZ      = 0.;
  fPx     = 0.;
  fPy     = 0.;
  fPz     = 0.;
  fMedium = 0.;
  fELoss  = 0.;
  fCRTh = 0.;
  fCRTMod = 0.;
  fCRTMag = 0.;
  fCRTRICH = 0.;
  fCRTTPC = 0.;

  fCopy = 0;
  for (Int_t i = 0; i < 5; i++ ) {
    fVolume[i] = 0;
  }

}


//____________________________________________________________________________
AliCRThit::AliCRThit(const AliCRThit & hit)
{
   //
   // copy ctor
   //

  fId     = hit.fId;
  fX      = hit.fX;
  fY      = hit.fY;
  fZ      = hit.fZ;
  fPx     = hit.fPx;
  fPy     = hit.fPy;
  fPz     = hit.fPz;
  fMedium = hit.fMedium;
  fELoss  = hit.fELoss;
  fCRTh = hit.fCRTh;
  fCRTMod = hit.fCRTMod;
  fCRTMag = hit.fCRTMag;
  fCRTRICH = hit.fCRTRICH;
  fCRTTPC = hit.fCRTTPC;

  //fCopy = hit.fCopy;
  //fVolume = hit.fVolume;

}

//_____________________________________________________________________________
AliCRThit& AliCRThit::operator= (const AliCRThit & hit)
{
   //
   // aisngment operator.
   //

  //fId     = hit.fId;
  //fX      = hit.fX;
  //fY      = hit.fY;
  //fZ      = hit.fZ;
  //fPx     = hit.fPx;
  //fPy     = hit.fPy;
  //fPz     = hit.fPz;
  //fMedium = hit.fMedium;
  //fELoss  = hit.fELoss;
  //fCRTh = hit.fCRTh;
  //fCRTMod = hit.fCRTMod;
  //fCRTMag = hit.fCRTMag;
  //fCRTRICH = hit.fCRTRICH;
  //fCRTTPC = hit.fCRTTPC;

  //fCopy = hit.fCopy;
  //fVolume = hit.fVolume;

  return *this;
}

//_____________________________________________________________________________
AliCRThit::AliCRThit(Int_t shunt, Int_t track, Int_t *vol,
                     Float_t *hits) :AliHit(shunt, track)
{
//
// Constructor of hit object
//

  fId     = hits[0];
  fX      = hits[1];
  fY      = hits[2];
  fZ      = hits[3];
  fPx     = hits[4];
  fPy     = hits[5];
  fPz     = hits[6];
  fMedium = hits[7];
  fELoss  = hits[8];
  fCRTh = hits[9];
  fCRTMod = hits[10];
  fCRTMag = hits[11];
  fCRTRICH = hits[12];
  fCRTTPC = hits[13];

  //fTrack = (Int_t)hits[9];

  for (Int_t i = 0; i < 5 ; i++ ) {
    fVolume[i] = vol[i];
  }

}

