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
Revision 1.6.4.1  2002/11/26 16:32:24  hristov
Merging NewIO with v3-09-04

Revision 1.6  2002/10/14 14:55:34  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.1.4.3  2002/10/10 14:40:31  hristov
Updating VirtualMC to v3-09-02

Revision 1.5  2002/10/09 15:34:23  gamez
Bad Copy implementation on gcc-3.2 (Yves SCHUTZ)

Revision 1.4  2002/10/09 15:00:17  gamez
Bad operator= implementation on gcc-3.2 (Yves SCHUTZ)

Revision 1.3  2002/10/07 11:19:18  gamez
Changes requested by coding conventions

Revision 1.2  2002/07/25 21:27:22  gamez
Variables renamed to avoid floating exceptions

*/

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

