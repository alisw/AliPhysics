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


#include "AliVZEROhit.h"

ClassImp(AliVZEROhit)
 
//_____________________________________________________________________________
AliVZEROhit::AliVZEROhit(Int_t shunt, Int_t track, Int_t* vol, Float_t* hits):
  AliHit(shunt, track)
{
  //
  // Add a VZERO hit
  //
  
  fVolume          = vol[2];
  fCopy            = vol[3];
  fX      	   = hits[0];
  fY               = hits[1];
  fZ               = hits[2];
  fXloc            = hits[3];
  fYloc            = hits[4];
  fZloc            = hits[5];
  fEdep            = hits[6];
  fEtot            = hits[7];
  fTrackPiD        = hits[8];
  fParticle        = hits[9];
  fTof             = hits[10];
  fIsTrackEntering = hits[11];
  fIsTrackExiting  = hits[12];
  fCharge          = hits[13];
  fIsCerenkov      = hits[14];
  fMulti           = hits[15];
  fTheta           = hits[16];
  fPhi             = hits[17];
  fNGCerenkovs     = hits[18];
}
