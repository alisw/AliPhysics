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
  
  fVolume          = vol[0];
  fCopy            = vol[1];
  fX      	   = hits[0];
  fY               = hits[1];
  fZ               = hits[2];
  fTrackPiD        = hits[3];
  fTof             = hits[4];
  fCharge          = hits[5];
  fTheta           = hits[6];
  fPhi             = hits[7];
  fRingNumber      = hits[8];
  
  fPt              = hits[9];
  fPmom            = hits[10];
  fPx              = hits[11];
  fPy              = hits[12];
  fPz              = hits[13];
  
  fVx              = hits[14];
  fVy              = hits[15];
  fVz              = hits[16];
  
  fEloss           = hits[17];
  fTleng           = hits[18];
  
}
