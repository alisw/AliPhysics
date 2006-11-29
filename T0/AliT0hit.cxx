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
// AliT0hit is the hit class for the T0. Hits are the information
// that comes from a Monte Carlo at each step as a particle mass through
// sensitive detector elements as particles are transported through a
// detector.
//
// Data members:
//
// Int_t fTrack
//     See AliHit for a full description. The track number of the track
// that made this hit.
//
// Float_t fX
//     See AliHit for a full description. The global x position of the
// hit (in the standard units of the Monte Carlo).
//
// Float_t fY
//     See AliHit for a full description. The global y position of the
// hit (in the standard units of the Monte Carlo).
//
// Float_t fZ
//     See AliHit for a full description. The global z position of the
// hit (in the standard units of the Monte Carlo).
//
// Int_t fStatus
//     The track status flag. This flag indicates the track status
// at the time of creating this hit. It is made up of the following 8
// status bits from highest order to lowest order bits
// 0           :  IsTrackAlive():    IsTrackStop():IsTrackDisappeared():
// IsTrackOut():IsTrackExiting():IsTrackEntering():IsTrackInside()     .
// See AliMC for a description of these functions. If the function is
// true then the bit is set to one, otherwise it is zero.
//
// Int_t fVolume
//     The number of the T0 detector that contains this hit.
//     0 - right array; 1 - left array 
// Int_t fPmt 
// the number of PMT tube that contains hit
// Float_t fEdep
//     The energy lost by the particle during the step ending in this
// hit. The units are those determined by the Monte Carlo.
//
// Float_t fTime
//     The time of flight associated with the particle  in this
// hit. The time is typically measured from the point of creation of the
// original particle (if this particle is a daughter).  The units
// are those determined by the Monte Carlo.
///////////////////////////////////////////////////////////////////////
  

#include "AliT0hit.h"

ClassImp(AliT0hit)

AliT0hit::AliT0hit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
//Normal T0 hit ctor
  
   fVolume = vol[0];
   fPmt=vol[1];
   fX=hits[0];
   fY=hits[1];
   fZ=hits[2];
   fEtot=Double_t (hits[3]);
   fParticle=Int_t (hits[4]);
   fTime=hits[5];
}
 
