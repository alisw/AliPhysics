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

//_________________________________________________________________________
//
//      Hit class for VZERO detector   
//  
//_________________________________________________________________________


#include "AliVZEROhit.h"

ClassImp(AliVZEROhit)
 
//_____________________________________________________________________________
AliVZEROhit::AliVZEROhit(Int_t shunt, Int_t track, Int_t* vol, Float_t* hits):
  AliHit(shunt, track)
{
  //
  // Adds a VZERO hit
  //
  
  fVolume          = vol[0];     // Volume ID
  fCopy            = vol[1];     // Copy number
  fX      	   = hits[0];    // X position of hit
  fY               = hits[1];    // Y position of hit
  fZ               = hits[2];    // Z position of hit
  fTrackPiD        = hits[3];    // Track PiD
  fTof             = hits[4];    // Particle time of flight
  fCharge          = hits[5];    // Particle charge
  fTheta           = hits[6];    // Incident theta angle in degrees 
  fPhi             = hits[7];    // Incident phi angle in degrees
  fRingNumber      = hits[8];    // Ring number 
  
  fPt              = hits[9];    // Local transverse momentum of the particle
  fPmom            = hits[10];   // Local P momentum of the particle
  fPx              = hits[11];   // Local Px momentum of the particle
  fPy              = hits[12];   // Local Py momentum of the particle
  fPz              = hits[13];   // Local Pz momentum of the particle
  
  fVx              = hits[14];   // Vertex x coordinate
  fVy              = hits[15];   // Vertex y coordinate
  fVz              = hits[16];   // Vertex z coordinate
  
  fEloss           = hits[17];   // Energy deposited inside volume
  fTleng           = hits[18];   // Track length inside volume
  
}
