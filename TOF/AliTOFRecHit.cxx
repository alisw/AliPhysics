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


////////////////////////////////////////////////////////////////////////
//  Dummy hit for TOF reconstruction : member variables description
//
//  fTrack   :      track number of the particle that produced the hit
//  fPdgCode :      GEANT code of the particle that produced the hit
//  fX       :      x-coordinate of the hit 
//  fY       :      y-coordinate of the hit 
//  fZ       :      z-coordinate of the hit
//  fP       :      momentum
//  fVrho    :      rho-coordinate of the Vertex
//  fFirst   :      =1 for the first hit of the track, =0 otherwise
//  fNoise   :      =1 for the noise hit (Rvtx>200 or second, ... hit), 
//                  =0 otherwise
//  fRmin    :      distance to the nearest TOFhit
//
// For more detailed informations about the meaning of the hit 
// for TOF reconstruction member variable look at 
// http://bogrid1.bo.infn.it/~pierella/TOFWEB/index.php3
//
// -- Authors: Bologna-ITEP-Salerno Group
//
// Description: dummy hit class used in reconstruction (derived from AliHit)
// For a given TOF hit, the class contains:
// - the distance to the nearest hit
// - flag for first or second track crossing
// - number of the track which produced the hit
// - flag for noise
////////////////////////////////////////////////////////////////////////////

#include "AliTOFRecHit.h"

ClassImp(AliTOFRecHit)

//____________________________________________________________________________
AliTOFRecHit::AliTOFRecHit(const AliTOFRecHit & hit)
:AliHit(hit)
{
  //
  // copy ctor for AliTOFRecHit object
  //
  fTrack  = hit.fTrack;
  fPdgCode= hit.fPdgCode;
  fX      = hit.fX;
  fY      = hit.fY;
  fZ      = hit.fZ;
  fP      = hit.fP;
  fVrho   = hit.fVrho;
  fFirst  = hit.fFirst; 
  fNoise  = hit.fNoise;
  fRmin   = hit.fRmin;

}
 
//______________________________________________________________________________
AliTOFRecHit::AliTOFRecHit(Int_t shunt, Int_t track)
:AliHit(shunt, track)
{
  //
  // ctor for hit object
  //
  fTrack=0;
  fPdgCode=0;
  fX=0;
  fY=0;
  fZ=0;
  fP=-1;
  fVrho=-1;
  fFirst=0;
  fNoise=0;
  fRmin=-1;
}

//______________________________________________________________________________
void AliTOFRecHit::SetHit(Int_t track, Int_t pdgCode, Float_t* mrfpos, Float_t mom, Float_t vtxRadius, Int_t isFirstHit)
{
  // Setter for
  // track number, pdg code, hit position in master reference frame, 
  // momentum, vertex radius and flag to check if it is the first hit
  //
  fTrack  =track;
  fPdgCode=pdgCode;
  fX=mrfpos[0];
  fY=mrfpos[1];
  fZ=mrfpos[2];
  fP=mom;
  fVrho=vtxRadius;
  fFirst=isFirstHit;
}
