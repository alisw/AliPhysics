//_________________________________________________________________________
//  TOF hit  : member variables
//  fTrack   :
//  fX       : X coordinate of the hit in the Master Reference Frame (LAB Frame)
//  fY       : Y coordinate of the hit in the Master Reference Frame (LAB Frame)
//  fZ       : Z coordinate of the hit in the Master Reference Frame (LAB Frame)
//  fSector  : Number of the TOF Sector which belongs the hit 
//  fPlate   : Number of the TOF Plate or Module which belongs the hit 
//  fStrip   : Number of the TOF Strip which belongs the hit 
//  fPadx    : Number of the pad in the strip along the x-axis - in the strip reference frame
//             - where hit is produced 
//  fPadz    : Number of the pad in the strip along the z-axis - in the strip reference frame
//             - where hit is produced
//  fPx      : x-director cosine of the Charged Particle Momentum when hit is
//             produced - expressed in the Master Reference Frame (LAB Frame) -
//  fPy      : y-director cosine of the Charged Particle Momentum when hit is
//             produced - expressed in the Master Reference Frame (LAB Frame) -
//  fPz      : z-director cosine of the Charged Particle Momentum when hit is
//             produced - expressed in the Master Reference Frame (LAB Frame) -
//  fPmom    : Modulus of the Charged Particle Momentum when hit is produced
//  fTof     : Time of Flight i.e. the time between the charged particle is produced and this
//             particle produce the hit on the TOF sensible volume (pad)
//  fDx      : Distance of the hit from the pad edge along x-axis
//  fDy      : y coordinate of the hit in the pad refernce frame  
//  fDz      : Distance of the hit from the pad edge along z-axis
//  fIncA    : Incidence Angle between the Normal to the sensible volume where hit
//             is produced (pad) and the Momentum Direction of the Charged Particle which
//             produces the hit
//  fEdep    : Energy released by charged particle on the sensible TOF volume where hit is
//             produced
// For more detailed informations about the meaning of the TOF-hit member
// variable look at 
// http://www.bo.infn.it/alice/alice-doc/TOFWEB/variables-hits.html
//
//  Getters, setters and member functions  defined here
//
//*-- Authors: F. Pierella, A. Seganti, D. Vicinanza



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

#include "AliTOFhit.h"

ClassImp(AliTOFhit)

//____________________________________________________________________________
  AliTOFhit::AliTOFhit():
  AliHit(),
  fSector(-1),
  fPlate(-1),
  fStrip(-1),
  fPadx(-1),
  fPadz(-1),
  fPx(0),
  fPy(0),
  fPz(0),
  fPmom(0),
  fTof(0),
  fDx(0),
  fDy(0),
  fDz(0),
  fIncA(0),
  fEdep(0)
{
}

//____________________________________________________________________________
AliTOFhit::AliTOFhit(const AliTOFhit & hit)
  : AliHit(hit),
  fSector(hit.fSector),
  fPlate(hit.fPlate),
  fStrip(hit.fStrip),
  fPadx(hit.fPadx),
  fPadz(hit.fPadz),
  fPx(hit.fPx),
  fPy(hit.fPy),
  fPz(hit.fPz),
  fPmom(hit.fPmom),
  fTof(hit.fTof),
  fDx(hit.fDx),
  fDy(hit.fDy),
  fDz(hit.fDz),
  fIncA(hit.fIncA),
  fEdep(hit.fEdep)
{
   //
   // copy ctor for AliTOFhit object
   //
  fTrack = hit.fTrack;
}
 
//______________________________________________________________________________
AliTOFhit::AliTOFhit(Int_t shunt, Int_t track, Int_t *vol,
                     Float_t *hits)
  :AliHit(shunt, track),
  fSector(-1),
  fPlate(-1),
  fStrip(-1),
  fPadx(-1),
  fPadz(-1),
  fPx(0),
  fPy(0),
  fPz(0),
  fPmom(0),
  fTof(0),
  fDx(0),
  fDy(0),
  fDz(0),
  fIncA(0),
  fEdep(0)
{
//
// Constructor of hit object
//
  //
  // Hit Volume
  // 
  fSector= vol[0];
  fPlate = vol[1];
  fStrip = vol[2];
  fPadx = vol[3];
  fPadz = vol[4];
  //
  //Position of the hit
  fX = hits[0];
  fY = hits[1];
  fZ = hits[2];
  //
  // Momentum components of the particle in the ALICE frame when hit is produced
  fPx  = hits[3];
  fPy  = hits[4];
  fPz  = hits[5];
  fPmom= hits[6];
  //
  // Time Of Flight for the particle that produces hit
  fTof = hits[7];   //TOF[s]
  //
  // Other Data
  fDx  = hits[8];   //Distance from the edge along x axis
  fDy  = hits[9];   //Y cohordinate of the hit
  fDz  = hits[10];  //Distance from the edge along z axis
  fIncA= hits[11];  //Incidence angle
  fEdep= hits[12];  //Energy loss in TOF pad
}

