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

// **************************************************************
//
//  		Hits classes for ZDC                  
//
// **************************************************************

#include "AliZDCHit.h"

ClassImp(AliZDCHit)
  
//_____________________________________________________________________________
AliZDCHit::AliZDCHit() :
//  AliHit(shunt, track),
  fPrimKinEn(0.),
  fXImpact(0.),
  fYImpact(0.),
  fSFlag(0),
  fLightPMQ(0.),
  fLightPMC(0.),
  fEnergy(0.), 
  fPDGCode(0),
  fTrackTOF(0.)

{
  //
  // Default constructor
  //
  for(Int_t i=0; i<2; i++) fVolume[i] = 0;
}

//_____________________________________________________________________________
AliZDCHit::AliZDCHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits) :
  AliHit(shunt, track),
  fPrimKinEn(hits[3]),
  fXImpact(hits[4]),
  fYImpact(hits[5]),
  fSFlag(hits[6]),
  fLightPMQ(hits[7]),
  fLightPMC(hits[8]),
  fEnergy(hits[9]), 
  fPDGCode((Int_t) hits[10]),
  fTrackTOF(hits[11])

{
  //
  // Standard constructor
  //
  Int_t i;
  for(i=0; i<2; i++) fVolume[i] = vol[i];
  fX 		= hits[0];
  fY 		= hits[1];
  fZ 		= hits[2];
}
  
//_____________________________________________________________________________
AliZDCHit::AliZDCHit(const AliZDCHit &oldhit) :
  AliHit(0,oldhit.GetTrack()),
  fPrimKinEn(oldhit.GetPrimKinEn()),
  fXImpact(oldhit.GetXImpact()),  
  fYImpact(oldhit.GetYImpact()),  
  fSFlag(oldhit.GetSFlag()),
  fLightPMQ(oldhit.GetLightPMQ()), 
  fLightPMC(oldhit.GetLightPMC()),
  fEnergy(oldhit.GetEnergy()),
  fPDGCode(oldhit.GetPDGCode()),
  fTrackTOF(oldhit.GetTrackTOF())
{
  // Copy constructor
  fX = oldhit.X();
  fY = oldhit.Y();
  fZ = oldhit.Z();
  for(Int_t i=0; i<2; i++) fVolume[i] = oldhit.GetVolume(i);
}
  
  
//_____________________________________________________________________________
void AliZDCHit::Print(Option_t *) const 
{
   // Print method
   printf("\t AliZDCHit: track %d PDGcode %d TOF %f ns E_prim = %f GeV SFlag = %1.0f\n" 
	  "\t det =  %d tow =  %d  (X, Y)impact (%f, %f) cm\n"
          "\t PMQLight %1.0f, PMCLight %1.0f,  E_dep %f\n ", 
          fTrack,fPDGCode,fTrackTOF,fPrimKinEn,fSFlag,
	  fVolume[0],fVolume[1],fXImpact,fYImpact,
          fLightPMQ,fLightPMC,fEnergy);
}
