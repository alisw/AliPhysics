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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.6  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  Hit impact class for PHOS     
//  A hit impact in PHOS is a set of parameters of a track which
//  enters the detector for the first time.
//  Track parameters are:
//  - track number
//  - primary particle number
//  - type of a particle
//  - impact coordinate
//  - impact 4-momentum
//               
//*-- Author: Yuri Kharlov (IHEP, Protvino/SUBATECH, Nantes)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSImpact.h"

ClassImp(AliPHOSImpact)

//____________________________________________________________________________
AliPHOSImpact::AliPHOSImpact() : 
  fPid(0),
  fPrimary(0),
  fMomentum()
{
  //Def ctor.
}

//____________________________________________________________________________
AliPHOSImpact::AliPHOSImpact(const AliPHOSImpact & hit) : 
  AliHit(hit),
  fPid(hit.fPid),
  fPrimary(hit.fPrimary),
  fMomentum(hit.fMomentum)
{
  // copy ctor
  fTrack   = hit.fTrack ; 
  fX       = hit.fX ; 
  fY       = hit.fY ; 
  fZ       = hit.fZ ; 
} 

//____________________________________________________________________________
AliPHOSImpact::AliPHOSImpact(Int_t shunt, Int_t primary, Int_t track, Int_t pid, TLorentzVector p, Float_t *xyz): 
  AliHit(shunt, track),
  fPid(pid),
  fPrimary(primary),
  fMomentum(p)
{
  //
  // Create a PHOS hit impact object
  //
  fX         = xyz[0];  //position of particle first entering cristall/pad
  fY         = xyz[1];
  fZ         = xyz[2];  
}

//____________________________________________________________________________
AliPHOSImpact & AliPHOSImpact::operator = (const AliPHOSImpact &)
{
  Fatal("operator =", "not implemented");
  return *this;
}

//____________________________________________________________________________
void AliPHOSImpact::Print(const Option_t *)const
{
  // Prints particle info
  printf("Impact Pid=%4d, p=(%6.2f,%6.2f,%6.2f,%6.2f) GeV, x=(%6.1f,%6.1f,%6.1f) cm\n",
	 fPid,fMomentum.Px(),fMomentum.Py(),fMomentum.Pz(),fMomentum.E(),
	 fX,fY,fZ);
}
