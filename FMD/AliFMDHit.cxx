/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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

//____________________________________________________________________
//
//  Hits in the FMD 
//
// Latest changes by Christian Holm Christensen
//
#include "AliFMDHit.h"		// ALIFMDHIT_H
#include "AliLog.h"		// ALILOG_H
#include "Riostream.h"		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliFMDHit)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif


//____________________________________________________________________
AliFMDHit::AliFMDHit()
  : fDetector(0), 
    fRing(0), 
    fSector(0), 
    fStrip('\0'), 
    fPx(0),
    fPy(0),
    fPz(0),
    fPdg(0),
    fEdep(0), 
    fTime(0)
{
  fX = fY = fZ = 0;
}
  

//____________________________________________________________________
AliFMDHit::AliFMDHit(Int_t    shunt, 
		     Int_t    track, 
		     UShort_t detector, 
		     Char_t   ring, 
		     UShort_t sector, 
		     UShort_t strip, 
		     Float_t  x, 
		     Float_t  y, 
		     Float_t  z,
		     Float_t  px, 
		     Float_t  py, 
		     Float_t  pz,
		     Float_t  edep,
		     Int_t    pdg,
		     Float_t  t)
  : AliHit(shunt, track),
    fDetector(detector), 
    fRing(ring), 
    fSector(sector), 
    fStrip(strip), 
    fPx(px),
    fPy(py),
    fPz(pz),
    fPdg(pdg),
    fEdep(edep), 
    fTime(t)
{
  // Normal FMD hit ctor
  // 
  // Parameters:
  // 
  //    shunt     ???
  //    track	  Track #
  //    detector  Detector # (1, 2, or 3)                      
  //    ring	  Ring ID ('I' or 'O')
  //    sector	  Sector # (For inner/outer rings: 0-19/0-39)
  //    strip	  Strip # (For inner/outer rings: 0-511/0-255)
  //    x	  Track's X-coordinate at hit
  //    y	  Track's Y-coordinate at hit
  //    z	  Track's Z-coordinate at hit
  //    px	  X-component of track's momentum 
  //    py	  Y-component of track's momentum
  //    pz	  Z-component of track's momentum
  //    edep	  Energy deposited by track
  //    pdg	  Track's particle Id #
  //    t	  Time when the track hit 
  // 
  fX = x;
  fY = y;
  fZ = z;
}

//____________________________________________________________________
void
AliFMDHit::Print(Option_t* /* option */) const 
{
  // Print Hit to standard out 
  cout << "AliFMDHit: FMD" 
       << fDetector << fRing << "[" 
       << setw(3) << fSector << ","
       << setw(3) << fStrip << "] = " 
       << fEdep << endl;
}

//____________________________________________________________________
//
// EOF
//
