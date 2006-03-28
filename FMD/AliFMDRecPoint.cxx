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
/** @file    AliFMDRecPoint.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:46:26 2006
    @brief   Pseudo reconstructed charged particle multiplicity 
*/
//____________________________________________________________________
//
// Base class for reconstructed charged particle multiplicty in the
// FMD.  
// The class contains the field fMethod which is a flag set by
// the AliFMDRecPointAlgorithm that created the object. The flag tells us
// which algorithm was used to create the data stored in the object. 
//
#include "AliFMDRecPoint.h"	// ALIFMDRECPOINT_H
#include <TString.h>            // ROOT_TString 
#include <Riostream.h>		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliFMDRecPoint)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDRecPoint::AliFMDRecPoint()
  : fDetector(0),
    fRing('\0'),
    fSector(0),
    fStrip(0),
    fEta(0),
    fPhi(0),
    fEdep(0)
{}

//____________________________________________________________________
AliFMDRecPoint::AliFMDRecPoint(UShort_t detector,  Char_t   ring, 
			       UShort_t sector,    UShort_t strip,
			       Float_t  eta,       Float_t  phi,
			       Float_t  edep,      Float_t  particles)
  : fDetector(detector),
    fRing(ring),
    fSector(sector),
    fStrip(strip),
    fEta(eta),
    fPhi(phi),
    fEdep(edep),
    fParticles(particles)
{
}

//____________________________________________________________________
const char*
AliFMDRecPoint::GetName() const 
{ 
  static TString n;
  n = Form("FMD%d%c[%2d,%3d]", fDetector,fRing,fSector,fStrip);
  return n.Data();
}

//____________________________________________________________________
const char*
AliFMDRecPoint::GetTitle() const 
{ 
  static TString t;
  t = Form("%f (%f,%f)", fParticles, fEta, fPhi);
  return t.Data();
}

//____________________________________________________________________
void
AliFMDRecPoint::Print(Option_t* option) const
{
  // Print information 
  // 
  // Options:
  //    D:           Detector (default)
  //    E:           Eta range (default)
  //    P:           Phi range (default)
  //
  TString opt(option);
  cout << "FMD RecPoint in a strip: " << fParticles << endl;
  if (opt.Contains("D", TString::kIgnoreCase))
    cout << "  Detector:      FMD" << fDetector << fRing 
	 << "[" << fSector << "," << fStrip << "]" << endl;
  if (opt.Contains("E", TString::kIgnoreCase))
    cout << "  Eta range:     " << fEta << endl;
  if (opt.Contains("P", TString::kIgnoreCase))
    cout << "  Phi range:     " << fPhi << endl;
}

    
//____________________________________________________________________
//
// EOF
//
