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
//  FMD reconstructed multiplicity in a strip.   It contains the
//  information of which strip in what sector of what ring, in which
//  detector the information belongs, as well as the pseudo-rapidity,
//  and azimuthal angle the strip had in the event.  Note, that this
//  may change when the interaction points z--coordinate changes
//  (which it probably will - experience from RHIC says so).  Also,
//  the reconstructed energy deposited is stored. 
//
//
#include "AliFMDMultStrip.h"	// ALIFMDMULTSTRIP_H
#include <TString.h>            // ROOT_TString
#include <Riostream.h>		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliFMDMultStrip)


//____________________________________________________________________
AliFMDMultStrip::AliFMDMultStrip()
  : fDetector(0),
    fRing('\0'),
    fSector(0),
    fStrip(0),
    fEta(0),
    fPhi(0),
    fEdep(0)
{}

//____________________________________________________________________
AliFMDMultStrip::AliFMDMultStrip(UShort_t detector,  Char_t   ring, 
				 UShort_t sector,    UShort_t strip,
				 Float_t  eta,       Float_t  phi,
				 Float_t  edep,      Float_t  particles, 
				 UShort_t method)
  : AliFMDMult(particles, method),
    fDetector(detector),
    fRing(ring),
    fSector(sector),
    fStrip(strip),
    fEta(eta),
    fPhi(phi),
    fEdep(edep)
{}


//____________________________________________________________________
void
AliFMDMultStrip::Print(Option_t* option) const
{
  // Print information 
  // 
  // Options:
  //    D:           Detector (default)
  //    E:           Eta range (default)
  //    P:           Phi range (default)
  //
  TString opt(option);
  cout << "FMD Multiplicity in a strip: " << fParticles << endl;
  if (opt.Contains("D", TString::kIgnoreCase))
    cout << "  Detector:      FMD" << fDetector << fRing 
	 << "[" << fSector << "," << fStrip << "]" << endl;
  if (opt.Contains("E", TString::kIgnoreCase))
    cout << "  Eta range:     " << fEta << endl;
  if (opt.Contains("P", TString::kIgnoreCase))
    cout << "  Phi range:     " << fPhi << endl;
  AliFMDMult::Print(option);
}

    
//____________________________________________________________________
//
// EOF
//
