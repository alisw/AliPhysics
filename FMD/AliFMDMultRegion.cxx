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
//  Forward Multiplicity Detector have to be reconstructed number of
//  particles in fixed pseudorapidity interval from fNumOfMinRing
//  to fNumOfMaxRing and phi interval from fNumOfMinSector to
//  fNumOfMaxSector
//
#include "AliFMDMultRegion.h"	// ALIFMDPARTICLES_H
#include <TString.h>            // ROOT_TString
#include <Riostream.h>		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliFMDMultRegion);


//____________________________________________________________________
AliFMDMultRegion::AliFMDMultRegion()
  : fDetector(0),
    fRing('\0'),
    fMinSector(0),
    fMaxSector(0),
    fMinStrip(0),
    fMaxStrip(0),
    fMinEta(0),
    fMaxEta(0),
    fMinPhi(0),
    fMaxPhi(0)
{}

//____________________________________________________________________
AliFMDMultRegion::AliFMDMultRegion(UShort_t detector,  Char_t ring, 
				   UShort_t minSector, UShort_t maxSector, 
				   UShort_t minStrip,  UShort_t maxStrip, 
				   Float_t  minEta,    Float_t  maxEta, 
				   Float_t  minPhi,    Float_t  maxPhi,
				   Float_t  particles, UShort_t method)
  : AliFMDMult(particles, method),
    fDetector(detector),
    fRing(ring),
    fMinSector(minSector),
    fMaxSector(maxSector),
    fMinStrip(minStrip),
    fMaxStrip(maxStrip),
    fMinEta(minEta),
    fMaxEta(maxEta),
    fMinPhi(minPhi),
    fMaxPhi(maxPhi)
{}


//____________________________________________________________________
void
AliFMDMultRegion::Print(Option_t* option) const
{
  // Print information 
  // 
  // Options:
  //    D:           Detector (default)
  //    S:           Sector range 
  //    T:           Strip range 
  //    E:           Eta range (default)
  //    P:           Phi range (default)
  //
  TString opt(option);
  cout << "FMD Reconstructed particles: " << fParticles << endl;
  if (opt.Contains("D", TString::kIgnoreCase))
    cout << "  Detector:      FMD" << fDetector << fRing << endl;
  if (opt.Contains("S", TString::kIgnoreCase))
    cout << "  Sector range:  [" << fMinSector << "," << fMaxSector << endl;
  if (opt.Contains("T", TString::kIgnoreCase))
    cout << "  Strip range:   [" << fMinStrip << "," << fMaxStrip << endl;
  if (opt.Contains("E", TString::kIgnoreCase))
    cout << "  Eta range:     [" << fMinEta << "," << fMaxEta << endl;
  if (opt.Contains("P", TString::kIgnoreCase))
    cout << "  Phi range:     [" << fMinPhi << "," << fMaxPhi << endl;
  AliFMDMult::Print(option);
}

    
//____________________________________________________________________
//
// EOF
//
