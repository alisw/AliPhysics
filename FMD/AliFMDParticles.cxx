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
#include "AliFMDParticles.h"	// ALIFMDPARTICLES_H
#include "Riostream.h"		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliFMDParticles)

//____________________________________________________________________
AliFMDParticles::AliFMDParticles()
  : fDetector(0),
    fRing('\0'),
    fMinSector(0),
    fMaxSector(0),
    fMinStrip(0),
    fMaxStrip(0),
    fMinEta(0),
    fMaxEta(0),
    fMinPhi(0),
    fMaxPhi(0),
    fParticles(0),
    fMethod(kNaive)
{}

//____________________________________________________________________
AliFMDParticles::AliFMDParticles(UShort_t detector,  Char_t ring, 
				 UShort_t minSector, UShort_t maxSector, 
				 UShort_t minStrip,  UShort_t maxStrip, 
				 Float_t  minEta,    Float_t  maxEta, 
				 Float_t  minPhi,    Float_t  maxPhi,
				 Float_t  particles, UShort_t method)
  : fDetector(detector),
    fRing(ring),
    fMinSector(minSector),
    fMaxSector(maxSector),
    fMinStrip(minStrip),
    fMaxStrip(maxStrip),
    fMinEta(minEta),
    fMaxEta(maxEta),
    fMinPhi(minPhi),
    fMaxPhi(maxPhi),
    fParticles(particles),
    fMethod(method)
{
  switch (fMethod) {
  case kPoission: 
  case kIterative: 
  case kNaive:    
    break;    
  default:        
    Warning("AliFMDParticles", "unknown method: %d", method);
    break;
  }
}


//____________________________________________________________________
void
AliFMDParticles::Print(Option_t* /* opt*/) const
{
  cout << "FMD Reconstructed particles: " << fParticles << "\n" 
       << "  Detector:      FMD" << fDetector << fRing << "\n"
       << "  Sector range:  [" << fMinSector << "," << fMaxSector << "\n"
       << "  Strip range:   [" << fMinStrip << "," << fMaxStrip << "\n"
       << "  Eta range:     [" << fMinEta << "," << fMaxEta << "\n"
       << "  Phi range:     [" << fMinPhi << "," << fMaxPhi << "\n"
       << "  Method:        " << flush;
  switch (fMethod) {
  case kPoission:  cout << "Poission"  << endl; break;
  case kIterative: cout << "Iterative" << endl; break;
  case kNaive:     cout << "Naive"     << endl; break; 
  default:         cout << "Unknown"   << endl; break;
  }
}

    
//____________________________________________________________________
//
// EOF
//
