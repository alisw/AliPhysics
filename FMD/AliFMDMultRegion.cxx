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
// FMD reconstructed multiplicity in a region of a ring.  The region
// is identified by (strip_min,sector_min)x(strip_max,sector_max) or
// by (eta_max,phi_min),(eta_min,phi_max).   It's also possible to
// get the peudorapidity of the center-of-mass of the region, as well
// as the mean azimuthal angle.
//
// [Note, that minStrip corresponds to maxEta, and maxStrip
// corresponds to minEta]
// 
// These objects are usually created by the Poisson reconstruction
// method. 
//
#include "AliFMDMultRegion.h"	// ALIFMDPARTICLES_H
#include <TString.h>            // ROOT_TString
#include <Riostream.h>		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliFMDMultRegion)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif


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
Float_t
AliFMDMultRegion::Eta() const 
{
  // Return the center-of-mass eta of the region.   This is calculated
  // as the weighted mean of min and max eta, where the weights are
  // the length of the arcs of the upper and lower edge of the region: 
  // 
  //            (maxPhi - minPhi)
  //  f       = -----------------
  //                   360 	
  //
  //  w_max   = ds*minStrip*2*pi*f
  //                                   
  //  w_min   = ds*maxStrip*2*pi*f
  // 
  //  1/w^2   = 1/w_max^2 + 1/w_min^2 
  //
  //                      1                        1  	      
  //          = ---------------------- + ----------------------
  //            (ds*minStrip*2*pi*f)^2	 (ds*maxStrip*2*pi*f)^2
  // 
  //            (ds*maxStrip*2*pi*f)^2 + (ds*minStrip*2*pi*f)^2
  //          = -----------------------------------------------
  //            (ds*maxStrip*2*pi*f)^2 * (ds*minStrip*2*pi*f)^2
  //
  //          
  //             4 * pi^2 * ds^2 * f^2 (maxStrip^2  + minStrip^2)
  //          = -------------------------------------------------
  //             16 * pi^4 * ds^4 * f^4 minStrip^2 * maxStrip^2 
  //
  //                     (maxStrip^2  + minStrip^2)
  //          = -------------------------------------------------
  //             4 * pi^2 * ds^2 * f^2 minStrip^2 * maxStrip^2 
  //
  //  <eta> = (maxEta/w_max^2 + minEta /w_min^2) / (1/w^2)
  //      
  //           w_min^2 * maxEta + w_max^2 * minEta
  //        = ------------------------------------ / (1/w^2)
  //                  w_max^2 * w_min^2
  // 
  //           4*pi^2*ds^2*(maxStrip*maxEta + minStrip*minEta)
  //        =  ----------------------------------------------- / (1/w^2)
  //               16*pi^4*ds^4*f^4*minStrip^2*maxStrip^2 
  //
  //           maxStrip * maxEta + minStrip * minEta   
  //        =  ------------------------------------- 
  //           4*pi^2*ds^2*f^2*minStrip^2*maxStrip^2 
  //
  //               4*pi^2*ds^2*f^2*minStrip^2*maxStrip^2
  //             * ------------------------------------- 		  
  //		         maxStrip^2+minStrip^2		  
  // 
  //           maxStrip * maxEta + minStrip * minEta   
  //        =  ------------------------------------- 
  //                  maxStrip^2 + minStrip^2
  //                                                     _
  //                                                    |_|
  //
  Float_t eta = (fMaxStrip * fMaxEta + fMinStrip * fMinEta) 
    / (fMaxStrip * fMaxStrip + fMinStrip * fMinStrip);
  return eta;
}

  
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
  cout << "FMD Multiplicity in a region: " << fParticles << endl;
  if (opt.Contains("D", TString::kIgnoreCase))
    cout << "  Detector:      FMD" << fDetector << fRing 
	 << "[" << fMinSector << "-" << fMaxSector 
	 << "," <<  fMinStrip << "-" << fMaxStrip << "]" << endl;
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
