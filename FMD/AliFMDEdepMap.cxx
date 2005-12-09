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

//____________________________________________________________________
//                                                                          
//
//
#include "AliFMDEdepMap.h"		// ALIFMDEDEPMAP_H

//____________________________________________________________________
ClassImp(AliFMDEdepMap)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDEdepMap::AliFMDEdepMap(const AliFMDEdepMap& other)
  : AliFMDMap(other.fMaxDetectors, other.fMaxRings, other.fMaxSectors, 
	      other.fMaxStrips), 
    fData(0)
{
  // Copy constructor 
  fData = new AliFMDEdepHitPair[fMaxDetectors * fMaxRings * 
				fMaxSectors * fMaxStrips];
  for (size_t i = 0; i < fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
       i++) fData[i] = other.fData[i];
}

  

//____________________________________________________________________
AliFMDEdepMap::AliFMDEdepMap(size_t maxDet, 
			     size_t maxRing, 
			     size_t maxSec, 
			     size_t maxStr)
  : AliFMDMap(maxDet, maxRing, maxSec, maxStr), 
    fData(0)
{
  // Construct a map
  //
  // Parameters:
  //     maxDet       Maximum # of detectors
  //     maxRinf      Maximum # of rings
  //     maxSec       Maximum # of sectors
  //     maxStr       Maximum # of strips
  fData = new AliFMDEdepHitPair[fMaxDetectors * fMaxRings * 
				fMaxSectors * fMaxStrips];
}

//____________________________________________________________________
AliFMDEdepMap&
AliFMDEdepMap::operator=(const AliFMDEdepMap& other) 
{
  // Assignment operator
  fMaxDetectors = other.fMaxDetectors;
  fMaxRings     = other.fMaxRings;
  fMaxSectors   = other.fMaxSectors;
  fMaxStrips    = other.fMaxStrips;
  if (fData) delete [] fData;
  fData = new AliFMDEdepHitPair[fMaxDetectors * fMaxRings * 
				fMaxSectors * fMaxStrips];
  for (size_t i = 0; i < fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
       i++) fData[i] = other.fData[i];
  return *this;
}

//____________________________________________________________________
void
AliFMDEdepMap::Reset() 
{
  // Reset to zero
  for (size_t i = 0; i < fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
       i++) { fData[i].fEdep = 0; fData[i].fN = 0; };
}

//____________________________________________________________________
void
AliFMDEdepMap::Reset(const AliFMDEdepHitPair& val) 
{
  // Reset to val
  for (size_t i = 0; i < fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
       i++) { fData[i].fEdep = val.fEdep; fData[i].fN = val.fN; };
}

//____________________________________________________________________
AliFMDEdepHitPair& 
AliFMDEdepMap::operator()(UShort_t det, Char_t ring, UShort_t sec, UShort_t str) 
{
  // Get data 
  // 
  // Parameters: 
  //     det       Detector #
  //     ring      Ring ID
  //     sec       Sector # 
  //     str       Strip # 
  //
  // Returns appropriate data
  //
  return fData[CalcIndex(det, ring, sec, str)];
}

//____________________________________________________________________
const AliFMDEdepHitPair& 
AliFMDEdepMap::operator()(UShort_t det, Char_t ring, UShort_t sec, UShort_t str) const
{
  // Get data 
  // 
  // Parameters: 
  //     det       Detector #
  //     ring      Ring ID
  //     sec       Sector # 
  //     str       Strip # 
  //
  // Returns appropriate data
  //
  return fData[CalcIndex(det, ring, sec, str)];
}


//___________________________________________________________________
//
// EOF
//
