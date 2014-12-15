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
/** @file    AliFMDEdepMap.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:39:50 2006
    @brief   Per strip map of energy deposited and number of hits 
    @ingroup FMD_sim
*/
//____________________________________________________________________
//                                                                          
// Contains a pair of energy deposited fEdep and number of hits  
// fN, fEdep is the summed energy deposition, and fN is the
// number of hits.  The map contains one such object or each strip.
// It is used to cache the data in the digitization classes
// AliFMDBaseDigitizer and so on. 
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
    fTotal(fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips),
    fData(0)
{
  // Copy constructor 
  if (fTotal == 0) fTotal = 51200;
  fData  = new AliFMDEdepHitPair[fTotal];
  for (Int_t i = 0; i < fTotal; i++) fData[i] = other.fData[i];
}

  
//____________________________________________________________________
AliFMDEdepMap::AliFMDEdepMap()
  : AliFMDMap(), 
    fTotal(0),
    fData(0)
{
  // Construct a map
  //
  // Parameters:
  //   None
}

//____________________________________________________________________
AliFMDEdepMap::AliFMDEdepMap(UShort_t maxDet, 
			     UShort_t maxRing, 
			     UShort_t maxSec, 
			     UShort_t maxStr)
  : AliFMDMap(maxDet, maxRing, maxSec, maxStr), 
    fTotal(fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips),
    fData(0)
{
  // Construct a map
  //
  // Parameters:
  //     maxDet       Maximum # of detectors
  //     maxRinf      Maximum # of rings
  //     maxSec       Maximum # of sectors
  //     maxStr       Maximum # of strips
  if (fTotal == 0) fTotal = 51200;
  fData  = new AliFMDEdepHitPair[fTotal];
}

//____________________________________________________________________
AliFMDEdepMap&
AliFMDEdepMap::operator=(const AliFMDEdepMap& other) 
{
  // Assignment operator
  if (&other == this) return *this; 
  fMaxDetectors = other.fMaxDetectors;
  fMaxRings     = other.fMaxRings;
  fMaxSectors   = other.fMaxSectors;
  fMaxStrips    = other.fMaxStrips;
  if (fData) delete [] fData;
  fTotal = fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
  if (fTotal == 0) fTotal = 51200;
  fData  = new AliFMDEdepHitPair[fTotal];
  for (Int_t i = 0; i < fTotal; i++) fData[i] = other.fData[i];
  return *this;
}

//____________________________________________________________________
void
AliFMDEdepMap::Reset() 
{
  // Reset to zero
  for (Int_t i = 0; i < fTotal; i++) { 
    fData[i].fEdep  = 0; 
    fData[i].fN     = 0; 
    fData[i].fNPrim = 0;
    fData[i].fLabels.Reset();
  };
}

//____________________________________________________________________
void
AliFMDEdepMap::Reset(const AliFMDEdepHitPair& val) 
{
  // Reset to val
  for (Int_t i = 0; i < fTotal; i++) { 
    fData[i].fEdep   = val.fEdep; 
    fData[i].fN      = val.fN; 
    fData[i].fNPrim  = val.fNPrim;
    fData[i].fLabels = val.fLabels;
  };
}

//____________________________________________________________________
AliFMDEdepHitPair& 
AliFMDEdepMap::operator()(UShort_t det, Char_t ring, 
			  UShort_t sec, UShort_t str) 
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
AliFMDEdepMap::operator()(UShort_t det, Char_t ring, 
			  UShort_t sec, UShort_t str) const
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
