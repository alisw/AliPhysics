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
/** @file    AliFMDUShortMap.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:48:18 2006
    @brief   Per strip of unisgned shorts (16 bit) data 
*/
//____________________________________________________________________
//                                                                          
// A map of per strip UShort_t information (for example ADC values,
// number of hits and so on). 
// Used for various calib information, and the like, 
// as well as in reconstruction. 
// Can be used elsewhere too.
//
#include "AliFMDUShortMap.h"		// ALIFMDUSHORTMAP_H

//____________________________________________________________________
ClassImp(AliFMDUShortMap)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDUShortMap::AliFMDUShortMap(const AliFMDUShortMap& other)
  : AliFMDMap(other.fMaxDetectors, other.fMaxRings, other.fMaxSectors, 
	      other.fMaxStrips), 
    fTotal(0),
    fData(0)
{
  // CTOR
  fTotal = fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
  fData  = new UShort_t[fTotal];
  for (Int_t i = 0; i < fTotal; i++) fData[i] = other.fData[i];
}

  

//____________________________________________________________________
AliFMDUShortMap::AliFMDUShortMap(UShort_t maxDet, 
				 UShort_t maxRing, 
				 UShort_t maxSec, 
				 UShort_t maxStr)
  : AliFMDMap(maxDet, maxRing, maxSec, maxStr), 
    fTotal(0),
    fData(0)
{
  // Construct a map
  //
  // Parameters:
  //     maxDet       Maximum # of detectors
  //     maxRinf      Maximum # of rings
  //     maxSec       Maximum # of sectors
  //     maxStr       Maximum # of strips
  fTotal = fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
  fData  = new UShort_t[fTotal];
}

//____________________________________________________________________
AliFMDUShortMap&
AliFMDUShortMap::operator=(const AliFMDUShortMap& other) 
{
  // Assignment operator
  fMaxDetectors = other.fMaxDetectors;
  fMaxRings     = other.fMaxRings;
  fMaxSectors   = other.fMaxSectors;
  fMaxStrips    = other.fMaxStrips;
  if (fData) delete [] fData;
  fTotal = fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
  fData  = new UShort_t[fTotal];
  for (Int_t i = 0; i < fTotal; i++) fData[i] = other.fData[i];
  return *this;
}

//____________________________________________________________________
void
AliFMDUShortMap::Reset(const UShort_t& val) 
{
  // Reset to val
  for (Int_t i = 0; i < fTotal; i++) fData[i] = val;
}

//____________________________________________________________________
UShort_t& 
AliFMDUShortMap::operator()(UShort_t det, Char_t ring, UShort_t sec, UShort_t str) 
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
const UShort_t& 
AliFMDUShortMap::operator()(UShort_t det, Char_t ring, UShort_t sec, UShort_t str) const
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
