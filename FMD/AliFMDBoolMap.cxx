/**************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN.          *
 * All rights reserved.                                       *
 *                                                            *
 * Author: The ALICE Off-line Project.                        *
 * Contributors are mentioned in the code where appropriate.  *
 *                                                            *
 * Permission to use, copy, modify and distribute this        *
 * software and its documentation strictly for non-commercial *
 * purposes is hereby granted without fee, provided that the  *
 * above copyright notice appears in all copies and that both *
 * the copyright notice and this permission notice appear in  *
 * the supporting documentation. The authors make no claims   *
 * about the suitability of this software for any purpose. It *
 * is provided "as is" without express or implied warranty.   *
 **************************************************************/
/* $Id$ */
/** @file    AliFMDBoolMap.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:28:42 2006
    @brief   Per strip Boolean map
*/
//__________________________________________________________
// 
// Map of Bool_t for each FMD strip
// Used in calibration and the like classes.
// Used amoung other things for dead-channel map
// Can also be used for other stuff too
// Created Mon Nov  8 12:51:51 2004 by Christian Holm Christensen
// 
#include "AliFMDBoolMap.h"	//ALIFMDBOOLMAP_H
//__________________________________________________________
ClassImp(AliFMDBoolMap)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif
//__________________________________________________________
AliFMDBoolMap::AliFMDBoolMap(const AliFMDBoolMap& other)
  : AliFMDMap(other.fMaxDetectors,
              other.fMaxRings,
              other.fMaxSectors,
              other.fMaxStrips),
    fTotal(0),
    fData(0)
{
  // Copy constructor
  fTotal = fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
  fData  = new Bool_t[fTotal];
  for (Int_t i = 0; i < fTotal; i++) fData[i] = other.fData[i];
}

//__________________________________________________________
AliFMDBoolMap::AliFMDBoolMap(UShort_t maxDet,
                         UShort_t maxRing,
                         UShort_t maxSec,
                         UShort_t maxStr)
  : AliFMDMap(maxDet, maxRing, maxSec, maxStr),
    fTotal(0),
    fData(0)
{
  // Constructor.
  // Parameters:
  //	maxDet	Maximum number of detectors
  //	maxRing	Maximum number of rings per detector
  //	maxSec	Maximum number of sectors per ring
  //	maxStr	Maximum number of strips per sector
  fTotal = fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
  fData  = new Bool_t[fTotal];
  Reset();
}

//__________________________________________________________
AliFMDBoolMap&
AliFMDBoolMap::operator=(const AliFMDBoolMap& other)
{
  // Assignment operator 
  fMaxDetectors = other.fMaxDetectors;
  fMaxRings     = other.fMaxRings;
  fMaxSectors   = other.fMaxSectors;
  fMaxStrips    = other.fMaxStrips;
  if (fData) delete [] fData;
  fTotal = fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips;
  fData  = new Bool_t[fTotal];
  for (Int_t i = 0; i < fTotal; i++) fData[i] = other.fData[i];
  return *this;
}

//__________________________________________________________
void
AliFMDBoolMap::Reset(const Bool_t& val)
{
  // Reset map to val
  for (Int_t i = 0; i < fTotal; i++) fData[i] = val;
}

//__________________________________________________________
Bool_t&
AliFMDBoolMap::operator()(UShort_t det, 
			  Char_t   ring, 
			  UShort_t sec, 
			  UShort_t str)
{
  // Get data
  // Parameters:
  //	det	Detector #
  //	ring	Ring ID
  //	sec	Sector #
  //	str	Strip #
  // Returns appropriate data
  return fData[CalcIndex(det, ring, sec, str)];
}

//__________________________________________________________
const Bool_t&
AliFMDBoolMap::operator()(UShort_t det, 
			  Char_t   ring, 
			  UShort_t sec, 
			  UShort_t str) const
{
  // Get data
  // Parameters:
  //	det	Detector #
  //	ring	Ring ID
  //	sec	Sector #
  //	str	Strip #
  // Returns appropriate data
  return fData[CalcIndex(det, ring, sec, str)];
}

//__________________________________________________________
// 
// EOF
// 

