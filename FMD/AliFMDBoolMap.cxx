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
//__________________________________________________________
// 
// Map of per strip Bool_t information
// 
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
    fData(0)
{
  // Copy constructor
  fData = new Bool_t[fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips];
  for (size_t i = 0; i < fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips; i++)
    fData[i] = other.fData[i];
}

//__________________________________________________________
AliFMDBoolMap::AliFMDBoolMap(size_t maxDet,
                         size_t maxRing,
                         size_t maxSec,
                         size_t maxStr)
  : AliFMDMap(maxDet, maxRing, maxSec, maxStr),
    fData(0)
{
  // Constructor.
  // Parameters:
  //	maxDet	Maximum number of detectors
  //	maxRing	Maximum number of rings per detector
  //	maxSec	Maximum number of sectors per ring
  //	maxStr	Maximum number of strips per sector
  fData = new Bool_t[fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips];
  Clear();
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
  fData = new Bool_t[fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips];
  for (size_t i = 0; i < fMaxDetectors * fMaxRings 
	 * fMaxSectors * fMaxStrips; i++)
    fData[i] = other.fData[i];
  return *this;
}

//__________________________________________________________
void
AliFMDBoolMap::Clear(const Bool_t& val)
{
  // Reset map to val
  for (size_t i = 0; i < fMaxDetectors * fMaxRings 
	 * fMaxSectors * fMaxStrips; i++)
    fData[i] = val;
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

