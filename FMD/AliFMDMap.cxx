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
#include "AliFMDMap.h"		// ALIFMDMAP_H

//____________________________________________________________________
ClassImp(AliFMDMap)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDMap::AliFMDMap(size_t maxDet, 
		     size_t maxRing, 
		     size_t maxSec, 
		     size_t maxStr)
  : fMaxDetectors(maxDet), 
    fMaxRings(maxRing), 
    fMaxSectors(maxSec), 
    fMaxStrips(maxStr)
{
  // Construct a map
  //
  // Parameters:
  //     maxDet       Maximum # of detectors
  //     maxRinf      Maximum # of rings
  //     maxSec       Maximum # of sectors
  //     maxStr       Maximum # of strips
}


//____________________________________________________________________
size_t 
AliFMDMap::CalcIndex(size_t det, Char_t ring, size_t sec, size_t str) const
{
  // Calculate index into storage from arguments. 
  // 
  // Parameters: 
  //     det       Detector #
  //     ring      Ring ID
  //     sec       Sector # 
  //     str       Strip # 
  //
  // Returns appropriate index into storage 
  //
  size_t ringi = (ring == 'I' ||  ring == 'i' ? 0 : 1);
  size_t idx = 
    (det + fMaxDetectors * (ringi + fMaxRings * (sec + fMaxSectors * str)));
  if (idx >= fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips) {
    Fatal("CalcIndex", "Index (%d,'%c',%d,%d) out of bounds, "
	  "in particular the %s index", 
	  det, ring, sec, str, 
	  (det >= fMaxDetectors ? "Detector" : 
	   (ringi >= fMaxRings ? "Ring" : 
	    (sec >= fMaxSectors ? "Sector" : "Strip"))));
    return 0;
  }
  return idx;
}


//___________________________________________________________________
//
// EOF
//
