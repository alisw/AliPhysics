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
/** @file    AliFMDDigit.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:37:41 2006
    @brief   Digits for the FMD 
*/
//////////////////////////////////////////////////////////////////////
//
//  Class that holds an FMD index.  That is, it holds the detector
//  coordinates for a given strip:
//
//     Variable | Type     | Range   | Description
//     ---------+----------+---------+------------------
//     detector | UShort_t | 1-3     | Detector number 
//     ring     | Char_t   | 'I'/'O' | Ring identifier 
//     sector   | UShort_t | 0-39    | Sector number
//     strip    | UShort_t | 0-511   | Strip number
//
//////////////////////////////////////////////////////////////////////

#include "AliFMDIndex.h"	// ALIFMDINDEX_H
#include "Riostream.h"		// ROOT_Riostream
#include <TString.h>            // ROOT_TString
#include <AliFMDMap.h>

//====================================================================
ClassImp(AliFMDIndex)
#if 0
  ; // This is here to keep Emacs from indenting the next line
#endif

//____________________________________________________________________
AliFMDIndex::AliFMDIndex()
  : fDetector(0), 
    fRing('\0'), 
    fSector(0), 
    fStrip(0), 
    fName(""),
    fHash(-1) 
{
  // CTOR
}

//____________________________________________________________________
AliFMDIndex::AliFMDIndex(const AliFMDIndex& o)
  : fDetector(o.fDetector), 
    fRing(o.fRing), 
    fSector(o.fSector), 
    fStrip(o.fStrip), 
    fName(""),
    fHash(o.fHash)
{
  // Copy constructor 
}

//____________________________________________________________________
AliFMDIndex::AliFMDIndex(UShort_t detector, 
			 Char_t   ring, 
			 UShort_t sector, 
			 UShort_t strip)
  : fDetector(detector), 
    fRing(ring), 
    fSector(sector), 
    fStrip(strip), 
    fName(""),
    fHash(-1)
{
  //
  // Creates a base data digit object
  //
  // Parameters 
  //
  //    detector  Detector # (1, 2, or 3)                      
  //    ring	  Ring ID ('I' or 'O')
  //    sector	  Sector # (For inner/outer rings: 0-19/0-39)
  //    strip	  Strip # (For inner/outer rings: 0-511/0-255)
}

//____________________________________________________________________
AliFMDIndex& 
AliFMDIndex::operator=(const AliFMDIndex& o)
{
  // Assignment operator 
  if (&o == this) return *this; 
  fDetector = o.fDetector;
  fRing     = o.fRing;
  fSector   = o.fSector;
  fStrip    = o.fStrip;
  fHash     = o.fHash;
  return *this;
}

//____________________________________________________________________
Int_t
AliFMDIndex::Hash() const 
{
  // calculate hash value 
  if (fHash < 0) {
    size_t ringi = (fRing == 'I' ||  fRing == 'i' ? 0 : 1);
    fHash = (fStrip + 
	     AliFMDMap::kMaxStrips * 
	     (fSector + AliFMDMap::kMaxSectors * 
	      (ringi + AliFMDMap::kMaxRings * (fDetector-1))));
  }
  return fHash;
}


//____________________________________________________________________
void
AliFMDIndex::Print(Option_t* /* option*/) const 
{
  // Print digit to standard out 
  cout << Name() << flush;
}

//____________________________________________________________________
const char*
AliFMDIndex::Name() const 
{ 
  // GEt the name of the index 
  if (fName.IsNull()) 
    fName = Form("FMD%d%c[%2d,%3d]", fDetector, fRing, fSector, fStrip);
  return fName.Data();
}

//====================================================================
ClassImp(AliFMDObjIndex)
#if 0
  ; // This is here to keep Emacs from indenting the next line
#endif

//____________________________________________________________________
Int_t 
AliFMDObjIndex::Compare(const TObject* o) const
{
  // Compare to another index 
  const AliFMDObjIndex* a = dynamic_cast<const AliFMDObjIndex*>(o);
  if (!a) {
    Fatal("Compare", 
	  "trying to compare to something not a AliFMDObjIndex object, "
	  "but a %s object", o->ClassName());
    return 0;
  }
  if (this->operator<(*a)) return -1;
  if (this->operator==(*a)) return 0;
  return 1;
}

//____________________________________________________________________
//
// EOF
//
