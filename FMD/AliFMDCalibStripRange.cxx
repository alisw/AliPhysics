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
/** @file    AliFMDCalibStripRange.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:31:09 2006
    @brief   Per digitizer card pulser calibration 
*/
//____________________________________________________________________
//                                                                          
// This class stores which strips are read-out.  
// In principle this can be set for each half-ring.   
// However, in real life, all the detectors will probably read out all
// strips, and dead areas can be handled off-line. 
// This information comes from DCS or the like.
//
#include "AliFMDCalibStripRange.h"	// ALIFMDCALIBGAIN_H
#include "AliFMDParameters.h"           // ALIFMDPARAMETERS_H

//____________________________________________________________________
ClassImp(AliFMDCalibStripRange)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDCalibStripRange::AliFMDCalibStripRange()
  : fRanges(AliFMDMap::kMaxDetectors, AliFMDMap::kMaxRings, 2, 1)
  // fRanges(3)
{
  // CTOR 
  fRanges.Reset(1);
}

//____________________________________________________________________
AliFMDCalibStripRange::AliFMDCalibStripRange(const AliFMDCalibStripRange& o)
  : TObject(o), fRanges(o.fRanges)
{
  // Copy CTOR 
}

//____________________________________________________________________
AliFMDCalibStripRange&
AliFMDCalibStripRange::operator=(const AliFMDCalibStripRange& o)
{
  // Assignement operator
  fRanges     = o.fRanges;
  return (*this);
}

//____________________________________________________________________
void
AliFMDCalibStripRange::Set(UShort_t det, Char_t ring, 
			   UShort_t sector, UShort_t, UShort_t min, 
			   UShort_t max)
{
  // Set the min and max for a half-ring 
  UInt_t nSec  = (ring == 'I' ? 20 : 40);
  UInt_t board = sector / nSec;
  fRanges(det, ring, board, 0) = ((max & 0x7f) << 8) + (min & 0x7f);
}

//____________________________________________________________________
UShort_t
AliFMDCalibStripRange::Min(UShort_t det, Char_t ring, 
			   UShort_t sec, UShort_t) const
{
  // Get the min for a half-ring 
  UInt_t nSec  = (ring == 'I' ? 20 : 40);
  UInt_t board = sec / nSec;
  return (fRanges(det, ring, board, 0) & 0x7f);
}

//____________________________________________________________________
UShort_t
AliFMDCalibStripRange::Max(UShort_t det, Char_t ring, 
			   UShort_t sec, UShort_t) const
{
  // Get the max for a half-ring 
  UInt_t nSec  = (ring == 'I' ? 20 : 40);
  UInt_t board = sec / nSec;
  return ((fRanges(det, ring, board, 0) >> 8) & 0x7f);
}

//____________________________________________________________________
//
// EOF
//
