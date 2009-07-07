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
/** @file    AliFMDCalibGain.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:30:02 2006
    @brief   Per strip gain calibration 
*/
//____________________________________________________________________
//                                                                          
// Gain value and width for each strip in the FMD. 
// Foo 
// Bar 
// Baz
// Gnus
//
#include "AliFMDCalibGain.h"	// ALIFMDCALIBGAIN_H
//____________________________________________________________________
ClassImp(AliFMDCalibGain)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDCalibGain::AliFMDCalibGain()
  : fValue(0), // nDet == 0 mean 51200 slots
    fThreshold(-1.)
{
  // CTOR
  fValue.Reset(-1.);
  fThreshold = -1.;
}

//____________________________________________________________________
AliFMDCalibGain::AliFMDCalibGain(const AliFMDCalibGain& o)
  : TObject(o), 
    fValue(o.fValue), 
    fThreshold(o.fThreshold)
{
  // Copy CTOR 
}

//____________________________________________________________________
AliFMDCalibGain&
AliFMDCalibGain::operator=(const AliFMDCalibGain& o)
{
  // Assignment operator 
  fValue     = o.fValue;
  fThreshold = o.fThreshold;
  return (*this);
}

//____________________________________________________________________
void
AliFMDCalibGain::Set(UShort_t det, Char_t ring, UShort_t sec, 
		     UShort_t str, Float_t val)
{
  // Set the value for a strip 
  if (fValue.CheckIndex(det, ring, sec, str) < 0) return;
  fValue(det, ring, sec, str) = val;
}

//____________________________________________________________________
Float_t
AliFMDCalibGain::Value(UShort_t det, Char_t ring, UShort_t sec, 
		       UShort_t str)
{
  // Get the value for a strip 
  return fValue(det, ring, sec, str);
}

//____________________________________________________________________
//
// EOF
//
