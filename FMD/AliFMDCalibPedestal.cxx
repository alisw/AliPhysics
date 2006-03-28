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
/** @file    AliFMDCalibPedestal.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:30:36 2006
    @brief   Per strip pedestal calibration 
*/
//____________________________________________________________________
//                                                                          
//
//
#include "AliFMDCalibPedestal.h"	// ALIFMDCALIBPEDESTAL_H
//____________________________________________________________________
ClassImp(AliFMDCalibPedestal)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDCalibPedestal::AliFMDCalibPedestal()
{
  fValue.Reset(-1.);
  fWidth.Reset(-1.);
}

//____________________________________________________________________
AliFMDCalibPedestal::AliFMDCalibPedestal(const AliFMDCalibPedestal& o)
  : TObject(o), fValue(o.fValue), fWidth(o.fWidth)
{}

//____________________________________________________________________
AliFMDCalibPedestal&
AliFMDCalibPedestal::operator=(const AliFMDCalibPedestal& o)
{
  fValue = o.fValue;
  fWidth = o.fWidth;
  return (*this);
}

//____________________________________________________________________
void
AliFMDCalibPedestal::Set(UShort_t det, Char_t ring, UShort_t sec, 
			 UShort_t str, Float_t ped, Float_t pedW)
{
  if (fValue.CheckIndex(det, ring, sec, str) < 0) return;
  fValue(det, ring, sec, str) = ped;
  fWidth(det, ring, sec, str) = pedW;
}

//____________________________________________________________________
Float_t
AliFMDCalibPedestal::Value(UShort_t det, Char_t ring, UShort_t sec, 
			   UShort_t str)
{
  return fValue(det, ring, sec, str);
}

//____________________________________________________________________
Float_t
AliFMDCalibPedestal::Width(UShort_t det, Char_t ring, UShort_t sec, 
			   UShort_t str)
{
  return fWidth(det, ring, sec, str);
}

//____________________________________________________________________
//
// EOF
//
