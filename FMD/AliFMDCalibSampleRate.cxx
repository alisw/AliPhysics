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
#include "AliFMDCalibSampleRate.h"	// ALIFMDCALIBGAIN_H
#include "AliFMDParameters.h"           // ALIFMDPARAMETERS_H

//____________________________________________________________________
ClassImp(AliFMDCalibSampleRate)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDCalibSampleRate::AliFMDCalibSampleRate()
  : fRates(3)
{
  fRates.Reset(0);
}

//____________________________________________________________________
AliFMDCalibSampleRate::AliFMDCalibSampleRate(const AliFMDCalibSampleRate& o)
  : TObject(o), fRates(o.fRates)
{}

//____________________________________________________________________
AliFMDCalibSampleRate&
AliFMDCalibSampleRate::operator=(const AliFMDCalibSampleRate& o)
{
  fRates     = o.fRates;
  return (*this);
}

//____________________________________________________________________
void
AliFMDCalibSampleRate::Set(UShort_t ddl, UShort_t rate)
{
  if (ddl - AliFMDParameters::kBaseDDL < 0) return;
  fRates[ddl - AliFMDParameters::kBaseDDL] = rate;
}

//____________________________________________________________________
UShort_t
AliFMDCalibSampleRate::Rate(UShort_t ddl) const
{
  if (ddl - AliFMDParameters::kBaseDDL < 0) return 0;
  return fRates[ddl - AliFMDParameters::kBaseDDL];
}

//____________________________________________________________________
//
// EOF
//
