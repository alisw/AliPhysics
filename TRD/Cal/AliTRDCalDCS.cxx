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

/* $Id: AliTRDCalDCS.cxx 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS parameters                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDCS.h"

ClassImp(AliTRDCalDCS)

//_____________________________________________________________________________
AliTRDCalDCS::AliTRDCalDCS()
  :TNamed()
  ,fNumberOfTimeBins(0)
  ,fTailCancelationTau1(0)
  ,fTailCancelationTau2(0)
  ,fTailCancelationAmp(0)
  ,fPedestal(0)
  ,fConfigID(0)
  ,fGainTableID(0)
  ,fFEEArr(new TObjArray(540))
  ,fPTRArr(new TObjArray(6))
  ,fGTUArr(new TObjArray(19))
{
  //
  // AliTRDCalDCS default constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCS::AliTRDCalDCS(const Text_t *name, const Text_t *title)
  :TNamed(name,title)
  ,fNumberOfTimeBins(0)
  ,fTailCancelationTau1(0)
  ,fTailCancelationTau2(0)
  ,fTailCancelationAmp(0)
  ,fPedestal(0)
  ,fConfigID(0)
  ,fGainTableID(0)
  ,fFEEArr(new TObjArray(540))
  ,fPTRArr(new TObjArray(6))
  ,fGTUArr(new TObjArray(19))
{
  //
  // AliTRDCalDCS constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCS::AliTRDCalDCS(const AliTRDCalDCS &cd)
  :TNamed(cd)
  ,fNumberOfTimeBins(0)
  ,fTailCancelationTau1(0)
  ,fTailCancelationTau2(0)
  ,fTailCancelationAmp(0)
  ,fPedestal(0)
  ,fConfigID(0)
  ,fGainTableID(0)
  ,fFEEArr(0)
  ,fPTRArr(0)
  ,fGTUArr(0)
{
  //
  // AliTRDCalDCS copy constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCS &AliTRDCalDCS::operator=(const AliTRDCalDCS &cd)
{
  //
  // Assignment operator
  //
  if (&cd == this) return *this;

  new (this) AliTRDCalDCS(cd);
  return *this;
}

