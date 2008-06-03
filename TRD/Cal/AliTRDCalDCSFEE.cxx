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

/* $Id: AliTRDCalDCSFEE.cxx 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS FEE configuration parameters           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDCSFEE.h"

ClassImp(AliTRDCalDCSFEE)

//_____________________________________________________________________________
AliTRDCalDCSFEE::AliTRDCalDCSFEE()
  :TNamed()
  ,fDCSID(0)
  ,fSM(0)
  ,fStack(0)
  ,fLayer(0)
  ,fNumberOfTimeBins(0)
  ,fTailCancelationTau1(0)
  ,fTailCancelationTau2(0)
  ,fTailCancelationAmp(0)
  ,fPedestal(0)
  ,fConfigID(0)
  ,fGainTableID(0)
{
  //
  // AliTRDCalDCSFEE default constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCSFEE::AliTRDCalDCSFEE(const char *name, const char *title)
  :TNamed(name,title)
  ,fDCSID(0)
  ,fSM(0)
  ,fStack(0)
  ,fLayer(0)
  ,fNumberOfTimeBins(0)
  ,fTailCancelationTau1(0)
  ,fTailCancelationTau2(0)
  ,fTailCancelationAmp(0)
  ,fPedestal(0)
  ,fConfigID(0)
  ,fGainTableID(0)
{
  //
  // AliTRDCalDCSFEE constructor
  //
}

