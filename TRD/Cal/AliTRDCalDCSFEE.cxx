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

// fStatusBit:
// 0: no errors
// 1: invalid data received by online fxsproxy
// 2: ROC was not in state CONFIGURED or STANDBY_INIT (most probably it was OFF or STANDBY)
// 3: expected and received DCS-ID do not match. This is a serious communication error!
// 4: DCS id from XML attributes <DCS> and <ack> and the one calculated from SM, S, L do not match
//    This should not happen since the DNR flag should have set the status bit to 0 or 1

#include "AliTRDCalDCSFEE.h"

ClassImp(AliTRDCalDCSFEE)
  
//_____________________________________________________________________________
AliTRDCalDCSFEE::AliTRDCalDCSFEE()
  :TNamed()
  ,fStatusBit(0)
  ,fDCSID(0)
  ,fSM(0)
  ,fStack(0)
  ,fLayer(0)
  ,fNumberOfTimeBins(0)
  ,fPedestal(0)
  ,fConfigTag(0)
  ,fSingleHitThres(0)
  ,fThrPdClsThres(0)
  ,fSelNoZS(0)
  ,fFastStatNoise(0)
  ,fTCFilterWeight(0)
  ,fTCFilterShortDecPar(0)
  ,fTCFilterLongDecPar(0)
  ,fFilterType(0)
  ,fReadoutParam(0)
  ,fTestPattern(0)
  ,fTrackletMode(0)
  ,fTrackletDef(0)
  ,fTriggerSetup(0)
  ,fAddOptions(0) 
  ,fConfigName(0)
  ,fConfigVersion(0)
  ,fGainTableID(0)
{
  //
  // AliTRDCalDCSFEE default constructor
  //
}


//_____________________________________________________________________________
AliTRDCalDCSFEE::AliTRDCalDCSFEE(const char *name, const char *title)
  :TNamed(name,title)
  ,fStatusBit(0)
  ,fDCSID(0)
  ,fSM(0)
  ,fStack(0)
  ,fLayer(0)
  ,fNumberOfTimeBins(0)
  ,fPedestal(0)
  ,fConfigTag(0)
  ,fSingleHitThres(0)
  ,fThrPdClsThres(0)
  ,fSelNoZS(0)
  ,fFastStatNoise(0)
  ,fTCFilterWeight(0)
  ,fTCFilterShortDecPar(0)
  ,fTCFilterLongDecPar(0)
  ,fFilterType(0)
  ,fReadoutParam(0)
  ,fTestPattern(0)
  ,fTrackletMode(0)
  ,fTrackletDef(0)
  ,fTriggerSetup(0)
  ,fAddOptions(0) 
  ,fConfigName(0)
  ,fConfigVersion(0)
  ,fGainTableID(0)
{
  //
  // AliTRDCalDCSFEE constructor
  //
}

