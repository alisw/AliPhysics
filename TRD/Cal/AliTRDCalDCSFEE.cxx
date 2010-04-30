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
// 0: no errors for that ROC
// 1: ROC sent invalid or corrupted data. 
// 2: ROC was not in state CONFIGURED or STANDBY_INIT (most probably it was in STANDBY)
// 3: No new data received from that ROC.
// 4: DCS id from XML attributes <DCS> and <ack> and the one calculated from SM, S, L do not match
// 5: ROC has not responded at all, most probably it was off.

#include "AliTRDCalDCSFEE.h"

ClassImp(AliTRDCalDCSFEE)
  
//_____________________________________________________________________________
AliTRDCalDCSFEE::AliTRDCalDCSFEE()
  :TObject()
  ,fStatusBit(0)
  ,fSM(-1)
  ,fStack(-1)
  ,fLayer(-1)
  ,fGainTableRocSerial(0)
  ,fDCSID(-1)
  ,fNumberOfTimeBins(-1)
  ,fConfigTag(-1)
  ,fSingleHitThres(-1)
  ,fThrPdClsThres(-1)
  ,fSelNoZS(-1)
  ,fTCFilterWeight(-1)
  ,fTCFilterShortDecPar(-1)
  ,fTCFilterLongDecPar(-1)
  ,fFastStatNoise(-1)
  ,fGainTableRocType("")
  ,fFilterType("")
  ,fReadoutParam("")
  ,fTestPattern("")
  ,fTrackletMode("")
  ,fTrackletDef("")
  ,fTriggerSetup("")
  ,fAddOptions("") 
  ,fConfigName("")
  ,fConfigVersion("")
  ,fGainTableName("")
  ,fGainTableDesc("")
{
  //
  // AliTRDCalDCSFEE default constructor
  //
  for(Int_t i=0; i<(Int_t)fgkROB; i++) {
    for(Int_t j=0; j<(Int_t)fgkMCM; j++) {
      fRStateGSM[i][j]  = -1;
      fRStateNI[i][j]   = -1;
      fRStateEV[i][j]   = -1;
      fRStatePTRG[i][j] = -1;
      fGainTableAdcdac[i][j] = -1;
      for(Int_t k=0; k<(Int_t)fgkADC; k++) {
	fGainTableFgfn[i][j][k] = -1;
	fGainTableFgan[i][j][k] = -1;
      }
    }
  }
}


//_____________________________________________________________________________
AliTRDCalDCSFEE::AliTRDCalDCSFEE(const AliTRDCalDCSFEE &c)
  :TObject(c)
  ,fStatusBit(c.fStatusBit)
  ,fSM(c.fSM)
  ,fStack(c.fStack)
  ,fLayer(c.fLayer)
  ,fGainTableRocSerial(c.fGainTableRocSerial)
  ,fDCSID(c.fDCSID)
  ,fNumberOfTimeBins(c.fNumberOfTimeBins)
  ,fConfigTag(c.fConfigTag)
  ,fSingleHitThres(c.fSingleHitThres)
  ,fThrPdClsThres(c.fThrPdClsThres)
  ,fSelNoZS(c.fSelNoZS)
  ,fTCFilterWeight(c.fTCFilterWeight)
  ,fTCFilterShortDecPar(c.fTCFilterShortDecPar)
  ,fTCFilterLongDecPar(c.fTCFilterLongDecPar)
  ,fFastStatNoise(c.fFastStatNoise)
  ,fGainTableRocType(c.fGainTableRocType)
  ,fFilterType(c.fFilterType)
  ,fReadoutParam(c.fReadoutParam)
  ,fTestPattern(c.fTestPattern)
  ,fTrackletMode(c.fTrackletMode)
  ,fTrackletDef(c.fTrackletDef)
  ,fTriggerSetup(c.fTriggerSetup)
  ,fAddOptions(c.fAddOptions) 
  ,fConfigName(c.fConfigName)
  ,fConfigVersion(c.fConfigVersion)
  ,fGainTableName(c.fGainTableName)
  ,fGainTableDesc(c.fGainTableDesc)
{
  //
  // AliTRDCalDCSFEE copy constructor
  //
  for(Int_t i=0; i<(Int_t)fgkROB; i++) {
    for(Int_t j=0; j<(Int_t)fgkMCM; j++) {
      fRStateGSM[i][j]  = c.fRStateGSM[i][j];
      fRStateNI[i][j]   = c.fRStateNI[i][j];
      fRStateEV[i][j]   = c.fRStateEV[i][j];
      fRStatePTRG[i][j] = c.fRStatePTRG[i][j];
      fGainTableAdcdac[i][j] = c.fGainTableAdcdac[i][j];
      for(Int_t k=0; k<(Int_t)fgkADC; k++) {
	fGainTableFgfn[i][j][k] = c.fGainTableFgfn[i][j][k];
	fGainTableFgan[i][j][k] = c.fGainTableFgan[i][j][k];
      }
    }
  }
}


//_____________________________________________________________________________
AliTRDCalDCSFEE &AliTRDCalDCSFEE::operator=(const AliTRDCalDCSFEE &c)
{
  //
  // Assignment operator
  //
  if (&c == this) return *this;

  new (this) AliTRDCalDCSFEE(c);
  return *this;
}

