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

/* $Id: AliTRDCalDCSPTRCba.cxx 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS GTU parameters                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDCSPTRCba.h"

ClassImp(AliTRDCalDCSPTRCba)

//_____________________________________________________________________________
AliTRDCalDCSPTRCba::AliTRDCalDCSPTRCba()
  :TNamed()
    ,fSide("")
    ,fPrimary(0)
    ,fChDelayT0(0)
    ,fChDelayV0(0)
    ,fChDelayV1(0)
    ,fChDelayV2(0)
    ,fChDelayV3(0)
    ,fChDisableT0(0)
    ,fChDisableV0(0)
    ,fChDisableV1(0)
    ,fChDisableV2(0)
    ,fChDisableV3(0)
    ,fTo27ParralelLb(0)
    ,fTo27ParralelHb(0)
    ,fTo28ParralelLb(0)
    ,fTo28ParralelHb(0)
    ,fTo29ParralelLb(0)
    ,fTo29ParralelHb(0)
    ,fTo30ParralelLb(0)
    ,fTo30ParralelHb(0)
    ,fTo31ParralelLb(0)
    ,fTo31ParralelHb(0)
    ,fTo32ParralelLb(0)
    ,fTo32ParralelHb(0)
    ,fTo33ParralelLb(0)
    ,fTo33ParralelHb(0)
    ,fTo34ParralelLb(0)
    ,fTo34ParralelHb(0)
    ,fTo35ParralelLb(0)
    ,fTo35ParralelHb(0)
    ,fTo36ParralelLb(0)
    ,fTo36ParralelHb(0)
    ,fClkLb(0)
    ,fClkHb(0)
    ,fBitsToCbB42Lb(0)
    ,fBitsToCbB42Hb(0)
    ,fBitsToCbB43Lb(0)
    ,fBitsToCbB43Hb(0)
    ,fBitsToCbB44Lb(0)
    ,fBitsToCbB44Hb(0)
    ,fBitsToCbB45Lb(0)
    ,fBitsToCbB45Hb(0)
    
{
  //
  // AliTRDCalDCSGTU default constructor
  //
  fSegmentsArr->SetOwner();
}

//_____________________________________________________________________________
AliTRDCalDCSPTRCba::AliTRDCalDCSPTRCba(const char *name, const char *title)
  :TNamed(name,title)
    ,fSide("")
    ,fPrimary(0)
    ,fChDelayT0(0)
    ,fChDelayV0(0)
    ,fChDelayV1(0)
    ,fChDelayV2(0)
    ,fChDelayV3(0)
    ,fChDisableT0(0)
    ,fChDisableV0(0)
    ,fChDisableV1(0)
    ,fChDisableV2(0)
    ,fChDisableV3(0)
    ,fTo27ParralelLb(0)
    ,fTo27ParralelHb(0)
    ,fTo28ParralelLb(0)
    ,fTo28ParralelHb(0)
    ,fTo29ParralelLb(0)
    ,fTo29ParralelHb(0)
    ,fTo30ParralelLb(0)
    ,fTo30ParralelHb(0)
    ,fTo31ParralelLb(0)
    ,fTo31ParralelHb(0)
    ,fTo32ParralelLb(0)
    ,fTo32ParralelHb(0)
    ,fTo33ParralelLb(0)
    ,fTo33ParralelHb(0)
    ,fTo34ParralelLb(0)
    ,fTo34ParralelHb(0)
    ,fTo35ParralelLb(0)
    ,fTo35ParralelHb(0)
    ,fTo36ParralelLb(0)
    ,fTo36ParralelHb(0)
    ,fClkLb(0)
    ,fClkHb(0)
    ,fBitsToCbB42Lb(0)
    ,fBitsToCbB42Hb(0)
    ,fBitsToCbB43Lb(0)
    ,fBitsToCbB43Hb(0)
    ,fBitsToCbB44Lb(0)
    ,fBitsToCbB44Hb(0)
    ,fBitsToCbB45Lb(0)
    ,fBitsToCbB45Hb(0)
    
{
  //
  // AliTRDCalDCSGTU constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCSPTRCba::AliTRDCalDCSPTRCba(const AliTRDCalDCSPTRCba &)
  :TNamed("","")
    ,fSide("")
    ,fPrimary(0)
    ,fChDelayT0(0)
    ,fChDelayV0(0)
    ,fChDelayV1(0)
    ,fChDelayV2(0)
    ,fChDelayV3(0)
    ,fChDisableT0(0)
    ,fChDisableV0(0)
    ,fChDisableV1(0)
    ,fChDisableV2(0)
    ,fChDisableV3(0)
    ,fTo27ParralelLb(0)
    ,fTo27ParralelHb(0)
    ,fTo28ParralelLb(0)
    ,fTo28ParralelHb(0)
    ,fTo29ParralelLb(0)
    ,fTo29ParralelHb(0)
    ,fTo30ParralelLb(0)
    ,fTo30ParralelHb(0)
    ,fTo31ParralelLb(0)
    ,fTo31ParralelHb(0)
    ,fTo32ParralelLb(0)
    ,fTo32ParralelHb(0)
    ,fTo33ParralelLb(0)
    ,fTo33ParralelHb(0)
    ,fTo34ParralelLb(0)
    ,fTo34ParralelHb(0)
    ,fTo35ParralelLb(0)
    ,fTo35ParralelHb(0)
    ,fTo36ParralelLb(0)
    ,fTo36ParralelHb(0)
    ,fClkLb(0)
    ,fClkHb(0)
    ,fBitsToCbB42Lb(0)
    ,fBitsToCbB42Hb(0)
    ,fBitsToCbB43Lb(0)
    ,fBitsToCbB43Hb(0)
    ,fBitsToCbB44Lb(0)
    ,fBitsToCbB44Hb(0)
    ,fBitsToCbB45Lb(0)
    ,fBitsToCbB45Hb(0)
    
{
  //
  // AliTRDCalDCSGTU constructor
  //
}



