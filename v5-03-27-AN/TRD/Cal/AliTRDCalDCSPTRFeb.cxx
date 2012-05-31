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

/* $Id: AliTRDCalDCSPTRFeb.cxx 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS GTU parameters                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDCSPTRFeb.h"

ClassImp(AliTRDCalDCSPTRFeb)

//_____________________________________________________________________________
AliTRDCalDCSPTRFeb::AliTRDCalDCSPTRFeb()
  :TNamed()
    ,fSide("")
    ,fDetName("")
    ,fPrimary(0)
    ,fPtChDelay0(0)
    ,fPtChDelay1(0)
    ,fPtChDelay2(0)
    ,fPtChDelay3(0)
    ,fPtChDelay4(0)
    ,fPtChDelay5(0)
    ,fPtChDelay6(0)
    ,fPtChDelay7(0)
    ,fPtChDelay8(0)
    ,fPtChDelay9(0)
    ,fPtChDelay10(0)
    ,fPtChDelay11(0)
    ,fPtChThr0(0)
    ,fPtChThr1(0)
    ,fPtChThr2(0)
    ,fPtChThr3(0)
    ,fPtChThr4(0)
    ,fPtChThr5(0)
    ,fPtChThr6(0)
    ,fPtChThr7(0)
    ,fPtChThr8(0)
    ,fPtChThr9(0)
    ,fPtChThr10(0)
    ,fPtChThr11(0)
    ,fClkLb(0)
    ,fClkHb(0)
    ,fCh0CountLb(0)
    ,fCh0CountHb(0)
    ,fCh1CountLb(0)
    ,fCh1CountHb(0)
    ,fCh2CountLb(0)
    ,fCh2CountHb(0)
    ,fCh3CountLb(0)
    ,fCh3CountHb(0)
    ,fCh4CountLb(0)
    ,fCh4CountHb(0)
    ,fCh5CountLb(0)
    ,fCh5CountHb(0)
    ,fCh6CountLb(0)
    ,fCh6CountHb(0)
    ,fCh7CountLb(0)
    ,fCh7CountHb(0)
    ,fCh8CountLb(0)
    ,fCh8CountHb(0)
    ,fCh9CountLb(0)
    ,fCh9CountHb(0)
    ,fCh10CountLb(0)
    ,fCh10CountHb(0)
    ,fCh11CountLb(0)
    ,fCh11CountHb(0)
    ,fTrigParallel0Lb(0)
    ,fTrigParallel0Hb(0)
    ,fTrigParallel1Lb(0)
    ,fTrigParallel1Hb(0)
    
{
  //
  // AliTRDCalDCSGTU default constructor
  //
  fSegmentsArr->SetOwner();
}

//_____________________________________________________________________________
AliTRDCalDCSPTRFeb::AliTRDCalDCSPTRFeb(const char *name, const char *title)
  :TNamed(name,title)
    ,fSide("")
    ,fDetName("")
    ,fPrimary(0)
    ,fPtChDelay0(0)
    ,fPtChDelay1(0)
    ,fPtChDelay2(0)
    ,fPtChDelay3(0)
    ,fPtChDelay4(0)
    ,fPtChDelay5(0)
    ,fPtChDelay6(0)
    ,fPtChDelay7(0)
    ,fPtChDelay8(0)
    ,fPtChDelay9(0)
    ,fPtChDelay10(0)
    ,fPtChDelay11(0)
    ,fPtChThr0(0)
    ,fPtChThr1(0)
    ,fPtChThr2(0)
    ,fPtChThr3(0)
    ,fPtChThr4(0)
    ,fPtChThr5(0)
    ,fPtChThr6(0)
    ,fPtChThr7(0)
    ,fPtChThr8(0)
    ,fPtChThr9(0)
    ,fPtChThr10(0)
    ,fPtChThr11(0)
    ,fClkLb(0)
    ,fClkHb(0)
    ,fCh0CountLb(0)
    ,fCh0CountHb(0)
    ,fCh1CountLb(0)
    ,fCh1CountHb(0)
    ,fCh2CountLb(0)
    ,fCh2CountHb(0)
    ,fCh3CountLb(0)
    ,fCh3CountHb(0)
    ,fCh4CountLb(0)
    ,fCh4CountHb(0)
    ,fCh5CountLb(0)
    ,fCh5CountHb(0)
    ,fCh6CountLb(0)
    ,fCh6CountHb(0)
    ,fCh7CountLb(0)
    ,fCh7CountHb(0)
    ,fCh8CountLb(0)
    ,fCh8CountHb(0)
    ,fCh9CountLb(0)
    ,fCh9CountHb(0)
    ,fCh10CountLb(0)
    ,fCh10CountHb(0)
    ,fCh11CountLb(0)
    ,fCh11CountHb(0)
    ,fTrigParallel0Lb(0)
    ,fTrigParallel0Hb(0)
    ,fTrigParallel1Lb(0)
    ,fTrigParallel1Hb(0)
{
  //
  // AliTRDCalDCSGTU constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCSPTRFeb::AliTRDCalDCSPTRFeb(const AliTRDCalDCSPTRFeb &)
  :TNamed("","")
    ,fSide("")
    ,fDetName("")
    ,fPrimary(0)
    ,fPtChDelay0(0)
    ,fPtChDelay1(0)
    ,fPtChDelay2(0)
    ,fPtChDelay3(0)
    ,fPtChDelay4(0)
    ,fPtChDelay5(0)
    ,fPtChDelay6(0)
    ,fPtChDelay7(0)
    ,fPtChDelay8(0)
    ,fPtChDelay9(0)
    ,fPtChDelay10(0)
    ,fPtChDelay11(0)
    ,fPtChThr0(0)
    ,fPtChThr1(0)
    ,fPtChThr2(0)
    ,fPtChThr3(0)
    ,fPtChThr4(0)
    ,fPtChThr5(0)
    ,fPtChThr6(0)
    ,fPtChThr7(0)
    ,fPtChThr8(0)
    ,fPtChThr9(0)
    ,fPtChThr10(0)
    ,fPtChThr11(0)
    ,fClkLb(0)
    ,fClkHb(0)
    ,fCh0CountLb(0)
    ,fCh0CountHb(0)
    ,fCh1CountLb(0)
    ,fCh1CountHb(0)
    ,fCh2CountLb(0)
    ,fCh2CountHb(0)
    ,fCh3CountLb(0)
    ,fCh3CountHb(0)
    ,fCh4CountLb(0)
    ,fCh4CountHb(0)
    ,fCh5CountLb(0)
    ,fCh5CountHb(0)
    ,fCh6CountLb(0)
    ,fCh6CountHb(0)
    ,fCh7CountLb(0)
    ,fCh7CountHb(0)
    ,fCh8CountLb(0)
    ,fCh8CountHb(0)
    ,fCh9CountLb(0)
    ,fCh9CountHb(0)
    ,fCh10CountLb(0)
    ,fCh10CountHb(0)
    ,fCh11CountLb(0)
    ,fCh11CountHb(0)
    ,fTrigParallel0Lb(0)
    ,fTrigParallel0Hb(0)
    ,fTrigParallel1Lb(0)
    ,fTrigParallel1Hb(0)
{
  //
  // AliTRDCalDCSGTU constructor
  //
}



