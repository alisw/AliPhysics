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

/* $Id: AliTRDCalDCSGTUSegment.cxx 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS GTU parameters                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDCSGTUSegment.h"
#include <TObjArray.h>

ClassImp(AliTRDCalDCSGTUSegment)


//_____________________________________________________________________________
AliTRDCalDCSGTUSegment::AliTRDCalDCSGTUSegment()
  :TNamed()
    ,fId(0)
    ,fFromRunNumber(0)
    ,fFromSORFlag(0)
    ,fChild(0)
    ,fTmuArr(new TObjArray(5))
    ,fSmuStackMask(0)
    ,fSmuTracklets(0)
    ,fSmuTracks(0)
    ,fSmuIdelay(0)
    ,fSmuTriggerWindowL1Low(0)
    ,fSmuTriggerWindowL1High(0)
    ,fSmuTriggerWindowL2Low(0)
    ,fSmuTriggerWindowL2High(0)
    ,fSmuTtcEmulatorEnable(0)
    ,fSmuBoardInfo()
{
  //
  // AliTRDCalDCSGTU default constructor
  //
  fTmuArr->SetOwner();


}

//_____________________________________________________________________________
AliTRDCalDCSGTUSegment::AliTRDCalDCSGTUSegment(const char *name, const char *title)
  :TNamed(name,title)
    ,fId(0)
    ,fFromRunNumber(0)
    ,fFromSORFlag(0)
    ,fChild(0)
    ,fTmuArr(new TObjArray(5))
    ,fSmuStackMask(0)
    ,fSmuTracklets(0)
    ,fSmuTracks(0)
    ,fSmuIdelay(0)
    ,fSmuTriggerWindowL1Low(0)
    ,fSmuTriggerWindowL1High(0)
    ,fSmuTriggerWindowL2Low(0)
    ,fSmuTriggerWindowL2High(0)
    ,fSmuTtcEmulatorEnable(0)
    ,fSmuBoardInfo()
{
  //
  // AliTRDCalDCSGTU constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCSGTUSegment::AliTRDCalDCSGTUSegment(const AliTRDCalDCSGTUSegment&)
  :TNamed("","")
    ,fId(0)
    ,fFromRunNumber(0)
    ,fFromSORFlag(0)
    ,fChild(0)
    ,fTmuArr(new TObjArray(5))
    ,fSmuStackMask(0)
    ,fSmuTracklets(0)
    ,fSmuTracks(0)
    ,fSmuIdelay(0)
    ,fSmuTriggerWindowL1Low(0)
    ,fSmuTriggerWindowL1High(0)
    ,fSmuTriggerWindowL2Low(0)
    ,fSmuTriggerWindowL2High(0)
    ,fSmuTtcEmulatorEnable(0)
    ,fSmuBoardInfo()
{
  //
  // AliTRDCalDCSGTUSegment constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCSGTUSegment& AliTRDCalDCSGTUSegment::operator=(const AliTRDCalDCSGTUSegment& sh)
{
  //
  // AliTRDCalDCSGTUSegment constructor
  //
  if (&sh == this) return *this;
  
  new (this) AliTRDCalDCSGTUSegment(sh);
  return *this;
}


