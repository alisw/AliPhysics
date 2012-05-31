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

#include "AliTPCeventInfo.h"

/*
   Placeholder for the TPC event specific information
   

*/

ClassImp(AliTPCeventInfo)

AliTPCeventInfo::AliTPCeventInfo()
  :fNSamples(0),
   fSamplingRate(0),
   fSamplingPhase(0),
   fBC(kFALSE),
   fBCMode1(0),
   fBCMode2(0),
   fBCPre(0),
   fBCPost(0),
   fNPreTrigger(0),
   fTC(kFALSE),
   fZS(kFALSE),
   fZSGlitchFilter(0),
   fZSPre(0),
   fZSPost(0),
   fZSOffset(0) {
  //
  // Constructor
  //
}

AliTPCeventInfo::~AliTPCeventInfo() {
  //
  // destructor
  //
}
