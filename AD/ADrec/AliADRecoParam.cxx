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
 
#include "AliLog.h"

#include "AliADRecoParam.h"

ClassImp(AliADRecoParam)

//_____________________________________________________________________________
AliADRecoParam::AliADRecoParam() : AliDetectorRecoParam(),
  fNSigmaPed(2.0),
  fStartClock(8),
  fEndClock(12),
  fNPreClocks(1),
  fNPostClocks(2),
  fTailBegin(16),
  fTailEnd(20),
  fAdcThresHold(1.0),
  fTimeWindowBBALow(-9.5),
  fTimeWindowBBAUp(22.5),
  fTimeWindowBGALow(-2.5),
  fTimeWindowBGAUp(5.0),
  fTimeWindowBBCLow(-2.5),
  fTimeWindowBBCUp(22.5),
  fTimeWindowBGCLow(-2.5),
  fTimeWindowBGCUp(2.5),
  fMaxResid(5.0),
  fResidRise(0.02)		

{
  //
  // constructor
  //
  SetName("AD");
  SetTitle("AD");
}

//_____________________________________________________________________________
AliADRecoParam::~AliADRecoParam() 
{
  //
  // destructor
  //  
}
