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
 
#include "AliADQAParam.h"

ClassImp(AliADQAParam)

//_____________________________________________________________________________
AliADQAParam::AliADQAParam():
  fNTdcTimeBins(3062),
  fTdcTimeMin(0.976562),
  fTdcTimeMax(300.0),
  fNTdcWidthBins(153),
  fTdcWidthMin(2.343750),
  fTdcWidthMax(121.875000),
  fNChargeChannelBins(1000),
  fChargeChannelMin(1),
  fChargeChannelMax(1001),
  fNChargeSideBins(500),
  fChargeSideMin(1),
  fChargeSideMax(5000),
  fNChargeCorrBins(101),
  fChargeCorrMin(1),
  fChargeCorrMax(1001),
  fNPairTimeCorrBins(614), 
  fPairTimeCorrMin(70.019531),
  fPairTimeCorrMax(129.980469),
  fNPairTimeDiffBins(154), 
  fPairTimeDiffMin(-15.039062),
  fPairTimeDiffMax(15.039062),
  fNMeanTimeCorrBins(614), 
  fMeanTimeCorrMin(70.019531),
  fMeanTimeCorrMax(129.980469),
  fSatMed(0.1),
  fSatHigh(0.3),
  fSatHuge(0.5),
  fMaxPedDiff(1) 	

{
  //
  // constructor
  //
}

//_____________________________________________________________________________
AliADQAParam::~AliADQAParam() 
{
  //
  // destructor
  //  
}
