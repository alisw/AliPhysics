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
  fNTdcTimeBinsFlag(410),
  fTdcTimeMinBGFlag(50.0),
  fTdcTimeMaxBGFlag(90.039062),
  fTdcTimeMinBBFlag(170.0),
  fTdcTimeMaxBBFlag(210.039062),
  fNTdcTimeRatioBins(300),
  fTdcTimeRatioMin(1),
  fTdcTimeRatioMax(301),
  fNTdcWidthBins(153),
  fTdcWidthMin(2.343750),
  fTdcWidthMax(121.875000),
  fNChargeChannelBins(1000),
  fChargeChannelMin(0),
  fChargeChannelMax(1000),
  fChargeChannelZoomMin(0),
  fChargeChannelZoomMax(50),
  fNChargeSideBins(500),
  fChargeSideMin(1),
  fChargeSideMax(5000),
  fNChargeCorrBins(101),
  fChargeCorrMin(1),
  fChargeCorrMax(1001),
  fNPairTimeCorrBins(306), 
  fPairTimeCorrMin(0.976562),
  fPairTimeCorrMax(299.804688),
  fNPairTimeDiffBins(154), 
  fPairTimeDiffMin(-15.039062),
  fPairTimeDiffMax(15.039062),
  fNMeanTimeCorrBins(306), 
  fMeanTimeCorrMin(50.976562),
  fMeanTimeCorrMax(349.804688),
  fChargeTrendMin(0),
  fChargeTrendMax(1000),
  fSatMed(0.1),
  fSatHigh(0.3),
  fSatHuge(0.5),
  fMaxPedDiff(1),
  fMaxPedWidth(1.5),
  fMaxNoTimeRate(10e-4),
  fMaxNoFlagRate(10e-3),
  fMaxBBVariation(0.1),
  fMaxBGVariation(0.1),
  fAsynchronBB(0.5),
  fAsynchronBG(0.5)
  
   	

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
