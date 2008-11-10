// $Id$

 /**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** 
 * @file   AliHLTPHOSClusterizer.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Clusterizer for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSEmcCalibrationHistogramProducer.h"
#include "AliHLTPHOSBase.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TH2F.h"

AliHLTPHOSEmcCalibrationHistogramProducer::AliHLTPHOSEmcCalibrationHistogramProducer() :
  AliHLTPHOSBase()
  //fChannelEnergyHistogramPtr(0),
  //fMeanEnergyHistogramPtr(0),
  //fMaxSignalHistogramPtr(0),
  //fNHitsHistogramPtr(0),
  //fLiveChannelHistogramPtr(0),
  // fBadChannelHistogramPtr(0)
{
  //See header file for documentation
}

AliHLTPHOSEmcCalibrationHistogramProducer::~AliHLTPHOSEmcCalibrationHistogramProducer()
{
  //See header file for documentation
}

void
AliHLTPHOSEmcCalibrationHistogramProducer::Reset()
{
  //See header file for documentation
}

void
AliHLTPHOSEmcCalibrationHistogramProducer::FillChannelEnergy(Int_t module, Int_t xCol, Int_t zRow, Int_t gain)
{
  //See header file for documentation
}

void
AliHLTPHOSEmcCalibrationHistogramProducer::FillRCUEnergies(TH1F* rcuChannelEnergyPtr[][N_ZROWS_RCU][N_GAINS], Int_t module, Int_t rcuX, Int_t rcuZ)
{
  //See header file for documentation

  Char_t histName[128];

  for(UInt_t x = 0; x < N_XCOLUMNS_RCU; x++)
    {
      for(UInt_t z = 0; z < N_ZROWS_RCU; z++)
	{
	  for(UInt_t gain = 0; gain < N_GAINS; gain++)
	    {
	      sprintf(histName, "mod%dcol%drow%dgain%d", module, x, z, gain);
	      *fChannelEnergyHistogramPtr[module][x+rcuX*N_XCOLUMNS_RCU][z+rcuZ*N_ZROWS_RCU][gain] = (TH1F*)rcuChannelEnergyPtr[x][z][gain]->Clone(histName);
	    }
	}
    }
}

void
AliHLTPHOSEmcCalibrationHistogramProducer::FillBadChannels(TH2F* rcuBadChannelHistPtr[], Int_t module)
{
  //See header file for documentation
  //TODO!
}

TH1F* 
AliHLTPHOSEmcCalibrationHistogramProducer::GetChannelEnergyHistogram(Int_t module, Int_t xCol, Int_t zRow, Int_t gain)
{
  //See header file for documentation
  return fChannelEnergyHistogramPtr[module][xCol][zRow][gain];
}

TH2F*
AliHLTPHOSEmcCalibrationHistogramProducer::GetMeanEnergyHistogram(Int_t module, Int_t gain)
{
  //See header file for documentation
  for(UInt_t x = 0; x < N_XCOLUMNS_RCU; x++)
    {
      for(UInt_t z = 0; z < N_ZROWS_RCU; z++)
	{
	  fMeanEnergyHistogramPtr->SetBinContent(x, z, fChannelEnergyHistogramPtr[module][x][z][gain]->GetMean());
	}
    }
}

TH2F*
AliHLTPHOSEmcCalibrationHistogramProducer::GetMaxSignalHistogram(Int_t module, Int_t gain)
{
  //See header file for documentation
  for(UInt_t x = 0; x < N_XCOLUMNS_RCU; x++)
    {
      for(UInt_t z = 0; z < N_ZROWS_RCU; z++)
	{
	  fMaxSignalHistogramPtr->SetBinContent(x, z, fChannelEnergyHistogramPtr[module][x][z][gain]->GetMaximum());
	}
    }
}

TH2I*
AliHLTPHOSEmcCalibrationHistogramProducer::GetNHitsHistogram(Int_t module, Int_t gain)
{
  //See header file for documentation
  for(UInt_t x = 0; x < N_XCOLUMNS_RCU; x++)
    {
      for(UInt_t z = 0; z < N_ZROWS_RCU; z++)
	{
	  // TODO!
	}
    }
}

TH2I*
AliHLTPHOSEmcCalibrationHistogramProducer::GetLiveChannelHistogram(Int_t module, Int_t gain, Float_t minSignal = 30)
{
  //See header file for documentation
  //TODO!
}

TH2I*
AliHLTPHOSEmcCalibrationHistogramProducer::GetBadChannelHistogram(Int_t module, Int_t gain)
{
  //See header file for documentation
  //TODO!
}

