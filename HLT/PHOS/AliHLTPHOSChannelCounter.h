
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
#ifndef ALIHLTPHOSCHANNELCOUNTER_H
#define ALIHLTPHOSCHANNELCOUNTER_H

#include "AliHLTPHOSBase.h"

class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSConstants;
class TH2I;
class TH2F;
//class AliHLTPHOSBase;

using namespace PhosHLTConst;

class AliHLTPHOSChannelCounter : public AliHLTPHOSBase
{
public:
  AliHLTPHOSChannelCounter();
  virtual ~AliHLTPHOSChannelCounter();

  void CountChannels(AliHLTPHOSRcuCellEnergyDataStruct* cellEnergy);
  void PrintOutOfSyncChannels(Int_t nEvents);
  void FillHistograms(Int_t nEvents);
  void WriteHistograms(const char* name);

private: 
  
  UInt_t fChannelArrayPtr[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS]; //comment
  TH2I *fHistHighGainPtr; //comment
  TH2I *fHistLowGainPtr; //comment
  TH2F *fHistHighRatioPtr; //comment
  TH2F *fHistLowRatioPtr; //comment

  //  ClassDef(AliHLTPHOSChannelCounter, 1);
};

#endif
