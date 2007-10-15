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
  ~AliHLTPHOSChannelCounter();

  void CountChannels(AliHLTPHOSRcuCellEnergyDataStruct*);
  void PrintOutOfSyncChannels(Int_t);
  void FillHistograms(Int_t);
  void WriteHistograms(const char*);

private: 
  
  UInt_t fChannelArrayPtr[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];
  TH2I *fHistHighGainPtr;
  TH2I *fHistLowGainPtr;
  TH2F *fHistHighRatioPtr;
  TH2F *fHistLowRatioPtr;

  //  ClassDef(AliHLTPHOSChannelCounter, 1);
};

#endif
