//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSONLINEDISPLAYFOURIERTAB_H
#define ALIHLTPHOSONLINEDISPLAYFOURIERTAB_H

#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include "AliHLTPHOSOnlineDisplayTab.h"

#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>
#include "AliHLTPHOSOnlineDisplayTH2D.h"
#include "AliHLTPHOSConstants.h"
//#include "AliHLTPHOSFourier.h"
#include "AliHLTPHOSRcuFFTDataStruct.h"
#define N_ZRCU_COORD 2
#define N_XRCU_COORD 2

using namespace PhosHLTConst;

class AliHLTPHOSGetEventButton;
class AliHLTHOMERReader;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSOnlineDisplay;
class AliHLTPHOSSharedMemoryInterface;
class AliHLTPHOSFourier;

class AliHLTPHOSOnlineDisplayFourierTab : public AliHLTPHOSOnlineDisplayTab
{
 public:
  virtual ~AliHLTPHOSOnlineDisplayFourierTab();
  AliHLTPHOSOnlineDisplayFourierTab(AliHLTPHOSOnlineDisplay *onlineDisplayPtr, TGTab *tabPtr, AliHLTHOMERReader *fgHomerReaderPtr, AliHLTHOMERReader *fgHomerReadersPtr[MAX_HOSTS], int nHosts);
  Int_t GetRawData(TH1D *histPtr, int x, int z, int gain);
  void UpdateDisplay();
  int GetNextEvent();
  virtual void ReadBlockData(AliHLTHOMERReader *homeReaderPtr);
  void FindFourierBlocks(AliHLTHOMERReader *homeReaderPtr);

  void ResetDisplay();
  TGTab               *fTab;
  TGTab               *fSubTab1;
  TRootEmbeddedCanvas *fEc1, *fEc2, *fEc3, *fEc4, *fEc5, *fEc6;
  //  TRootEmbeddedCanvas *fEc1, *fEc2, *fEc3, *fEc4;
  //  TGCompositeFrame    *fSubF1, *fSubF2;
  TGCompositeFrame    *fSubF1, *fSubF2, *fSubF3;
  TCanvas *fgCanvasPtr[N_GAINS];
  AliHLTPHOSOnlineDisplayTH2D *fgLegoPlotPtr[N_GAINS];

  TH1D *fFourierHistoNew[N_GAINS];
  TH1D *fFourierHistoOld[N_GAINS];
  TH1D *fFourierHistoAccumulated[N_GAINS];

  // TRootEmbeddedCanvas *fFourierHistoAccumulatedEC[N_GAINS];
  // TRootEmbeddedCanvas *fFourierHistoOldEC[N_GAINS];
  // TRootEmbeddedCanvas *fFourierHistoAccumulatedEC[N_GAINS];

  //  int *fChannelData[N_MODULES][N_XRCU_COORD][N_ZRCU_COORD][N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  // Int_t fNChannelSamples[N_MODULES][N_XRCU_COORD][N_ZRCU_COORD][N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  // Int_t fChannelEnergy[N_MODULES][N_XRCU_COORD][N_ZRCU_COORD][N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  const  char* Gain2Text(const int gain, const char delimeter);

 protected:
  Bool_t fgAccumulate;

 private:
  void FillHistograms(const AliHLTPHOSRcuFFTDataStruct psd, const int size);

  AliHLTPHOSOnlineDisplayFourierTab();
  AliHLTPHOSGetEventButton* fgEventButtPtr; 
  void InitDisplay(TGTab *tabPtr);
  AliHLTPHOSOnlineDisplay *fOnlineDisplayPtr;
  AliHLTPHOSSharedMemoryInterface *fShmPtr;   

  AliHLTPHOSFourier *fFourierPtr;
  char fGainText[256];

  unsigned long fEvtCnt;
};


#endif
