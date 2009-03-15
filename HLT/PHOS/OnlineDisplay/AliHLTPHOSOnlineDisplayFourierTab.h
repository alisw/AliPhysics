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
  AliHLTPHOSOnlineDisplayFourierTab(AliHLTPHOSOnlineDisplay *onlineDisplayPtr, TGTab *tabPtr, AliHLTHOMERReader *fgHomerReaderPtr, AliHLTHOMERReader *fgHomerReadersPtr[MAXHOSTS], int nHosts);
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
  TCanvas *fgCanvasPtr[NGAINS];
  AliHLTPHOSOnlineDisplayTH2D *fgLegoPlotPtr[NGAINS];

  TH1D *fFourierHistoNew[NGAINS];
  TH1D *fFourierHistoOld[NGAINS];
  TH1D *fFourierHistoAccumulated[NGAINS];

  // TRootEmbeddedCanvas *fFourierHistoAccumulatedEC[NGAINS];
  // TRootEmbeddedCanvas *fFourierHistoOldEC[NGAINS];
  // TRootEmbeddedCanvas *fFourierHistoAccumulatedEC[NGAINS];

  //  int *fChannelData[NMODULES][NXRCUCOORD][NZRCUCOORD][NXCOLUMNSRCU][NZROWSRCU][NGAINS];
  // Int_t fNChannelSamples[NMODULES][NXRCUCOORD][NZRCUCOORD][NXCOLUMNSRCU][NZROWSRCU][NGAINS];
  // Int_t fChannelEnergy[NMODULES][NXRCUCOORD][NZRCUCOORD][NXCOLUMNSRCU][NZROWSRCU][NGAINS];
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
