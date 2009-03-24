//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Per Thomas Hille for the ALICE                                *
 * offline/HLT Project. Contributors are mentioned in the code where      *
 * appropriate.                                                           *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSONLINEDISPLAYFOURIERTAB_H
#define ALIHLTPHOSONLINEDISPLAYFOURIERTAB_H

#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include "AliHLTPHOSOnlineDisplayTab.h"

// #include <TCanvas.h>
// #include <TH2D.h>
// #include <TH1D.h>
// #include "AliHLTPHOSOnlineDisplayTH2D.h"
// #include "AliHLTPHOSConstants.h"
// //#include "AliHLTPHOSFourier.h"
// #include "AliHLTPHOSRcuFFTDataStruct.h"
#define NZRCUCOORD 2
#define NXRCUCOORD 2

using namespace PhosHLTConst;

class TH1D;
class TH2D;
class TCanvas;
class TRootEmbeddedCanvas;
class TGTab;
class AliHLTPHOSRcuFFTDataStruct;
class AliHLTPHOSConstants;
class AliHLTPHOSOnlineDisplayTH2D;
class AliHLTPHOSGetEventButton;
class AliHLTHOMERReader;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSOnlineDisplay;
class AliHLTPHOSSharedMemoryInterface;
class AliHLTPHOSFourier;

class AliHLTPHOSOnlineDisplayFourierTab : public AliHLTPHOSOnlineDisplayTab
{
 public:
  virtual ~AliHLTPHOSOnlineDisplayFourierTab(); // destructor 
  AliHLTPHOSOnlineDisplayFourierTab(AliHLTPHOSOnlineDisplay * const onlineDisplayPtr, TGTab *tabPtr, const AliHLTHOMERReader * fgHomerReaderPtr, const AliHLTHOMERReader * const fgHomerReadersPtr[MAXHOSTS], int nHosts); // constructor
  Int_t GetRawData(TH1D *histPtr, int x, int z, int gain); // GetRawData
  void UpdateDisplay(); //UpdateDisplay
  int GetNextEvent(); //GetNextEvent
  virtual void ReadBlockData(AliHLTHOMERReader * const homeReaderPtr); //ReadBlockData
  void FindFourierBlocks(AliHLTHOMERReader *homeReaderPtr);//FindFourierBlocks

  void ResetDisplay(); //ResetDisplay
  TGTab               *fTab;  //fTab
  TGTab               *fSubTab1; //fSubTab1
  TRootEmbeddedCanvas *fEc1, *fEc2, *fEc3, *fEc4, *fEc5, *fEc6; //Canvases
  //  TRootEmbeddedCanvas *fEc1, *fEc2, *fEc3, *fEc4;
  //  TGCompositeFrame    *fSubF1, *fSubF2;
  TGCompositeFrame    *fSubF1, *fSubF2, *fSubF3; //frames
  TCanvas *fgCanvasPtr[NGAINS]; // canvas
  AliHLTPHOSOnlineDisplayTH2D *fgLegoPlotPtr[NGAINS]; //legoplot

  TH1D *fFourierHistoNew[NGAINS]; //histogram
  TH1D *fFourierHistoOld[NGAINS]; //histogram 
  TH1D *fFourierHistoAccumulated[NGAINS]; //histogram

  // TRootEmbeddedCanvas *fFourierHistoAccumulatedEC[NGAINS];
  // TRootEmbeddedCanvas *fFourierHistoOldEC[NGAINS];
  // TRootEmbeddedCanvas *fFourierHistoAccumulatedEC[NGAINS];

  //  int *fChannelData[NMODULES][NXRCUCOORD][NZRCUCOORD][NXCOLUMNSRCU][NZROWSRCU][NGAINS];
  // Int_t fNChannelSamples[NMODULES][NXRCUCOORD][NZRCUCOORD][NXCOLUMNSRCU][NZROWSRCU][NGAINS];
  // Int_t fChannelEnergy[NMODULES][NXRCUCOORD][NZRCUCOORD][NXCOLUMNSRCU][NZROWSRCU][NGAINS];
  const  char* Gain2Text(const int gain, const char delimeter); //gain2text

 protected:
  Bool_t fgAccumulate; //fgAccumulate

 private:
  void FillHistograms(const AliHLTPHOSRcuFFTDataStruct psd, const int size); //FillHistograms

  AliHLTPHOSOnlineDisplayFourierTab(); // default constructor
  AliHLTPHOSGetEventButton* fgEventButtPtr;  // fgEventButtPtr
  void InitDisplay(TGTab *tabPtr); //InitDisplay
  AliHLTPHOSOnlineDisplay *fOnlineDisplayPtr; //fOnlineDisplayPtr
  AliHLTPHOSSharedMemoryInterface *fShmPtr;   //fShmPtr

  AliHLTPHOSFourier *fFourierPtr; //fFourierPtr
  char fGainText[256];  //fGainText

  unsigned long fEvtCnt; //fEvtCnt

  AliHLTPHOSOnlineDisplayFourierTab(const AliHLTPHOSOnlineDisplayFourierTab&); //copy constructor
  AliHLTPHOSOnlineDisplayFourierTab & operator=(const AliHLTPHOSOnlineDisplayFourierTab); //assignement operator


};


#endif
