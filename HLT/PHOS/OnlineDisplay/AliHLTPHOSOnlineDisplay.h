#ifndef ALIHLTPHOSONLINEDISPLAY
#define ALIHLTPHOSONLINEDISPLAY

/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#include "HOMERData.h"
#include "HOMERReader.h"
#include "HOMERWriter.h"
#include "Rtypes.h"
#include <TGFrame.h>
#include "TH2.h"
#include "AliHLTPHOSGetEventButton.h" 


class TCanvas;
class AliHLTPHOSRcuCellEnergyDataStruct;


class AliHLTPHOSOnlineDisplay : public  TGMainFrame
{
 public:
  AliHLTPHOSOnlineDisplay();
  ~AliHLTPHOSOnlineDisplay();
  static int GetNextEvent();
  static AliHLTPHOSOnlineDisplay* Instance();  

 private:
  static AliHLTPHOSGetEventButton* fgEventButtPtr; 
  static AliHLTPHOSOnlineDisplay* fgInstancePtr;
  static HOMERReader* fgHomerReaderPtr;
  static TH2S *legoPlotLGPtr;
  static TH2S *legoPlotHGPtr;
  static char *fgDefaultDet;        
  static char *fgDefaultDataType;   
  static int fgEvntCnt;
  static TCanvas *fgCanvasHGPtr;
  static TCanvas *fgCanvasLGPtr;
  static Bool_t fgAccumulate;
};


#endif
