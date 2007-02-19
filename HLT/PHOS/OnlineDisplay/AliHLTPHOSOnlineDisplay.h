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
#include "TGTab.h"
#include <TRootEmbeddedCanvas.h>
#include "TGFrame.h"

class TCanvas;
class AliHLTPHOSRcuCellEnergyDataStruct;

class AliHLTPHOSOnlineDisplay : public  TGMainFrame
{
 public:
  AliHLTPHOSOnlineDisplay();
  AliHLTPHOSOnlineDisplay(char *hosname, int port);
  ~AliHLTPHOSOnlineDisplay();
  //  static int GetNextEvent();
  int GetNextEvent();
  void InitDisplay();
  void UpdateDisplay();
  static AliHLTPHOSOnlineDisplay* Instance(char *hostname, int port);  

 private:
  static TGCompositeFrame    *fFrame1, *fF1, *fF2, *fF3, *fF4, *fF5;
  static TGTab               *fTab;
  static TRootEmbeddedCanvas *fEc1, *fEc2, *fEc3, *fEc4, *fEc5, *fEc6;
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
  static Bool_t test[17920][2];
};


#endif
