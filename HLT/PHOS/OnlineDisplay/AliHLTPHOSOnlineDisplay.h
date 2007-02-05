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
  static HOMERReader* homerReaderPtr;
  static TH2S *legoPlotLGPtr;
  static TH2S *legoPlotHGPtr;
  static char *fgDefaultDet;         //= "SOHP"; //PHOS written backwards
  static char *fgDefaultDataType;   //= "RENELLEC"; //CELLENER written backwards  
  static TCanvas *fgCanvasHGPtr[100];
  static TCanvas *fgCanvasLGPtr[100];
  static int fgEvntCnt;
};


#endif
