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
#include <TCanvas.h>
#include "TGFrame.h"
//#include "AliHLTPHOSCommonDefs.h";
//#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSRcuChannelDataStruct.h"
#include "AliHLTPHOSOnlineDisplayEventTab.h"
#include "AliHLTPHOSOnlineDisplayCalibTab.h"
#include "AliHLTPHOSOnlineDisplayRawTab.h"
#include "AliHLTPHOSOnlineDisplayRawDataMenu.h"


//#include "AliHLTPHOSCommonDefs.h"

#include "AliHLTPHOSConstants.h"
using namespace PhosHLTConst;

#define N_SAMPLES 70

//#define MAX_HISTOGRAMS 25
#define MAX_HISTOGRAMS 320

class TCanvas;
class AliHLTPHOSRcuCellEnergyDataStruct;

class AliHLTPHOSOnlineDisplay : public  TGMainFrame
{
 public:
  ~AliHLTPHOSOnlineDisplay();
  int GetNextEvent();
  int GetNextEventRaw();
  int GetHistogram();
  void InitDisplay();
  void EvaluateAverage();
  int ScanArguments(int argc, char** argv);

  //  static void ShowRawData();

  void ShowRawData();

  static AliHLTPHOSOnlineDisplay* Instance(int argc, char** argv);  
 private:
  AliHLTPHOSOnlineDisplay();
  AliHLTPHOSOnlineDisplay(int argc, char** argv);
  static AliHLTPHOSOnlineDisplayRawDataMenu *fgRawMenuPtr;
  static AliHLTPHOSOnlineDisplayEventTab  *fgEventTabPtr;
  static AliHLTPHOSOnlineDisplayCalibTab  *fgCalibTabPtr;
  static AliHLTPHOSOnlineDisplayRawTab    *fgRawTabPtr;
  static TGTab               *fTab;
  static AliHLTPHOSOnlineDisplay* fgInstancePtr;
  static unsigned int fgNHosts;
  static unsigned int fgNPorts;
  static HOMERReader* fgHomerReaderPtr;
  static HOMERReader* fgHomerReadersPtr[MAX_HOSTS];
  static char  *fgHosts[MAX_HOSTS];
  static short unsigned    *fgPorts;
  static Bool_t fgAccumulate;
  static Bool_t fgSyncronize;



  TCanvas  *fgRawDataCanvas;
  TH1D     *fgRawDataPlotsPtr[MAX_HISTOGRAMS];

  //  static TCanvas  *fgRawDataCanvas;
  //  static TH1D     *fgRawDataPlotsPtr[MAX_HISTOGRAMS];

  //  static TH1D     *fgRawDataPlotsPtr[25];

};


#endif
