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

#define MAX_HOSTS 10
#define MAX_HOSTNAME_LENGTH 64
#define DEFAULT_PORT 42001 
//#define MAX_PORTS_PER_HOST

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
  int GetNextEvent2();
  void InitDisplay();
  void UpdateDisplay();
  static int ScanArguments(int argc, char** argv);
  //  static AliHLTPHOSOnlineDisplay* Instance(char *hostname, int port);  
  static AliHLTPHOSOnlineDisplay* Instance();  
 private:
  static TGCompositeFrame    *fFrame1, *fF1, *fF2, *fF3, *fF4, *fF5, *fSubF1, *fSubF2, *fSubF3;
  static TGTab               *fTab;
  static TGTab               *fSubTab;
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

  static  char *host;
  //  char *port;
  static int port;

  static unsigned int fgNHosts;
  static unsigned int fgNPorts;
  //  static const char  **fgHosts;

  static HOMERReader* fgHomerReadersPtr[MAX_HOSTS];

  static char  *fgHosts[MAX_HOSTS];
  //  static char  **fgHosts;
 // static short unsigned    fgPorts[MAX_HOSTS];
  static short unsigned    *fgPorts;

  static Bool_t fgAccumulate;
  //  static Bool_t test[17920][2];
};


#endif
